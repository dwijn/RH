#!/usr/bin/env python3
"""Regression tests for RH: compare a candidate build against an unmodified
baseline build on small, fast runs.

usage: run_tests.py --baseline DIR --candidate DIR [options]

  --baseline DIR   the reference build (e.g. upstream RH built with its
                   makefiles). DIR is a source tree with the programs built
                   in place or in build/, or a CMake build directory.
  --candidate DIR  the build under test, given the same way. Example
                   inputs, Atoms/, Atmos/ and Molecules/ are taken from its
                   source tree.
  --asan DIR       also run the candidate cases with an AddressSanitizer
                   build, e.g. a CMake build directory configured with
                   -DRH_SANITIZE=address
  --exact          require byte-identical output between baseline and
                   candidate at 1 thread. Use for pure refactors / memory
                   management; leave off when summation order may change.
  --tol X          allowed relative difference otherwise. Defaults are set per
                   geometry from the noise of unmodified code against itself
                   at different thread counts: 1e-9 (1-D, 3-D) and 1e-7 (2-D,
                   whose shipped example is ill-conditioned, ~2e-9 noise)
  --threads "A B"  thread counts to check against the candidate's own
                   1-thread result (default "3 8": an uneven split of
                   the work and more threads than cores)
  --memtol X       candidate peak memory may be at most X times the
                   baseline's in the memory case (default 1.05)
  --only PATTERN   only run cases whose name contains PATTERN
  --skip PATTERN   skip cases whose name contains PATTERN
  --work DIR      scratch directory (default tests/work)

What is checked, per case (1-D, 2-D, 3-D and spherical geometries, with and without
Stokes, CRD only, two active atoms, ...):
  * candidate at 1 thread runs cleanly and matches the baseline
  * candidate at each thread count runs cleanly and matches its own
    1-thread result (up to floating-point summation order)
  * with --asan: no sanitizer reports at 1 thread and the largest count
    (3-D at the largest count only)
plus a memory-scaling case and a Kurucz line list with a leading blank line.

The code reads atomic data through hard-coded ../../Atoms paths, so runs
are made in WORK/runs/<name>, with WORK/Atoms etc. symlinked to the
candidate tree.

(Python 3.8 or newer.)
Standard library only. The tree logic is rhtest.tree, the comparison of
spectra rhtest.compare, the run directories rhtest.run.
"""
import sys

if sys.version_info < (3, 8):
    sys.exit("run_tests.py needs Python 3.8 or newer, this is %d.%d (%s).\n"
             "Run it with a newer python3."
             % (sys.version_info[0], sys.version_info[1], sys.executable))

import mmap
import os
import re
import shutil
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from rhtest import compare, keywords, run, tree, workdir      # noqa: E402

HELP = __doc__.split("\n\n(Python 3.8")[0] + "\n"
USAGE = "\n".join(HELP.split("\n")[:5]) + "\n"

GEOMDIR = tree.GEOMDIR


def fail(message):
    print(message)
    sys.exit(1)


def parse_args(argv):
    """The options, in order, as the old loop read them."""
    o = dict(baseline="", candidate="", asan="", exact=False, tol="", threads="3 8",
             memtol="1.05", only="", skip="", work=os.path.join(HERE, "work"))
    with_value = {"--baseline": "baseline", "--candidate": "candidate", "--asan": "asan",
                  "--tol": "tol", "--threads": "threads", "--memtol": "memtol",
                  "--only": "only", "--skip": "skip", "--work": "work"}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in with_value:
            if i + 1 >= len(argv):
                # the old loop never ended here: shift 2 failed and left the option
                fail("option needs a value: " + a)
            o[with_value[a]] = argv[i + 1]
            i += 2
        elif a == "--exact":
            o["exact"] = True
            i += 1
        elif a in ("-h", "--help"):
            sys.stdout.write(HELP)
            sys.exit(0)
        else:
            fail("unknown option: " + a)
    if not (o["baseline"] and o["candidate"]):
        sys.stdout.write(USAGE)
        sys.exit(1)
    for t in o["threads"].split():
        if not (t.isdigit() and int(t) > 0):
            fail("--threads needs positive integers, not '%s'" % t)
    for k in ("tol", "memtol"):
        if o[k]:
            try:
                float(o[k])
            except ValueError:
                fail("--%s needs a number, not '%s'" % (k, o[k]))
    return o


def absdir(path, what):
    """cd DIR && pwd"""
    if not os.path.isdir(path):
        fail("%s: not a directory: %s" % (what, path))
    return os.path.abspath(path)


class Suite:
    def __init__(self, o):
        self.o = o
        self.exact = o["exact"]
        self.tol = o["tol"]
        self.threads = o["threads"].split()
        self.npass = self.nfail = self.nskip = 0
        self.failed = []
        self.rundir = ""            # the run directory the last prepare() made

    # -- reporting ------------------------------------------------------
    def record(self, status, name, detail):
        if status == "PASS":
            self.npass += 1
        elif status == "FAIL":
            self.nfail += 1
            self.failed.append(name)
        else:
            self.nskip += 1
        print("[%s] %-34s %s" % (status, name, detail))
        # Keep the log of a failed run (the run directory is reused by later
        # cases), so that intermittent failures can be looked at afterwards
        log = os.path.join(self.rundir, "log.log")
        if status == "FAIL" and self.rundir and os.path.isfile(log):
            keep = os.path.join(self.work, "failed", "%s_%s" % (
                name.replace(" ", "_").replace("/", "_").replace("=", "_"),
                time.strftime("%H%M%S")))
            os.makedirs(keep, exist_ok=True)
            for f in ("log.log", "keyword.input", "time.txt"):
                try:
                    shutil.copy(os.path.join(self.rundir, f), keep)
                except OSError:
                    pass
            print("       log kept in %s" % keep)

    # -- run directories --------------------------------------------------
    def prepare(self, geom, name, label, n, mods):
        """Run directory WORK/runs/<label>_<geom>_<name>_n<n> from the
        candidate's example: keyword.input with N_THREADS, N_MAX_ITER = 4 and
        the modifications mods (KEY=VALUE or @H6_ACTIVE) set."""
        d = self.rundir = run.fresh_dir(self.wd, "%s_%s_%s_n%s" % (label, geom, name, n))
        t = run.template(self.cand_src, GEOMDIR[geom])
        for f in sorted(os.listdir(t)):
            if not f.startswith(".") and os.path.isfile(os.path.join(t, f)):
                shutil.copy(os.path.join(t, f), d)
        pairs = [("N_THREADS", n), ("N_MAX_ITER", "4")]
        for m in mods:
            if m == "@H6_ACTIVE":
                self.h6_active(os.path.join(d, "atoms.input"))
            elif "=" in m[1:]:
                pairs.append(tuple(m.split("=", 1)))
        keywords.edit_file(os.path.join(d, "keyword.input"), pairs)
        return d

    @staticmethod
    def h6_active(path):
        with open(path, newline="", errors="surrogateescape") as f:
            lines = f.read().split("\n")
        lines = [re.sub(r"(H_6\.atom[ \t\r\f\v]+)PASSIVE", r"\1ACTIVE", l, count=1) for l in lines]
        with open(path, "w", newline="", errors="surrogateescape") as f:
            f.write("\n".join(lines))

    def execute(self, bin_dir, geom):
        """Run the program of the tree in the last run directory; the exit
        status of the program (128 + signal). Sanitizers exit with 99 when
        they report an error (their default, 1, would look like "not
        converged")."""
        env = dict(os.environ, ASAN_OPTIONS="detect_leaks=0:exitcode=99",
                   UBSAN_OPTIONS="exitcode=99")
        prog = os.path.join(bin_dir, GEOMDIR[geom], GEOMDIR[geom])
        return run.execute_to_log(prog, self.rundir, env=env)

    def problem(self):
        return run.log_problem(os.path.join(self.rundir, "log.log"))

    def compare(self, a, b, geom, exact=False):
        """(ok, text) of the spectra of two run directories."""
        tol = self.tol or ("1e-7" if geom == "2d" else "1e-9")
        status, text = compare.compare_spectra(
            os.path.join(a, "spectrum.out"), os.path.join(b, "spectrum.out"),
            tol=float(tol), exact=exact,
            stokes=True if keywords.has_stokes_input(os.path.join(a, "keyword.input")) else None)
        return status == 0, text

    # -- the cases ---------------------------------------------------------
    def run_case(self, geom, name, mods):
        cid = "%s/%s" % (geom, name)
        threads_max = max(self.threads, key=int) if self.threads else ""
        self.prepare(geom, name, "base", "1", mods)
        B1 = self.rundir
        rc = self.execute(self.base_bin, geom)
        if not run.completed(rc):
            self.record("SKIP", cid, "baseline itself failed (rc=%d), nothing to compare against" % rc)
            return

        self.prepare(geom, name, "cand", "1", mods)
        C1 = self.rundir
        rc = self.execute(self.cand_bin, geom)
        why = self.problem()
        if not run.completed(rc) or why:
            self.record("FAIL", "%s candidate N=1" % cid, "rc=%d %s" % (rc, why))
            return
        ok, text = self.compare(B1, C1, geom, exact=self.exact)
        self.record("PASS" if ok else "FAIL", "%s candidate N=1 vs baseline" % cid, text)

        for n in self.threads:
            self.prepare(geom, name, "cand", n, mods)
            CN = self.rundir
            rc = self.execute(self.cand_bin, geom)
            why = self.problem()
            if not run.completed(rc) or why:
                self.record("FAIL", "%s candidate N=%s" % (cid, n), "rc=%d %s" % (rc, why))
                continue
            ok, text = self.compare(C1, CN, geom)
            self.record("PASS" if ok else "FAIL", "%s N=%s vs N=1" % (cid, n), text)

        if not self.asan_bin:
            return
        # 3-D only at the largest count: the code that runs only with one
        # thread (Gamma and the rates accumulated in place, the background
        # without workers) is shared by all geometries and checked by the
        # cheaper cases, and 3-D under the sanitizers at 1 thread took 105 s
        # of CI's slowest job
        counts = ["1"] + ([threads_max] if threads_max else [])
        if geom == "3d":
            counts = [threads_max] if threads_max else []
        for n in counts:
            self.prepare(geom, name, "asan", n, mods)
            rc = self.execute(self.asan_bin, geom)
            why = self.problem()
            if run.completed(rc) and not why:
                self.record("PASS", "%s ASan N=%s" % (cid, n), "clean")
            else:
                self.record("FAIL", "%s ASan N=%s" % (cid, n), "rc=%d %s" % (rc, why))

    @staticmethod
    def mapped_mb(bin_dir, d):
        """MB of the background file that the program in bin dir maps
        read-only (P-18, recognised by its warning text) in run dir d: the
        pages it reads count in the peak RSS but are clean, shared page
        cache, which pread used to copy from, not memory the program
        allocates. 0 for a program that reads the file with pread."""
        prog = os.path.join(bin_dir, GEOMDIR["2d"], GEOMDIR["2d"])
        bg = os.path.join(d, "background.dat")
        with open(prog, "rb") as f:
            maps = False
            if os.fstat(f.fileno()).st_size > 0:
                with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as m:
                    maps = m.find(b"Cannot map the background file") >= 0
        if maps and os.path.isfile(bg):
            return int("%.0f" % (os.path.getsize(bg) / 1048576))
        return 0

    def rss_mb(self):
        """Peak RSS of the last run in MB; OSError (a failure of the case)
        where it was not measured."""
        mb = run.peak_rss_mb(self.rundir)
        if mb is None:
            raise OSError("no peak memory in %s/time.txt" % self.rundir)
        return mb

    def run_mem_case(self):
        """Peak memory must not grow with the number of wavelengths on a
        large grid. (Per-wavelength caching of opacity buffers once added
        Nspect*Nspace*40 bytes.) Both runs use SOLVER = NG with the Ng
        settings the shipped 2-D keyword file had before S-11 made ANDERSON
        the 2-D default, so that a change of the default solver does not
        count as growth and Ng's history is still measured. 96 columns:
        about 115 MB, of which 8 MB do not grow with the grid; the old
        caching would add about 85% to the part that does, at any grid size.
        A program that maps the background file (P-18) has its pages in its
        RSS too (26 MB here); they are page cache, not allocations, and are
        not counted for either program (mapped_mb)."""
        if not run.have_time():
            self.record("SKIP", "2d/mem", "no /usr/bin/time")
            return
        tiled = os.path.join(self.work, "tiled_2d.atmos")
        r = subprocess.run([sys.executable, "-I", os.path.join(HERE, "make_tiled_2d_atmos.py"),
                            os.path.join(self.cand_src, "Atmos", "FALC_03x82.atmos"), tiled, "96"],
                           stdout=subprocess.DEVNULL)
        if r.returncode != 0:
            self.record("FAIL", "2d/mem", "tiling failed")
            return
        mods = ["ATMOS_FILE=" + tiled, "N_MAX_ITER=1", "SOLVER=NG", "NG_DELAY=6", "NG_ORDER=2",
                "NG_PERIOD=3", "PRD_ITER_LIMIT=1.0E-2"]
        self.prepare("2d", "mem", "base", "1", mods)
        B = self.rundir
        rc = self.execute(self.base_bin, "2d")
        if not run.completed(rc):
            self.record("SKIP", "2d/mem", "baseline failed (rc=%d)" % rc)
            return
        bm_raw = self.rss_mb()
        self.prepare("2d", "mem", "cand", "1", mods)
        C = self.rundir
        rc = self.execute(self.cand_bin, "2d")
        why = self.problem()
        if not run.completed(rc) or why:
            self.record("FAIL", "2d/mem candidate", "rc=%d %s" % (rc, why))
            return
        # Peak RSS without the mapped background file of either program (P-18)
        bmap = self.mapped_mb(self.base_bin, B)
        cmap = self.mapped_mb(self.cand_bin, C)
        bm = bm_raw - bmap
        cm = self.rss_mb() - cmap
        note = ("limit x%s; mapped background file not counted: baseline %d MB, candidate %d MB"
                % (self.o["memtol"], bmap, cmap))
        text = "candidate %d MB vs baseline %d MB (%s)" % (cm, bm, note)
        self.record("PASS" if cm <= bm * float(self.o["memtol"]) else "FAIL", "2d/mem peak memory", text)

    def run_kurucz_blank_case(self):
        """A Kurucz list with a leading blank line (as downloaded from the
        current kurucz.cfa.harvard.edu) must give the same result as the file
        without it."""
        mods = ["STOKES_MODE=NO_STOKES", "STOKES_INPUT=none"]
        A = self.prepare("1d", "kurucz", "cand", "1", mods)
        with open(os.path.join(A, "kurucz.input"), "w") as f:
            f.write("../../Atoms/Kurucz/gf0400.10\n")
        rc = self.execute(self.cand_bin, "1d")
        why = self.problem()
        if not run.completed(rc) or why:
            self.record("FAIL", "1d/kurucz_blank reference", "rc=%d %s" % (rc, why))
            return
        B = self.prepare("1d", "kurucz_blank", "cand", "1", mods)
        with open(os.path.join(self.cand_src, "Atoms", "Kurucz", "gf0400.10"), "rb") as f:
            data = f.read()
        with open(os.path.join(B, "kurucz_blank.dat"), "wb") as f:
            f.write(b"\n" + data)
        with open(os.path.join(B, "kurucz.input"), "w") as f:
            f.write("./kurucz_blank.dat\n")
        rc = self.execute(self.cand_bin, "1d")
        why = self.problem()
        if not run.completed(rc) or why:
            self.record("FAIL", "1d/kurucz_blank", "rc=%d %s" % (rc, why))
            return
        ok, text = self.compare(A, B, "1d", exact=True)
        if ok:
            self.record("PASS", "1d/kurucz_blank", "same as without the blank line")
        else:
            self.record("FAIL", "1d/kurucz_blank", text)

    def guarded(self, name, case, *args):
        try:
            case(*args)
        except OSError as e:
            self.record("FAIL", name, "error: %s" % e)

    def selected(self, name):
        return (not self.o["only"] or self.o["only"] in name) and \
               (not self.o["skip"] or self.o["skip"] not in name)

    def cases(self):
        w = self.work
        return [
            ("1d", "nostokes", ["STOKES_MODE=NO_STOKES", "STOKES_INPUT=none"]),
            ("1d", "fieldfree", ["STOKES_MODE=FIELD_FREE"]),
            ("1d", "fullstokes", ["STOKES_MODE=FULL_STOKES"]),
            ("1d", "crd", ["STOKES_MODE=NO_STOKES", "STOKES_INPUT=none", "PRD_N_MAX_ITER=0"]),
            ("1d", "twoactive", ["STOKES_MODE=NO_STOKES", "STOKES_INPUT=none", "@H6_ACTIVE"]),
            ("2d", "small", []),
            ("3d", "small", ["N_MAX_ITER=2"]),
            ("2d", "fieldfree", ["STOKES_INPUT=%s/field_2d.B" % w, "STOKES_MODE=FIELD_FREE"]),
            ("3d", "fieldfree", ["STOKES_INPUT=%s/field_3d.B" % w, "STOKES_MODE=FIELD_FREE",
                                 "N_MAX_ITER=2"]),
            ("2d", "lcbezier", ["STOKES_INPUT=%s/field_2d.B" % w, "STOKES_MODE=FIELD_FREE",
                                "S_INTERPOLATION_LONGCHAR=S_BEZIER3"]),
            ("2d", "fullstokes", ["STOKES_INPUT=%s/field_2d.B" % w, "STOKES_MODE=FULL_STOKES"]),
            ("sph", "small", []),
        ]

    def main(self):
        o = self.o
        base = absdir(o["baseline"], "--baseline")
        cand = absdir(o["candidate"], "--candidate")
        asan = absdir(o["asan"], "--asan") if o["asan"] else ""
        self.base_bin, self.cand_bin = tree.tree_bin(base), tree.tree_bin(cand)
        self.asan_bin = ""
        if asan:
            self.asan_bin = tree.tree_bin(asan)
            if not tree.tree_sanitized(asan) and os.path.isfile(
                    os.path.join(self.asan_bin, "CMakeCache.txt")):
                fail("--asan %s: %s is not configured with -DRH_SANITIZE" % (asan, self.asan_bin))
        self.cand_src = tree.tree_src(cand)
        self.wd = workdir.WorkDir(self.cand_src, o["work"], relink=True, announce=False)
        self.work = self.wd.path
        # Magnetic fields for the 2-D and 3-D Stokes cases. make_field.py is
        # the one of this checkout (as before), not the candidate's tree, which
        # may be an older one.
        try:
            for g, atm in (("2d", "FALC_03x82.atmos"), ("3d", "FALC_04x04x82.atmos")):
                run.make_field(os.path.dirname(HERE), g, os.path.join(self.cand_src, "Atmos", atm),
                               os.path.join(self.work, "field_%s.B" % g))
        except (OSError, subprocess.CalledProcessError):
            print("make_field.py failed: the 2d/3d fieldfree cases will fail")

        for t in (self.base_bin, self.cand_bin) + ((self.asan_bin,) if asan else ()):
            miss = tree.missing_programs(t)
            if miss:
                fail("missing binary: " + miss[0])

        if not run.have_time():
            print("warning: /usr/bin/time missing, memory case skipped")

        print("baseline : " + self.base_bin)
        print("candidate: " + self.cand_bin + ("   (ASan: %s)" % self.asan_bin if asan else ""))
        print("mode     : %s, threads: %s" % (
            "byte-identical" if self.exact else "tolerance " + (self.tol or "per geometry"),
            o["threads"]))
        print()

        # An OSError of a case (a directory that cannot be made, a program that
        # cannot be started, a file that cannot be read) fails that case; the
        # others run and the summary is printed
        for geom, name, mods in self.cases():
            if self.selected("%s/%s" % (geom, name)):
                self.guarded("%s/%s" % (geom, name), self.run_case, geom, name, mods)
        if self.selected("2d/mem"):
            self.guarded("2d/mem", self.run_mem_case)
        if self.selected("1d/kurucz_blank"):
            self.guarded("1d/kurucz_blank", self.run_kurucz_blank_case)

        print()
        print("passed %d, failed %d, skipped %d" % (self.npass, self.nfail, self.nskip))
        if self.nfail:
            for f in self.failed:
                print("FAILED: " + f)
            sys.exit(1)


def main(argv):
    sys.stdout.reconfigure(line_buffering=True)
    Suite(parse_args(argv)).main()


if __name__ == "__main__":
    main(sys.argv[1:])
