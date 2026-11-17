#!/usr/bin/env python3
"""Stress, convergence and performance cases for RH: NLTE problems run to
convergence (or to a fixed number of iterations for the timing cases) with
one or two builds, reporting iterations, final population change, wall time,
peak memory, and the difference between the builds.

usage: run_stress.py --tree NAME=DIR [--tree NAME=DIR ...] [options]

  --tree NAME=DIR  a build of RH: a source tree with rhf1d, rhsc2d, rhsc3d
                   and rhsphere built in place or in build/, or a CMake
                   build directory. The first tree is the reference for
                   comparisons.
  --only PATTERN   only run cases whose name matches PATTERN (grep -E).
                   Matching cases that --set or a missing --bifrost exclude
                   are listed as "not run" with the reason; when no case is
                   left the script stops with status 1 and says why
  --set conv|perf|all   which group of cases to run (default conv)
  --jobs N         run up to N cases at the same time (default 4). Timing
                   cases always run one at a time.
  --threads N      override N_THREADS of the multi-D cases
  --work DIR       scratch directory (default tests/work_stress)
  --kw 'K=V;...'   keyword.input settings applied to every case after its
                   own (e.g. the reference recipe of make_references.sh);
                   -K removes keyword K. The G solvers refuse keywords
                   they do not honour (docs/SOLVER.md, Anderson as the
                   default solver, item 2): e.g. -NG_DELAY;-NG_PERIOD for
                   the shipped 1-D file; a case that sets one, such as
                   1d_falxco_collrel (COLL_RELAX_FACTOR), stops at input
  --ref DIR        reference solutions (make_references.sh): every case with
                   a DIR/CASE directory is compared with it (population and
                   spectrum errors, analyze_run.py compare)
  --bifrost DIR    directory with the FITS cubes of a Bifrost snapshot
                   (en024048_hion, Carlsson et al. 2016, A&A 585, A4; not
                   distributed with RH, download it yourself); also taken
                   from $RH_BIFROST_DIR. Without it the bifrost cases are
                   skipped (an error if --only leaves nothing else)
  --list           print the cases and exit

Results go to WORK/results_<date-time>.tsv (one file per invocation, one line
per case and tree) and are printed as a table. Every run writes its

(Python 3.8 or newer. The tree logic is rhtest.tree, the runs rhtest.run; the analysis is
analyze_run.py, which needs numpy.)
"""
import sys

if sys.version_info < (3, 8):
    sys.exit("run_stress.py needs Python 3.8 or newer, this is %d.%d (%s).\n"
             "Run it with a newer python3."
             % (sys.version_info[0], sys.version_info[1], sys.executable))

import json
import os
import re
import shutil
import subprocess
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from rhtest import compare, keywords, run, tree, workdir      # noqa: E402

DOC = __doc__.split("\n\n(Python 3.8")[0] + "\n"
GEOMDIR = {"1d": "rhf1d", "2d": "rhsc2d", "3d": "rhsc3d", "sph": "rhsphere"}

# Cases: NAME | SET | GEOM | ATMOS | ACTIVE ATOMS | KEYWORDS
#
# ATMOS is a file in Atmos/ or tests/atmos/, gen:2d:NX / gen:3d:NX:NY for
# a perturbed FALC (make_perturbed_atmos.py), or a cut of the Bifrost snapshot
# (make_bifrost_atmos.py, z from -0.5 to 3 Mm at every second point):
# bifrost:1d:X (column X of row 0), bifrost:2d:S (row 0, every S-th column),
# bifrost:3d:S (the whole box, every S-th point; S divides 504, so the cut is
# periodic like the box), bifrost:3d:1:N (an N x N window at full resolution,
# not periodic: for timing only). The keyword STOKES_INPUT=@B uses the
# field of the cut. ACTIVE ATOMS is a comma-separated
# list of atom files in Atoms/ to make active; each replaces the line of the
# same element in the geometry's example atoms.input. KEYWORDS are
# ';'-separated keyword.input settings (K=V, or -K to remove K), plus
# @NOMOL=NAME to drop a molecule
# (e.g. MgH, which needs neutral Mg, when Mg II is the Mg model).
#
# 1-D cases run on 1 thread unless they set N_THREADS, so that two builds can
# be compared without summation-order noise (some cases, e.g. He_large, end at
# visibly different "converged" solutions for different thread counts).
#
# Cases run each program's default solver: ANDERSON everywhere (S-10; S-11
# in 2-D and 3-D, where it extrapolates POPULATIONS_RHO). Those that test
# the Ng iteration itself (no Ng, PRD Ng, collisional relaxation, tight
# convergence without Ng) pin it with $NG: SOLVER = NG and the Ng settings
# the shipped keyword files had, written out, since S-10 and S-11 removed
# them from those files. They
# keep testing Ng when the default changes. --kw settings come after a
# case's own; a G solver set there stops these cases at input, because they
# set keywords it does not honour.
#
# conv: run to convergence (N_MAX_ITER large), compare converged results.
# perf: fixed number of iterations, for wall time and peak memory.
# The comment after each group says what the case is known to expose; where
# it names Ng, that is what the case showed with NG, the default before
# S-10 (--kw with the $NG settings shows it again; MODERNIZATION.md 2.2 and
# 2.5 have the results of both solvers).
L = "ITER_LIMIT=1e-3;N_MAX_ITER=200"
NG = "SOLVER=NG;NG_ORDER=2;NG_DELAY=10;NG_PERIOD=3;PRD_ITER_LIMIT=1.0E-2"
CASES = [
    # Well-behaved references (converge in 15-30 iterations)
    f"1d_falc_h6_ca_mg  |conv|1d|FALC_82.atmos  |H_6,CaII,MgII     |{L}",
    f"1d_fala_h_ca      |conv|1d|FALA_80.atmos  |H,CaII            |{L}",
    f"1d_falc_h20       |conv|1d|FALC_82.atmos  |H_20              |{L}",
    f"1d_falc_mgi66     |conv|1d|FALC_82.atmos  |MgI_66            |{L}",
    f"1d_falc_fe33      |conv|1d|FALC_82.atmos  |Fe33              |{L}",
    f"1d_falc_nong      |conv|1d|FALC_82.atmos  |H_6,CaII,MgII     |{L};{NG};NG_ORDER=0",
    f"1d_falxco_h6crd   |conv|1d|FALXCO_80.atmos|H_6_CRD,CaII      |{L}",
    # Slow / non-convergence: PRD sub-iteration (H Lyman PRD), Ng without safeguards
    f"1d_falxco_h6_ca   |conv|1d|FALXCO_80.atmos|H_6,CaII          |{L}",
    f"1d_t5750_h6       |conv|1d|T5750_g4.5_p00_ext.atmos|H_6     |{L}",
    f"1d_falc_helarge   |conv|1d|FALC_82.atmos  |He_large          |{L}",
    f"1d_falc_h6_tight  |conv|1d|FALC_82.atmos  |H_6               |ITER_LIMIT=1e-8;N_MAX_ITER=300;{NG};NG_ORDER=0",
    f"1d_falc_civi      |conv|1d|FALC_82.atmos  |CI-VI_Judge       |{L}",
    # Moving atmospheres: linear vs Bezier (they differ by 5-16 % in the
    # chromosphere on this 82-point grid: depth resolution, not the scheme)
    f"1d_v10_ca_lin     |conv|1d|FALC_v10.atmos |CaII              |{L};S_INTERPOLATION=S_LINEAR",
    f"1d_v10_ca_bez     |conv|1d|FALC_v10.atmos |CaII              |{L};S_INTERPOLATION=S_BEZIER3",
    f"1d_shock_h6       |conv|1d|FALC_shock.atmos|H_6              |{L}",
    # Crashes / wrong answers from options
    f"1d_falc_cail      |conv|1d|FALC_82.atmos  |CaII_large        |{L}",
    f"1d_falc_solvene   |conv|1d|FALC_82.atmos  |H_6,CaII          |{L};SOLVE_NE=ITERATION",
    f"1d_shock_h6crd_ne |conv|1d|FALC_shock.atmos|H_6_CRD          |{L};SOLVE_NE=ITERATION",
    f"1d_falc_prdng     |conv|1d|FALC_82.atmos  |CaII              |{L};{NG};PRD_NG_ORDER=2",
    f"1d_falxco_collrel |conv|1d|FALXCO_80.atmos|H_6               |{L};{NG};COLL_RELAX_NSTEP=20;COLL_RELAX_FACTOR=0.01",
    # Multi-D
    f"2d_pert_h6_ca_mg  |conv|2d|gen:2d:32      |H_6,CaII,MgII     |{L};N_THREADS=8;S_INTERPOLATION=S_BEZIER3;@NOMOL=MgH",
    f"3d_pert_h6_ca     |conv|3d|gen:3d:8:8     |H_6,CaII          |{L};N_THREADS=8",
    f"2d_perf_many      |perf|2d|gen:2d:128     |H_6,CaII,MgII,Fe33|N_MAX_ITER=5;N_THREADS=16;S_INTERPOLATION=S_BEZIER3;@NOMOL=MgH",
    f"3d_perf_many      |perf|3d|gen:3d:24:24   |H_6,CaII,MgII,Fe33|N_MAX_ITER=3;N_THREADS=16;@NOMOL=MgH",
    f"3d_perf_mgi66     |perf|3d|gen:3d:16:16   |MgI_66            |N_MAX_ITER=3;N_THREADS=16",
    # Bifrost cuts (need --bifrost). Ca II in CRD: with PRD the iteration
    # diverged in these moving atmospheres (B-29, fixed); the _prd cases test that
    f"1d_bf_c260_ca     |conv|1d|bifrost:1d:260 |CaII_CRD          |{L}",
    f"1d_bf_c260_prd    |conv|1d|bifrost:1d:260 |CaII              |{L}",
    f"1d_bf_c100_prd    |conv|1d|bifrost:1d:100 |CaII              |{L}",
    f"2d_bf_s4_ca       |conv|2d|bifrost:2d:4   |CaII_CRD          |{L};N_THREADS=16",
    f"2d_bf_s4_prd      |conv|2d|bifrost:2d:4   |CaII              |{L};N_THREADS=16",
    f"3d_bf_s21_ca      |conv|3d|bifrost:3d:21  |CaII_CRD          |{L};N_THREADS=16",
    f"2d_bf_s1_perf     |perf|2d|bifrost:2d:1   |CaII_CRD,H_6      |N_MAX_ITER=5;N_THREADS=16",
    f"3d_bf_w32_perf    |perf|3d|bifrost:3d:1:32|CaII_CRD          |N_MAX_ITER=3;N_THREADS=16",
]


def trim(s):
    return s.strip(" \t\n\r\f\v")


def fields(case):
    """name, set, geometry, atmosphere, atoms, keywords of a case line (not
    trimmed); keywords keep any further bars."""
    parts = case.split("|", 5)
    return parts + [""] * (6 - len(parts))


def fail(message, stream=sys.stdout):
    print(message, file=stream)
    sys.exit(1)


def parse_args(argv):
    o = dict(trees=[], only="", set="conv", jobs="4", threads="", work=os.path.join(HERE, "work_stress"),
             list=False, kw="", ref="", bifrost=os.environ.get("RH_BIFROST_DIR", ""))
    with_value = {"--tree": "tree", "--only": "only", "--set": "set", "--jobs": "jobs",
                  "--threads": "threads", "--work": "work", "--kw": "kw", "--ref": "ref",
                  "--bifrost": "bifrost"}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in with_value:
            if i + 1 >= len(argv):
                # the old loop never ended here: shift 2 failed and left the option
                fail("option needs a value: " + a)
            v = argv[i + 1]
            if a == "--tree":
                o["trees"].append(v)
            elif a == "--ref":
                if not os.path.isdir(v):
                    fail("--ref: not a directory: " + v, sys.stderr)
                o["ref"] = os.path.abspath(v)
            else:
                o[with_value[a]] = v
            i += 2
        elif a == "--list":
            o["list"] = True
            i += 1
        elif a in ("-h", "--help"):
            sys.stdout.write(DOC)
            sys.exit(0)
        else:
            fail("unknown option: " + a)
    return o


def grep_matches(pattern, text):
    """grep -Eq -e PATTERN on text: the patterns of --only are grep -E's."""
    r = subprocess.run(["grep", "-Eq", "-e", pattern], input=text + "\n",
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       universal_newlines=True)
    return r.returncode == 0


def list_cases():
    for c in CASES:
        n, s, g, a, at, kw = fields(c)
        print("%-18s %-4s %-3s %-15s %-28s %s" % (trim(n), trim(s), trim(g), trim(a), trim(at),
                                                  trim(kw)))


def element(path):
    """First word (upper case) of the first non-comment line of an atom file."""
    with open(path, errors="replace") as f:
        for l in f:
            if l.strip() and not l.lstrip().startswith("#"):
                return l.split()[0].upper()


def set_atoms(path, atomdir, wanted):
    """All atoms passive, then the requested ones active, each replacing the
    entry of the same element. SystemExit with a message if there is none."""
    with open(path) as f:
        lines = f.read().splitlines()
    for i, l in enumerate(lines):
        lines[i] = l.replace(" ACTIVE ", " PASSIVE ").replace("\tACTIVE\t", "\tPASSIVE\t")
    for a in filter(None, wanted.split(",")):
        el = element("%s/%s.atom" % (atomdir, a))
        for i, l in enumerate(lines):
            w = l.split()
            if w and "Atoms/" in w[0] and element("%s/%s" % (atomdir, os.path.basename(w[0]))) == el:
                lines[i] = "  ../../Atoms/%s.atom   ACTIVE   ZERO_RADIATION" % a
                break
        else:
            sys.exit("no %s entry in atoms.input to replace with %s" % (el, a))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def drop_molecule(path, name):
    """The molecule name (the line containing /<name>.molecule) removed from
    molecules.input and the molecule count, the first non-comment line,
    decremented."""
    with open(path) as f:
        lines = f.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    lines = [l for l in lines if "/%s.molecule" % name not in l]
    out, done = [], False
    for l in lines:
        if done or l.lstrip(" \t\r\f\v").startswith("#") or not l.split():
            out.append(l)
        else:
            m = re.match(r"[-+]?\d+(\.\d*)?", l.split()[0])
            out.append("%g" % (float(m.group(0)) - 1 if m else -1))
            done = True
    with open(path, "w") as f:
        f.write("".join(l + "\n" for l in out))


class Stress:
    def __init__(self, o):
        self.o = o
        self.lock = threading.Lock()

    # -- atmospheres -------------------------------------------------------
    def atmos_path(self, spec, setname):
        """The path of the atmosphere file of spec (generated if needed) for
        a case of the set setname. As in the old run_stress.sh, the name of a
        generated 3-D atmosphere has the set of the case where the size
        should be: `local s=... f=...${s/:/x}...` in the bash used the s of
        the case loop. The cases of one set share one file, whichever is
        generated first (3d_perf_mgi66 runs on the 24 x 24 columns of
        3d_perf_many)."""
        w, ref = self.work, self.ref_src
        if spec.startswith("gen:2d:"):
            n = spec[len("gen:2d:"):]
            f = "%s/pert_2d_%s.atmos" % (w, n)
            if not os.path.isfile(f):
                subprocess.run([sys.executable, "-I", os.path.join(HERE, "make_perturbed_atmos.py"),
                                ref + "/Atmos/FALC_03x82.atmos", f, n], stdout=subprocess.DEVNULL)
            return f
        if spec.startswith("gen:3d:"):
            s = spec[len("gen:3d:"):]
            f = "%s/pert_3d_%s.atmos" % (w, setname)
            if not os.path.isfile(f):
                subprocess.run([sys.executable, "-I", os.path.join(HERE, "make_perturbed_atmos.py"),
                                ref + "/Atmos/FALC_03x82.atmos", f, s.split(":")[0],
                                s.split(":")[-1]], stdout=subprocess.DEVNULL)
            return f
        if spec.startswith("bifrost:"):
            body = spec[len("bifrost:"):]
            g = body.split(":")[0]
            s = body.split(":", 1)[1] if ":" in body else body
            f = "%s/bifrost_%s_%s" % (w, g, s.replace(":", "x", 1))
            opt = []
            if g == "1d":
                opt = ["--1d", "--x0", s]
            elif g == "2d":
                opt = ["--stride", s]
            elif g == "3d" and ":" in s:
                opt = ["--3d", "--stride", s.split(":")[0], "--nx", s.split(":")[-1],
                       "--ny", s.split(":")[-1]]
            elif g == "3d":
                opt = ["--3d", "--stride", s]
            if not os.path.isfile(f + ".atmos"):
                subprocess.run([sys.executable, "-I", os.path.join(HERE, "make_bifrost_atmos.py"),
                                self.o["bifrost"], f] + opt + ["--zstride", "2"],
                               stdout=subprocess.DEVNULL)
            return f + ".atmos"
        if os.path.isfile(os.path.join(HERE, "atmos", spec)):
            return os.path.join(HERE, "atmos", spec)
        return "../../Atmos/" + spec

    # -- one run --------------------------------------------------------------
    def prepare(self, name, geom, atm, atoms, kws, tname, setname):
        """The run directory of a case and tree; the error text if the atoms
        cannot be set (also left in prepare_error)."""
        d = os.path.join(self.work, "runs", "%s__%s" % (name, tname))
        if os.path.lexists(d):
            shutil.rmtree(d)
        os.makedirs(d)
        t = run.template(self.ref_src, GEOMDIR[geom])
        for f in sorted(os.listdir(t)):
            if not f.startswith(".") and os.path.isfile(os.path.join(t, f)):
                shutil.copy(os.path.join(t, f), d)
        kwf = os.path.join(d, "keyword.input")
        with open(kwf) as f:
            text = f.read()
        # Plain intensity output unless a case asks otherwise
        text = keywords.apply(text, [("ATMOS_FILE", self.atmos_path(atm, setname)),
                                     ("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
                                     ("PRINT_CPU", "FALSE")])
        if geom == "1d":
            text = keywords.setkw(text, "N_THREADS", "1")
        text = keywords.apply(text, [("RADRATE_OUTPUT", "radrate.out"),
                                     ("COLLRATE_OUTPUT", "collrate.out")])
        allkw = kws + (";" + self.o["kw"] if self.o["kw"] else "")
        for m in filter(None, allkw.split(";")):
            if m.startswith("@NOMOL="):
                drop_molecule(os.path.join(d, "molecules.input"), m[len("@NOMOL="):])
            elif m == "STOKES_INPUT=@B":
                b = self.atmos_path(atm, setname)
                text = keywords.setkw(text, "STOKES_INPUT", (b[:-6] if b.endswith(".atmos") else b) + ".B")
            elif m.startswith("-") and len(m) > 1:
                text = keywords.unsetkw(text, m[1:])
            elif "=" in m[1:]:
                text = keywords.setkw(text, m.split("=", 1)[0], m.split("=", 1)[1])
        if self.o["threads"] and geom != "1d":
            text = keywords.setkw(text, "N_THREADS", self.o["threads"])
        with open(kwf, "w") as f:
            f.write(text)
        try:
            set_atoms(os.path.join(d, "atoms.input"), self.ref_src + "/Atoms", atoms)
        except SystemExit as e:
            msg = str(e.code)
            with open(os.path.join(d, "prepare_error"), "w") as f:
                f.write(msg + "\n")
            return d, msg
        return d, None

    def iter_limit(self, d):
        """ITER_LIMIT of the run's keyword.input, as the awk of the old script
        found it (every line that starts with it; the text between the first
        two = signs, spaces removed)."""
        vals = []
        with open(os.path.join(d, "keyword.input"), errors="replace") as f:
            for l in f.read().split("\n"):
                if re.match(r"[ \t\r\f\v]*ITER_LIMIT", l):
                    p = l.split("=")
                    vals.append(p[1].replace(" ", "") if len(p) > 1 else "")
        return "".join(v + "\n" for v in vals).rstrip("\n")

    def run_one(self, case, i):
        n, s, g, a, at, kw = [trim(x) for x in fields(case)]
        tname, tdir = self.names[i], self.dirs[i]
        d, err = self.prepare(n, g, a, at, kw, tname, s)
        if err is not None:
            self.append(n, tname, '{"status": "%s"}' % err)
            return
        prog = os.path.join(tdir, GEOMDIR[g], GEOMDIR[g])
        run.execute_to_log(prog, d)
        lim = self.iter_limit(d) or "1e-2"
        analyze = os.path.join(HERE, "analyze_run.py")
        js = subprocess.run([sys.executable, "-I", analyze, "summary", d, lim],
                            stdout=subprocess.PIPE, universal_newlines=True).stdout.rstrip("\n")
        ref = self.o["ref"]
        if ref and os.path.isdir(os.path.join(ref, n)):
            cmp_js = subprocess.run([sys.executable, "-I", analyze, "compare", d,
                                     os.path.join(ref, n)],
                                    stdout=subprocess.PIPE, universal_newlines=True).stdout
            try:
                a_ = json.loads(js)
                a_["ref"] = json.loads(cmp_js)
                js = json.dumps(a_)
            except ValueError:
                js = ""
        self.append(n, tname, js)

    def append(self, n, tname, js):
        with self.lock:
            with open(self.results, "a") as f:
                f.write("%s\t%s\t%s\n" % (n, tname, js))

    # -- the selection and the main flow -------------------------------------
    def why_not(self, s, atmos):
        """Every reason the case is excluded, one per line (none: selected)."""
        why = []
        if self.o["set"] != "all" and s != self.o["set"]:
            why.append("in set %s, not selected by --set %s" % (s, self.o["set"]))
        if atmos.startswith("bifrost:") and not self.o["bifrost"]:
            why.append("a bifrost case, and no --bifrost DIR or RH_BIFROST_DIR")
        return why

    def select(self):
        o = self.o
        if o["set"] not in ("conv", "perf", "all"):
            fail("run_stress.py: --set must be conv, perf or all, not '%s'" % o["set"], sys.stderr)
        if o["only"]:
            r = subprocess.run(["grep", "-Eq", "-e", o["only"]], stdin=subprocess.DEVNULL,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True)
            if r.returncode > 1:
                fail("run_stress.py: invalid --only pattern '%s': %s"
                     % (o["only"], r.stdout.rstrip("\n")), sys.stderr)
        sel, missed = [], []
        for c in CASES:
            n, s, g, a, rest, _ = fields(c)
            n = trim(n)
            if o["only"] and not grep_matches(o["only"], n):
                continue
            why = self.why_not(trim(s), trim(a))
            if not why:
                sel.append(c)
            elif o["only"]:
                missed.append("%s: %s" % (n, "; ".join(why)))
        if not sel:
            print("run_stress.py: no case selected by --set %s%s"
                  % (o["set"], " --only '%s'" % o["only"] if o["only"] else ""), file=sys.stderr)
            if missed:
                for m in missed:
                    print("  " + m, file=sys.stderr)
            elif o["only"]:
                print("  no case name matches '%s' (--list shows them)" % o["only"], file=sys.stderr)
            sys.exit(1)
        for m in missed:
            print("not run: " + m)
        return sel

    def main(self):
        o = self.o
        sel = self.select()
        self.names, self.dirs, srcs = [], [], []        # dirs: programs, srcs: source trees
        for t in o["trees"]:
            self.names.append(t.split("=")[0])
            d = t.split("=", 1)[1] if "=" in t else t
            if not os.path.isdir(d):
                fail("run_stress.py: not a directory: " + d, sys.stderr)
            d = os.path.abspath(d)
            srcs.append(tree.tree_src(d))
            d = tree.tree_bin(d)
            self.dirs.append(d)
            miss = tree.missing_programs(d)
            if miss:
                fail("missing binary: " + miss[0])
        self.ref_src = srcs[0]
        os.makedirs(os.path.join(o["work"], "runs"), exist_ok=True)
        wd = workdir.WorkDir(self.ref_src, o["work"], relink=True, announce=False)
        self.work = wd.path
        if not run.have_time():
            fail("need GNU /usr/bin/time")
        self.results = os.path.join(self.work, "results_%s.tsv" % time.strftime("%Y%m%d-%H%M%S"))
        open(self.results, "w").close()

        print("trees: " + "".join("%s=%s " % (n, d) for n, d in zip(self.names, self.dirs)))
        print("work : " + self.work)
        if not o["bifrost"]:
            print("bifrost cases skipped (no --bifrost DIR or RH_BIFROST_DIR)")
        # Pre-generate atmospheres serially (avoids races between parallel jobs)
        for c in sel:
            self.atmos_path(trim(fields(c)[3]), trim(fields(c)[1]))

        jobs = int(o["jobs"])
        slots = threading.Semaphore(jobs)
        pool = []

        def job(c, i):
            try:
                self.run_one(c, i)
            finally:
                slots.release()

        for c in sel:
            s = trim(fields(c)[1])
            for i in range(len(self.names)):
                if s == "perf":
                    self.run_one(c, i)
                else:
                    slots.acquire()
                    th = threading.Thread(target=job, args=(c, i))
                    th.start()
                    pool.append(th)
        for th in pool:
            th.join()
        self.report()

    # -- the table ---------------------------------------------------------
    def report(self):
        runs, names = os.path.join(self.work, "runs"), self.names
        rows = {}
        with open(self.results) as f:
            for line in f:
                case, tname, js = line.rstrip("\n").split("\t", 2)
                rows.setdefault(case, {})[tname] = json.loads(js)
        print(f"\n{'case':18s} {'tree':10s} {'iters':>5s} {'delta':>9s} {'conv':>4s} "
              f"{'trend':11s} {'wall[s]':>8s} {'RSS[MB]':>8s}  status / diff vs {names[0]}")
        for case, by in rows.items():
            for t in names:
                r = by.get(t)
                if r is None:
                    continue
                extra = r.get("status", "?")
                if t != names[0] and names[0] in by:
                    a, b = f"{runs}/{case}__{names[0]}", f"{runs}/{case}__{t}"
                    p = subprocess.run([sys.executable, "-I", os.path.join(HERE, "analyze_run.py"),
                                        "pops", a, b], stdout=subprocess.PIPE,
                                       universal_newlines=True).stdout.strip()
                    s = compare.compare_spectra(f"{a}/spectrum.out", f"{b}/spectrum.out",
                                                tol=1.0)[1].strip()
                    extra += f"; pops {p}; spec {s.split('(')[-1].split(')')[-1].strip() or s}"
                fmt = lambda v, f: (f % v) if isinstance(v, (int, float)) and v is not None else "-"
                w, fse, ref = r.get("work"), r.get("fse"), r.get("ref")
                if w:
                    extra += f"; FS {w['fs']:.1f}, Background {w['background_calls']}"
                if isinstance(fse, dict):
                    extra += f"; F_SE max@1e-6 {fse['max_rel@1e-06']:.1e} rms {fse['rms']:.1e}"
                if isinstance(ref, dict):
                    pe = ref.get("pops")
                    extra += "; vs ref: pops@1e-6 " + (fmt(pe.get("1e-06"), "%.1e") if isinstance(pe, dict) else str(pe))
                    extra += f", spec {fmt(ref.get('spec_wl'), '%.1e')} (line {fmt(ref.get('spec_line'), '%.1e')})"
                print(f"{case:18s} {t:10s} {fmt(r.get('iters'), '%d'):>5s} {fmt(r.get('delta'), '%.2e'):>9s} "
                      f"{'yes' if r.get('converged') else 'no':>4s} {r.get('trend', '-'):11s} "
                      f"{fmt(r.get('wall'), '%.1f'):>8s} {fmt(r.get('rss'), '%.0f'):>8s}  {extra}")


def main(argv):
    sys.stdout.reconfigure(line_buffering=True)
    o = parse_args(argv)
    if o["list"]:
        list_cases()
        sys.exit(0)
    if not o["trees"]:
        sys.stdout.write("\n".join(DOC.split("\n")[:9]) + "\n")
        sys.exit(1)
    try:
        if int(o["jobs"]) < 1:
            raise ValueError
    except ValueError:
        fail("run_stress.py: --jobs needs a positive integer, not '%s'" % o["jobs"], sys.stderr)
    Stress(o).main()


if __name__ == "__main__":
    main(sys.argv[1:])
