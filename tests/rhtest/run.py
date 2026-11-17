"""Run directories and running the RH programs."""
import os
import shutil
import subprocess
import sys
import threading

from . import keywords

# Exit statuses of the RH programs (error.h): 0, RH_EXIT_NOT_CONVERGED = 1
# (the run ended but not converged, expected when a test stops after a few
# iterations), RH_EXIT_FATAL = 70 (a fatal error, after its message), above
# 128 a signal. run_tests.sh accepts every flag combination 0-63; the e2e
# tests have only these two and say so.
EXIT_OK = 0
EXIT_NOT_CONVERGED = 1
EXIT_FATAL = 70
RAN = (EXIT_OK, EXIT_NOT_CONVERGED)     # the run went to its end


def template(src, geom):
    """The shipped example of a geometry, e.g. template(src, "rhf1d")."""
    return os.path.join(src, geom, "run")


def _remove(d):
    """Remove the directory d if it exists; a failed removal raises."""
    if os.path.lexists(d):
        shutil.rmtree(d)


def fresh_dir(work, name, create=True):
    """<work>/runs/<name>, removed if it exists and (create) made again."""
    d = work.join("runs", name)
    _remove(d)
    if create:
        os.makedirs(d)
    return d


def setup_dir(d, src, geom, pairs, atoms=None, copy=("molecules.input", "ray.input"),
              required=()):
    """Fill the run directory d from the template of the geometry geom
    (rhf1d, rhsc2d, ...): keyword.input with pairs applied (keywords.apply),
    atoms.input from atoms (keywords.atoms_input) if given, the files of
    copy that the template has, and the files of required, which it must
    have (else FileNotFoundError)."""
    t = template(src, geom)
    with open(os.path.join(t, "keyword.input")) as f:
        text = f.read()
    with open(os.path.join(d, "keyword.input"), "w") as f:
        f.write(keywords.apply(text, pairs))
    if atoms is not None:
        keywords.write_atoms_input(d, atoms)
    for f in copy:
        if os.path.exists(os.path.join(t, f)):
            shutil.copy(os.path.join(t, f), d)
    for f in required:
        shutil.copy(os.path.join(t, f), d)
    return d


def prepare(work, name, geom, pairs, atoms=None, copy=("molecules.input", "ray.input"),
            required=()):
    """fresh_dir + setup_dir: a run directory of the template of geom."""
    return setup_dir(fresh_dir(work, name), work.src, geom, pairs, atoms, copy, required)


def copy_dir(work, name, source):
    """A fresh copy of the directory source as <work>/runs/<name>."""
    d = fresh_dir(work, name, create=False)
    shutil.copytree(source, d)
    return d


class Result:
    """What execute() returns: program, dir, returncode, stdout (stdout and
    stderr merged, str)."""

    def __init__(self, program, d, returncode, stdout):
        self.program, self.dir = program, d
        self.returncode, self.stdout = returncode, stdout

    @property
    def ran(self):
        """The status is one of RAN."""
        return self.returncode in RAN

    def tail(self, n=40):
        return "\n".join(self.stdout.splitlines()[-n:])

    def require(self, ok=RAN, bad_text=None, tail=40):
        """Exit the test (status 1, message with the end of out.log, which
        carries the backtrace of a trapped floating-point exception, B-33)
        unless the status is in ok (None: any status) and bad_text (e.g.
        "ERROR") is not in the output. Returns self."""
        if (ok is not None and self.returncode not in ok) or (bad_text and bad_text in self.stdout):
            sys.exit(f"{os.path.basename(self.program)} failed in {self.dir} "
                     f"(status {self.returncode}); end of out.log:\n{self.tail(tail)}")
        return self


def _result(program, d, returncode, out, log, timed_out):
    """The common end of execute() and execute_streaming(): a timeout ends
    the test, the output goes to d/<log> (not written with log=None)."""
    if timed_out:
        sys.exit(f"{os.path.basename(program)} timed out in {d}")
    if log:
        with open(os.path.join(d, log), "w") as f:
            f.write(out)
    return Result(program, d, returncode, out)


def execute(program, d, args=(), timeout=600, log="out.log"):
    """Run program in the directory d, its output (stdout and stderr merged)
    in d/<log> (not written with log=None). A timeout (seconds, None for
    none) ends the test."""
    program = os.path.abspath(program)
    try:
        r = subprocess.run([program, *args], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return _result(program, d, None, "", log, True)
    return _result(program, d, r.returncode, r.stdout, log, False)


def execute_streaming(program, d, args=(), timeout=600, log="out.log", on_line=None):
    """execute(), reading the output line by line as the program writes it
    and calling on_line(line) for each, e.g. to create a STOP_RH file at a
    given point of the run. The timeout (seconds, None for none) kills the
    program; whatever ends the loop (an exception in on_line included), the
    timer is cancelled and the program does not outlive the call."""
    program = os.path.abspath(program)
    proc = subprocess.Popen([program, *args], cwd=d, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    timed_out = threading.Event()

    def kill():
        timed_out.set()
        proc.kill()

    timer = threading.Timer(timeout, kill) if timeout is not None else None
    if timer:
        timer.start()
    lines = []
    try:
        for line in proc.stdout:
            lines.append(line)
            if on_line:
                on_line(line)
        rc = proc.wait()
    finally:
        if timer:
            timer.cancel()
        if proc.poll() is None:
            proc.kill()
        proc.wait()
        proc.stdout.close()
    return _result(program, d, rc, "".join(lines), log, timed_out.is_set())


def make_field(src, geom, atmos, dst):
    """Write a magnetic field for the 2-D or 3-D atmosphere file atmos
    (tests/make_field.py; geom "2d" or "3d") to dst."""
    subprocess.run([sys.executable, "-I", os.path.join(src, "tests", "make_field.py"),
                    geom, atmos, dst], check=True, stdout=subprocess.DEVNULL)
