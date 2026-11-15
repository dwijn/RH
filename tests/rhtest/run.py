"""Run directories and running the RH programs."""
import os
import shutil
import subprocess
import sys

from . import keywords

# Exit statuses of the RH programs (error.h): bits 0-5 are flags of a run
# that completed (bit 0, RH_EXIT_NOT_CONVERGED = 1: not converged, which is
# expected when a test stops after a few iterations), 70 = RH_EXIT_FATAL
# (a fatal error, after its message), above 128 a signal. Sanitizer builds
# are not run by the e2e tests, so there is no other status.
EXIT_NOT_CONVERGED = 1
EXIT_FATAL = 70


def completed(status):
    """The program ran to its end (run_tests.sh completed())."""
    return 0 <= status <= 63


def template(src, geom):
    """The shipped example of a geometry, e.g. template(src, "rhf1d")."""
    return os.path.join(src, geom, "run")


def fresh_dir(work, name, create=True):
    """<work>/runs/<name>, removed if it exists and (create) made again."""
    d = work.join("runs", name)
    shutil.rmtree(d, ignore_errors=True)
    if create:
        os.makedirs(d)
    return d


def setup_dir(d, src, geom, pairs, atoms=None, copy=("molecules.input", "ray.input")):
    """Fill the run directory d from the template of the geometry geom
    (rhf1d, rhsc2d, ...): keyword.input with pairs applied (keywords.apply),
    atoms.input from atoms (keywords.atoms_input) if given, and the files of
    copy that the template has."""
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
    return d


def prepare(work, name, geom, pairs, atoms=None, copy=("molecules.input", "ray.input")):
    """fresh_dir + setup_dir: a run directory of the template of geom."""
    return setup_dir(fresh_dir(work, name), work.src, geom, pairs, atoms, copy)


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

    def tail(self, n=40):
        return "\n".join(self.stdout.splitlines()[-n:])

    def require(self, ok=(0, EXIT_NOT_CONVERGED), bad_text=None, tail=40):
        """Exit the test (status 1, message with the end of out.log, which
        carries the backtrace of a trapped floating-point exception, B-33)
        unless the status is in ok and bad_text (e.g. "ERROR") is not in the
        output. Returns self."""
        if self.returncode not in ok or (bad_text and bad_text in self.stdout):
            sys.exit(f"{os.path.basename(self.program)} failed in {self.dir} "
                     f"(status {self.returncode}); end of out.log:\n{self.tail(tail)}")
        return self


def execute(program, d, args=(), timeout=600):
    """Run program in the directory d, its output in d/out.log. A timeout
    (seconds, None for none) ends the test."""
    program = os.path.abspath(program)
    try:
        r = subprocess.run([program, *args], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        sys.exit(f"{os.path.basename(program)} timed out in {d}")
    with open(os.path.join(d, "out.log"), "w") as f:
        f.write(r.stdout)
    return Result(program, d, r.returncode, r.stdout)
