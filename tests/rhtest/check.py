"""Pass/fail bookkeeping and comparisons."""
import os
import sys

from . import output


class Checker:
    """Records the checks of a test; finish() prints the summary and exits
    with status 0 (all passed) or 1."""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.n = 0
        self.failed = []

    def fail(self, label):
        """A check that failed."""
        self.n += 1
        self.failed.append(label)

    def check(self, ok, what, detail=""):
        """Record ok; print PASS/FAIL (with verbose, or when it failed)."""
        self.n += 1
        if not ok:
            self.failed.append(what)
        if self.verbose or not ok:
            print(f"{'PASS' if ok else 'FAIL'}  {what} {detail}".rstrip())
        return ok

    def report(self, ok, line, label=None):
        """Record ok and always print line; a failure is listed under label
        (default: line)."""
        self.n += 1
        if not ok:
            self.failed.append(label or line)
        print(line)
        return ok

    def value(self, name, value, tol, at_least=False):
        """value <= tol, or value >= tol with at_least."""
        ok = value >= tol if at_least else value <= tol
        self.n += 1
        if not ok:
            self.failed.append(name)
        if self.verbose or not ok:
            print(f"{'PASS' if ok else 'FAIL'}  {name}: {value:.3e} "
                  f"({'at least' if at_least else 'tol'} {tol:.1e})")
        return ok

    def finish(self, name, work=None):
        """Summary, then the end of the work directory (workdir.finish) and
        sys.exit."""
        ok = not self.failed
        if ok:
            print(f"{name}: all {self.n} checks passed")
        else:
            print(f"{name}: {len(self.failed)} of {self.n} checks FAILED: "
                  + "; ".join(self.failed))
        if work is not None:
            work.finish(ok)
        sys.exit(0 if ok else 1)


def _flat(x):
    for e in x:
        if isinstance(e, (int, float)):
            yield e
        else:
            yield from _flat(e)


def maxabs(x):
    """Largest |x| of a sequence, nested sequences included."""
    return max(abs(e) for e in _flat(x))


def maxdiff(x, y, sign=1.0):
    """Largest |x - sign * y| (sequences of equal shape, nesting allowed)."""
    return max(abs(u - sign * v) for u, v in zip(_flat(x), _flat(y)))


def maxrel(x, y):
    """Largest |x - y| / max(|x|, |y|, 1e-300)."""
    return max(abs(u - v) / max(abs(u), abs(v), 1e-300) for u, v in zip(x, y))


def same_file(d1, d2, f):
    """File f of the directories d1 and d2 is identical byte for byte."""
    return output.read_bytes(os.path.join(d1, f)) == output.read_bytes(os.path.join(d2, f))
