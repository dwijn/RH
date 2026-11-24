"""The command line the e2e scripts share."""
import argparse
import os

# tests/rhtest/cli.py -> the source tree
SRC = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def parse(required=(), optional=(), verbose=False, src_help=None, argv=None):
    """Parse --PROGRAM for the names in required and optional (paths of
    executables), --src (default: the source tree this file is in), --work,
    and, with verbose, -v. The paths and src come back absolute (optional
    programs not given stay None)."""
    p = argparse.ArgumentParser()
    for k in required:
        p.add_argument("--" + k, required=True)
    for k in optional:
        p.add_argument("--" + k)
    p.add_argument("--src", default=SRC,
                   help=src_help or "RH source tree (Atoms/, Atmos/, Molecules/, "
                                    "examples/)")
    p.add_argument("--work", help="run directory, kept (default: a temporary dir under "
                                  "$TMPDIR, removed when the test passes, kept and named "
                                  "on failure)")
    if verbose:
        p.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args(argv)
    for k in tuple(required) + tuple(optional) + ("src",):
        v = getattr(a, k)
        if v:
            setattr(a, k, os.path.abspath(v))
    return a
