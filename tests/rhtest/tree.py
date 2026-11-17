"""Locating the programs and the data of an RH tree (the logic of tests/tree.sh,
which make_references.sh still sources).

A tree is given as
  * a CMake build directory (it contains CMakeCache.txt): the programs are
    there, the data (Atoms/, example runs, ...) in its source directory;
  * a source tree with a CMake build in build/;
  * a source tree built in place with the old makefiles (rhf1d/rhf1d, ...),
    e.g. upstream RH as a baseline (tests/build_upstream.sh).
"""
import os
import re

PROGRAMS = ("rhf1d", "rhsc2d", "rhsc3d", "rhsphere")
# geometry (1d, 2d, 3d, sph) -> program and directory
GEOMDIR = {"1d": "rhf1d", "2d": "rhsc2d", "3d": "rhsc3d", "sph": "rhsphere"}


def _cache(d):
    return os.path.join(d, "CMakeCache.txt")


def tree_bin(d):
    """The directory holding rhf1d/rhf1d etc."""
    if os.path.isfile(_cache(d)):
        return d
    if os.path.isfile(_cache(os.path.join(d, "build"))):
        return os.path.join(d, "build")
    return d


def _cache_lines(d):
    with open(_cache(tree_bin(d)), errors="replace") as f:
        return f.read().split("\n")


def tree_sanitized(d):
    """A CMake build with RH_SANITIZE set (to anything)."""
    if not os.path.isfile(_cache(tree_bin(d))):
        return False
    return any(re.match(r"RH_SANITIZE:STRING=.+", l) for l in _cache_lines(d))


def tree_src(d):
    """The source tree (the one a build directory was configured from)."""
    if os.path.isfile(_cache(d)):
        prefix = "CMAKE_HOME_DIRECTORY:INTERNAL="
        for l in _cache_lines(d):
            if l.startswith(prefix):
                return l[len(prefix):]
        return ""
    return d


def missing_programs(bin_dir):
    """The programs bin_dir/<p>/<p> that are not executable."""
    return [os.path.join(bin_dir, p, p) for p in PROGRAMS
            if not os.access(os.path.join(bin_dir, p, p), os.X_OK)]
