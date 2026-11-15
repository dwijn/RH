"""Shared plumbing of the end-to-end test scripts (tests/e2e/*.py), Q-31.

Python standard library only (the CI build jobs have no numpy). A script
finds the package by putting the directory tests/ on sys.path:

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

Modules:
    cli       the common command line (--src, --work, -v, program paths)
    workdir   the work directory: lifecycle and the Atoms/Atmos/Molecules links
    keywords  keyword.input editing (setkw), atoms.input writing
    run       run directories from a geometry's template, executing programs,
              RH exit statuses
    output    readers of RH's binary outputs
    check     pass/fail bookkeeping, comparisons, same_file
"""
