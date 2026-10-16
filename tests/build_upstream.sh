#!/bin/bash
# Build an upstream RH checkout (han-uitenbroek/RH, which still uses its
# hand-written makefiles) as the baseline for run_tests.sh and
# run_stress.sh: rhf1d, rhsc2d, rhsc3d and rhsphere, built in place,
# without the Fortran Voigt routines so that no Fortran compiler is needed.
# Compare against it in tolerance mode: the Stokes cases differ by ~1e-10
# (complex arithmetic in the Voigt code). For byte-identical checks of a
# refactor, use a build of the branch's base as the baseline instead.
#
# usage: build_upstream.sh UPSTREAM_DIR
#
# e.g.  git clone https://github.com/han-uitenbroek/RH upstream
#       git -C upstream checkout c95bdf6
#       tests/build_upstream.sh upstream
# The makefiles are not safe under make -j, so this builds serially (about
# a minute).
set -e
[ -f "$1/Makefile" ] && [ -d "$1/rhf1d" ] ||
  { echo "usage: $0 UPSTREAM_DIR   (an upstream RH checkout)"; exit 1; }
cd "$1"
cc=${CC:-gcc}
cflags="-O2"; libs="-lrh -lm -lpthread"
if [ -d /usr/include/tirpc ]; then
  cflags="$cflags -I/usr/include/tirpc"; libs="$libs -ltirpc"
fi
make CC="$cc" CFLAGS="$cflags" librh.a > build_upstream.log 2>&1 ||
  { tail -20 build_upstream.log; exit 1; }
for g in rhf1d rhsc2d rhsc3d rhsphere; do
  (cd $g && make CC="$cc" LD="$cc" CFLAGS="$cflags -I../" LIBS="$libs" \
     -o librh $g >> ../build_upstream.log 2>&1) ||
    { tail -20 build_upstream.log; exit 1; }
done
echo "upstream build ready in $PWD"
