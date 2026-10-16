#!/bin/bash
# Build an upstream RH checkout (han-uitenbroek/RH, which still uses its
# hand-written makefiles) as the baseline for run_tests.sh and
# run_stress.sh: rhf1d, rhsc2d, rhsc3d and rhsphere, built in place,
# without the Fortran Voigt routines so that no Fortran compiler is needed
# and the results compare exactly with the C code of this tree.
#
# usage: build_upstream.sh UPSTREAM_DIR
#
# e.g.  git clone https://github.com/han-uitenbroek/RH upstream
#       git -C upstream checkout c95bdf6
#       tests/build_upstream.sh upstream
# The makefiles are not safe under make -j, so this builds serially (about
# a minute).
set -e
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
