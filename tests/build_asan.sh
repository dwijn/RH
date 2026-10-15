#!/bin/bash
# Build an AddressSanitizer copy of an RH source tree for run_tests.sh --asan.
#
# usage: build_asan.sh SRC_TREE DEST_DIR
#
# SRC_TREE is copied (without .git or run directories) so the original
# build is left alone. No Fortran compiler is needed: CFLAGS is overridden
# on the command line, which drops -DHAVE_F90 and uses the C Voigt code.
set -e
src=$(cd "$1" && pwd); dest=$2
[ -n "$dest" ] || { echo "usage: $0 SRC_TREE DEST_DIR"; exit 1; }
mkdir -p "$dest"; dest=$(cd "$dest" && pwd)
rsync -a --exclude .git --exclude 'run*' --exclude 'en024048*' \
      --exclude '*.o' --exclude 'librh*.a' "$src"/ "$dest"/
F="-O1 -g -fsanitize=address -fno-omit-frame-pointer"
cd "$dest"
make librh.a CFLAGS="$F" >build_asan.log 2>&1 || { tail -20 build_asan.log; exit 1; }
for d in rhf1d rhsc2d rhsc3d rhsphere; do
  (cd $d && rm -f $d && make -o librh $d CFLAGS="$F -I../" \
        LDFLAGS="-fsanitize=address -L../" >>../build_asan.log 2>&1) \
    || { tail -20 build_asan.log; exit 1; }
done
echo "ASan build ready in $dest"
