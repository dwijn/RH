# Locate the programs and the data of an RH tree; sourced by run_tests.sh
# and run_stress.sh. A tree argument can be
#   * a CMake build directory (contains CMakeCache.txt): programs there,
#     data (Atoms/, example runs, ...) in its source directory;
#   * a source tree with a CMake build in build/;
#   * a source tree built in place with the old makefiles (rhf1d/rhf1d,
#     ...), e.g. upstream RH as the baseline.

tree_bin() {   # DIR -> directory holding rhf1d/rhf1d etc.
  if [ -f "$1/CMakeCache.txt" ]; then echo "$1"
  elif [ -f "$1/build/CMakeCache.txt" ]; then echo "$1/build"
  else echo "$1"; fi
}

tree_sanitized() {   # DIR -> true if a CMake build with RH_SANITIZE set
  local c="$(tree_bin "$1")/CMakeCache.txt"
  [ -f "$c" ] && grep -Eq '^RH_SANITIZE:STRING=.+' "$c"
}

tree_src() {   # DIR -> source tree
  if [ -f "$1/CMakeCache.txt" ]; then
    sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$1/CMakeCache.txt"
  else echo "$1"; fi
}
