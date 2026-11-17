#!/bin/bash
# Regression tests for RH: compare a candidate build against an unmodified
# baseline build on small, fast runs.
#
# usage: run_tests.sh --baseline DIR --candidate DIR [options]
#
#   --baseline DIR   the reference build (e.g. upstream RH built with its
#                    makefiles). DIR is a source tree with the programs built
#                    in place or in build/, or a CMake build directory.
#   --candidate DIR  the build under test, given the same way. Example
#                    inputs, Atoms/, Atmos/ and Molecules/ are taken from its
#                    source tree.
#   --asan DIR       also run the candidate cases with an AddressSanitizer
#                    build, e.g. a CMake build directory configured with
#                    -DRH_SANITIZE=address
#   --exact          require byte-identical output between baseline and
#                    candidate at 1 thread. Use for pure refactors / memory
#                    management; leave off when summation order may change.
#   --tol X          allowed relative difference otherwise. Defaults are set per
#                    geometry from the noise of unmodified code against itself
#                    at different thread counts: 1e-9 (1-D, 3-D) and 1e-7 (2-D,
#                    whose shipped example is ill-conditioned, ~2e-9 noise)
#   --threads "A B"  thread counts to check against the candidate's own
#                    1-thread result (default "3 8": an uneven split of
#                    the work and more threads than cores)
#   --memtol X       candidate peak memory may be at most X times the
#                    baseline's in the memory case (default 1.05)
#   --only PATTERN   only run cases whose name contains PATTERN
#   --skip PATTERN   skip cases whose name contains PATTERN
#   --work DIR      scratch directory (default tests/work)
#
# What is checked, per case (1-D, 2-D, 3-D and spherical geometries, with and without
# Stokes, CRD only, two active atoms, ...):
#   * candidate at 1 thread runs cleanly and matches the baseline
#   * candidate at each thread count runs cleanly and matches its own
#     1-thread result (up to floating-point summation order)
#   * with --asan: no sanitizer reports at 1 thread and the largest count
#     (3-D at the largest count only)
# plus a memory-scaling case and a Kurucz line list with a leading blank line.
#
# The code reads atomic data through hard-coded ../../Atoms paths, so runs
# are made in WORK/runs/<name>, with WORK/Atoms etc. symlinked to the
# candidate tree.

HERE=$(cd "$(dirname "$0")" && pwd)
BASE=; CAND=; ASAN=; EXACT=0; TOL=; THREADS="3 8"; MEMTOL=1.05
ONLY=; SKIP=; WORK="$HERE/work"

while [ $# -gt 0 ]; do
  case "$1" in
    --baseline)  BASE=$2; shift 2;;
    --candidate) CAND=$2; shift 2;;
    --asan)      ASAN=$2; shift 2;;
    --exact)     EXACT=1; shift;;
    --tol)       TOL=$2; shift 2;;
    --threads)   THREADS=$2; shift 2;;
    --memtol)    MEMTOL=$2; shift 2;;
    --only)      ONLY=$2; shift 2;;
    --skip)      SKIP=$2; shift 2;;
    --work)      WORK=$2; shift 2;;
    -h|--help)   sed -n '2,43p' "$0"; exit 0;;
    *) echo "unknown option: $1"; exit 1;;
  esac
done
[ -n "$BASE" ] && [ -n "$CAND" ] || { sed -n '2,6p' "$0"; exit 1; }
. "$HERE/tree.sh"
BASE=$(cd "$BASE" && pwd); CAND=$(cd "$CAND" && pwd)
[ -z "$ASAN" ] || ASAN=$(cd "$ASAN" && pwd)
BASE_BIN=$(tree_bin "$BASE"); CAND_BIN=$(tree_bin "$CAND")
if [ -n "$ASAN" ]; then
  ASAN_BIN=$(tree_bin "$ASAN")
  tree_sanitized "$ASAN" || [ ! -f "$ASAN_BIN/CMakeCache.txt" ] ||
    { echo "--asan $ASAN: $ASAN_BIN is not configured with -DRH_SANITIZE"; exit 1; }
fi
BASE=$(tree_src "$BASE"); CAND=$(tree_src "$CAND")
mkdir -p "$WORK/runs"; WORK=$(cd "$WORK" && pwd)
for d in Atoms Atmos Molecules; do ln -sfn "$CAND/$d" "$WORK/$d"; done
# Magnetic fields for the 2-D and 3-D Stokes cases (make_field.py)
python3 -I "$HERE/make_field.py" 2d "$CAND/Atmos/FALC_03x82.atmos" "$WORK/field_2d.B" >/dev/null &&
python3 -I "$HERE/make_field.py" 3d "$CAND/Atmos/FALC_04x04x82.atmos" "$WORK/field_3d.B" >/dev/null ||
  echo "make_field.py failed: the 2d/3d fieldfree cases will fail"
THREADS_MAX=$(echo $THREADS | tr ' ' '\n' | sort -n | tail -1)

for t in "$BASE_BIN" "$CAND_BIN" ${ASAN:+"$ASAN_BIN"}; do
  for g in rhf1d rhsc2d rhsc3d rhsphere; do
    [ -x "$t/$g/$g" ] || { echo "missing binary: $t/$g/$g"; exit 1; }
  done
done

if [ "$(uname)" = Darwin ]; then TIMEFLAG=-l; else TIMEFLAG=-v; fi
[ -x /usr/bin/time ] || echo "warning: /usr/bin/time missing, memory case skipped"

NPASS=0; NFAIL=0; NSKIP=0; FAILED=()
record() {   # PASS|FAIL|SKIP  name  detail
  local keep
  case "$1" in
    PASS) NPASS=$((NPASS+1));;
    FAIL) NFAIL=$((NFAIL+1)); FAILED+=("$2");;
    SKIP) NSKIP=$((NSKIP+1));;
  esac
  printf '[%s] %-34s %s\n' "$1" "$2" "$3"
  # Keep the log of a failed run (the run directory is reused by later
  # cases), so that intermittent failures can be looked at afterwards
  if [ "$1" = FAIL ] && [ -n "$RUNDIR" ] && [ -f "$RUNDIR/log.log" ]; then
    keep="$WORK/failed/$(echo "$2" | tr ' /=' '___')_$(date +%H%M%S)"
    mkdir -p "$keep"
    cp "$RUNDIR/log.log" "$RUNDIR/keyword.input" "$keep/" 2>/dev/null
    cp "$RUNDIR/time.txt" "$keep/" 2>/dev/null
    printf '       log kept in %s\n' "$keep"
  fi
}

# setkw FILE KEY VALUE: set a keyword.input entry, appending it if absent
setkw() {
  if grep -Eq "^[[:space:]]*$2[[:space:]]*=" "$1"; then
    sed -E -i.bak "s|^[[:space:]]*$2[[:space:]]*=.*|  $2 = $3|" "$1" && rm -f "$1.bak"
  else
    printf '  %s = %s\n' "$2" "$3" >> "$1"
  fi
}

geomdir() { case $1 in 1d) echo rhf1d;; 2d) echo rhsc2d;; 3d) echo rhsc3d;; sph) echo rhsphere;; esac; }

# prepare GEOM CASE LABEL NTHREADS MODS  -> sets RUNDIR
prepare() {
  local geom=$1 name=$2 label=$3 n=$4 mods=$5 m
  RUNDIR="$WORK/runs/${label}_${geom}_${name}_n$n"
  rm -rf "$RUNDIR"; mkdir -p "$RUNDIR"
  cp "$CAND/$(geomdir $geom)/run/"* "$RUNDIR/"
  setkw "$RUNDIR/keyword.input" N_THREADS "$n"
  setkw "$RUNDIR/keyword.input" N_MAX_ITER 4
  local IFS=';'
  for m in $mods; do
    case "$m" in
      @H6_ACTIVE) sed -E -i.bak 's/(H_6\.atom[[:space:]]+)PASSIVE.*/\1ACTIVE     ZERO_RADIATION    pops.H.out/' "$RUNDIR/atoms.input"
                  rm -f "$RUNDIR/atoms.input.bak";;
      ?*=*)       setkw "$RUNDIR/keyword.input" "${m%%=*}" "${m#*=}";;
    esac
  done
}

# completed RC: did the program complete? 0 = converged; 1-63 = completed
# with problems (bit 1: populations not converged, which is expected here,
# since the cases stop after 4 iterations); 70 = fatal error; >128 = signal
completed() { [ "$1" -ge 0 ] && [ "$1" -le 63 ]; }
# Sanitizers exit with 99 when they report an error (their default, 1,
# would look like "not converged")
SAN_OPTIONS=detect_leaks=0:exitcode=99   # ASan; UBSan gets exitcode=99

# execute TREE GEOM -> runs in RUNDIR, returns the program's exit status
execute() {
  local bin="$1/$(geomdir $2)/$(geomdir $2)"
  if [ -x /usr/bin/time ]; then
    ( cd "$RUNDIR" && ASAN_OPTIONS=$SAN_OPTIONS UBSAN_OPTIONS=exitcode=99 \
        /usr/bin/time $TIMEFLAG -o time.txt "$bin" > log.log 2>&1 )
  else
    ( cd "$RUNDIR" && ASAN_OPTIONS=$SAN_OPTIONS UBSAN_OPTIONS=exitcode=99 \
        "$bin" > log.log 2>&1 )
  fi
}

log_problem() {   # prints a reason if the log shows a failure
  grep -Eqi 'TERMINATING_ERROR|AddressSanitizer|runtime error|Segmentation' "$RUNDIR/log.log" \
    && { grep -Ei -m1 'TERMINATING_ERROR|AddressSanitizer|runtime error' -A1 "$RUNDIR/log.log" | tr '\n' ' ' | cut -c1-150; return; }
  grep -Eq '(^|[^A-Za-z])NAN([^A-Za-z]|$)' "$RUNDIR/log.log" && echo "NaN in log"
}

has_stokes() {    # does this run read a magnetic field (=> Q,U,V in output)?
  grep -Eq '^[[:space:]]*STOKES_INPUT[[:space:]]*=' "$1/keyword.input" &&
    ! grep -Eq '^[[:space:]]*STOKES_INPUT[[:space:]]*=[[:space:]]*none' "$1/keyword.input"
}

tol_for() {
  if [ -n "$TOL" ]; then echo "$TOL"; else case $1 in 2d) echo 1e-7;; *) echo 1e-9;; esac; fi
}

compare() {       # A_DIR B_DIR GEOM [exact]  -> sets CMPOUT, returns status
  local flags="--tol $(tol_for $3)"
  has_stokes "$1" && flags="$flags --stokes"
  [ "$4" = exact ] && flags="--exact ${flags}"
  CMPOUT=$(python3 -I "$HERE/compare_spectra.py" "$1/spectrum.out" "$2/spectrum.out" $flags 2>&1)
}

# MB of the background file that the program in bin dir $1 maps read-only
# (P-18, recognised by its warning text) in run dir $2: the pages it reads
# count in the peak RSS but are clean, shared page cache, which pread used
# to copy from, not memory the program allocates. 0 for a program that
# reads the file with pread.
mapped_mb() {
  local bin="$1/$(geomdir 2d)/$(geomdir 2d)"
  if grep -aq 'Cannot map the background file' "$bin" && [ -f "$2/background.dat" ]; then
    awk -v s="$(wc -c < "$2/background.dat")" 'BEGIN{printf "%.0f", s/1048576}'
  else
    echo 0
  fi
}

rss_mb() {
  if [ "$TIMEFLAG" = -l ]; then
    awk '/maximum resident set size/{printf "%.0f", $1/1048576}' "$1/time.txt"
  else
    awk '/Maximum resident set size/{printf "%.0f", $NF/1024}' "$1/time.txt"
  fi
}

# ---------------------------------------------------------------------------
CASES=(
  "1d:nostokes:STOKES_MODE=NO_STOKES;STOKES_INPUT=none"
  "1d:fieldfree:STOKES_MODE=FIELD_FREE"
  "1d:fullstokes:STOKES_MODE=FULL_STOKES"
  "1d:crd:STOKES_MODE=NO_STOKES;STOKES_INPUT=none;PRD_N_MAX_ITER=0"
  "1d:twoactive:STOKES_MODE=NO_STOKES;STOKES_INPUT=none;@H6_ACTIVE"
  "2d:small:"
  "3d:small:N_MAX_ITER=2"
  "2d:fieldfree:STOKES_INPUT=$WORK/field_2d.B;STOKES_MODE=FIELD_FREE"
  "3d:fieldfree:STOKES_INPUT=$WORK/field_3d.B;STOKES_MODE=FIELD_FREE;N_MAX_ITER=2"
  "2d:lcbezier:STOKES_INPUT=$WORK/field_2d.B;STOKES_MODE=FIELD_FREE;S_INTERPOLATION_LONGCHAR=S_BEZIER3"
  "2d:fullstokes:STOKES_INPUT=$WORK/field_2d.B;STOKES_MODE=FULL_STOKES"
  "sph:small:"
)

run_case() {
  local geom=$1 name=$2 mods=$3 id="$1/$2" n rc why
  prepare $geom $name base 1 "$mods"; local B1=$RUNDIR
  execute "$BASE_BIN" $geom; rc=$?
  if ! completed $rc; then
    record SKIP "$id" "baseline itself failed (rc=$rc), nothing to compare against"; return
  fi

  prepare $geom $name cand 1 "$mods"; local C1=$RUNDIR
  execute "$CAND_BIN" $geom; rc=$?; why=$(log_problem)
  if ! completed $rc || [ -n "$why" ]; then
    record FAIL "$id candidate N=1" "rc=$rc $why"; return
  fi
  if [ $EXACT -eq 1 ]; then compare "$B1" "$C1" $geom exact; else compare "$B1" "$C1" $geom; fi
  [ $? -eq 0 ] && record PASS "$id candidate N=1 vs baseline" "$CMPOUT" \
               || record FAIL "$id candidate N=1 vs baseline" "$CMPOUT"

  for n in $THREADS; do
    prepare $geom $name cand $n "$mods"; local CN=$RUNDIR
    execute "$CAND_BIN" $geom; rc=$?; why=$(log_problem)
    if ! completed $rc || [ -n "$why" ]; then
      record FAIL "$id candidate N=$n" "rc=$rc $why"; continue
    fi
    compare "$C1" "$CN" $geom
    [ $? -eq 0 ] && record PASS "$id N=$n vs N=1" "$CMPOUT" \
                 || record FAIL "$id N=$n vs N=1" "$CMPOUT"
  done

  [ -n "$ASAN" ] || return
  # 3-D only at the largest count: the code that runs only with one thread
  # (Gamma and the rates accumulated in place, the background without
  # workers) is shared by all geometries and checked by the cheaper cases,
  # and 3-D under the sanitizers at 1 thread took 105 s of CI's slowest job
  local counts="1 $THREADS_MAX"
  [ "$geom" = 3d ] && counts=$THREADS_MAX
  for n in $counts; do
    prepare $geom $name asan $n "$mods"
    execute "$ASAN_BIN" $geom; rc=$?; why=$(log_problem)
    completed $rc && [ -z "$why" ] && record PASS "$id ASan N=$n" "clean" \
                                   || record FAIL "$id ASan N=$n" "rc=$rc $why"
  done
}

# Peak memory must not grow with the number of wavelengths on a large grid.
# (Per-wavelength caching of opacity buffers once added Nspect*Nspace*40 bytes.)
# Both runs use SOLVER = NG with the Ng settings the shipped 2-D keyword file
# had before S-11 made ANDERSON the 2-D default, so that a change of the
# default solver does not count as growth and Ng's history is still measured.
# 96 columns: about 115 MB, of which 8 MB do not grow with the grid; the old
# caching would add about 85% to the part that does, at any grid size.
# A program that maps the background file (P-18) has its pages in its RSS
# too (26 MB here); they are page cache, not allocations, and are not
# counted for either program (mapped_mb).
run_mem_case() {
  [ -x /usr/bin/time ] || { record SKIP "2d/mem" "no /usr/bin/time"; return; }
  python3 -I "$HERE/make_tiled_2d_atmos.py" "$CAND/Atmos/FALC_03x82.atmos" \
          "$WORK/tiled_2d.atmos" 96 >/dev/null || { record FAIL "2d/mem" "tiling failed"; return; }
  local mods="ATMOS_FILE=$WORK/tiled_2d.atmos;N_MAX_ITER=1;SOLVER=NG;NG_DELAY=6;NG_ORDER=2;NG_PERIOD=3;PRD_ITER_LIMIT=1.0E-2" rc why bm cm
  prepare 2d mem base 1 "$mods"; local B=$RUNDIR; execute "$BASE_BIN" 2d; rc=$?
  completed $rc || { record SKIP "2d/mem" "baseline failed (rc=$rc)"; return; }
  prepare 2d mem cand 1 "$mods"; local C=$RUNDIR; execute "$CAND_BIN" 2d; rc=$?; why=$(log_problem)
  completed $rc && [ -z "$why" ] || { record FAIL "2d/mem candidate" "rc=$rc $why"; return; }
  # Peak RSS without the mapped background file of either program (P-18)
  local bmap cmap
  bmap=$(mapped_mb "$BASE_BIN" "$B"); cmap=$(mapped_mb "$CAND_BIN" "$C")
  bm=$(( $(rss_mb "$B") - bmap )); cm=$(( $(rss_mb "$C") - cmap ))
  local note="limit x$MEMTOL; mapped background file not counted: baseline ${bmap} MB, candidate ${cmap} MB"
  if awk -v c="$cm" -v b="$bm" -v t="$MEMTOL" 'BEGIN{exit !(c <= b*t)}'; then
    record PASS "2d/mem peak memory" "candidate ${cm} MB vs baseline ${bm} MB ($note)"
  else
    record FAIL "2d/mem peak memory" "candidate ${cm} MB vs baseline ${bm} MB ($note)"
  fi
}

# A Kurucz list with a leading blank line (as downloaded from the current
# kurucz.cfa.harvard.edu) must give the same result as the file without it.
run_kurucz_blank_case() {
  local mods="STOKES_MODE=NO_STOKES;STOKES_INPUT=none" rc why
  prepare 1d kurucz cand 1 "$mods"; local A=$RUNDIR
  echo "../../Atoms/Kurucz/gf0400.10" > "$A/kurucz.input"
  execute "$CAND_BIN" 1d; rc=$?; why=$(log_problem)
  completed $rc && [ -z "$why" ] || { record FAIL "1d/kurucz_blank reference" "rc=$rc $why"; return; }
  prepare 1d kurucz_blank cand 1 "$mods"; local B=$RUNDIR
  { echo; cat "$CAND/Atoms/Kurucz/gf0400.10"; } > "$B/kurucz_blank.dat"
  echo "./kurucz_blank.dat" > "$B/kurucz.input"
  execute "$CAND_BIN" 1d; rc=$?; why=$(log_problem)
  completed $rc && [ -z "$why" ] || { record FAIL "1d/kurucz_blank" "rc=$rc $why"; return; }
  compare "$A" "$B" 1d exact
  [ $? -eq 0 ] && record PASS "1d/kurucz_blank" "same as without the blank line" \
               || record FAIL "1d/kurucz_blank" "$CMPOUT"
}

selected() {
  { [ -z "$ONLY" ] || [[ "$1" == *"$ONLY"* ]]; } &&
    { [ -z "$SKIP" ] || [[ "$1" != *"$SKIP"* ]]; }
}

echo "baseline : $BASE_BIN"
echo "candidate: $CAND_BIN${ASAN:+   (ASan: $ASAN_BIN)}"
echo "mode     : $([ $EXACT -eq 1 ] && echo byte-identical || echo "tolerance ${TOL:-per geometry}"), threads: $THREADS"
echo

for c in "${CASES[@]}"; do
  IFS=: read -r geom name mods <<<"$c"
  selected "$geom/$name" && run_case "$geom" "$name" "$mods"
done
selected 2d/mem && run_mem_case
selected 1d/kurucz_blank && run_kurucz_blank_case

echo
echo "passed $NPASS, failed $NFAIL, skipped $NSKIP"
[ $NFAIL -eq 0 ] || { printf 'FAILED: %s\n' "${FAILED[@]}"; exit 1; }
