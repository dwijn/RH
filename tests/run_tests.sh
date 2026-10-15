#!/bin/bash
# Regression tests for RH: compare a candidate build against an unmodified
# baseline build on small, fast runs.
#
# usage: run_tests.sh --baseline DIR --candidate DIR [options]
#
#   --baseline DIR   source tree of the reference build (e.g. a checkout of
#                    upstream master) with rhf1d/rhf1d, rhsc2d/rhsc2d and
#                    rhsc3d/rhsc3d already built
#   --candidate DIR  source tree under test, built the same way. Example
#                    inputs, Atoms/, Atmos/ and Molecules/ are taken from here.
#   --asan DIR       also run the candidate cases with an AddressSanitizer
#                    build of the candidate (create it with build_asan.sh)
#   --exact          require byte-identical output between baseline and
#                    candidate at 1 thread. Use for pure refactors / memory
#                    management; leave off when summation order may change.
#   --tol X          allowed relative difference otherwise. Defaults are set per
#                    geometry from the noise of unmodified code against itself
#                    at different thread counts: 1e-9 (1-D, 3-D) and 1e-7 (2-D,
#                    whose shipped example is ill-conditioned, ~2e-9 noise)
#   --threads "A B"  thread counts to check against the candidate's own
#                    1-thread result (default "2 3 8")
#   --memtol X       candidate peak memory may be at most X times the
#                    baseline's in the memory case (default 1.05)
#   --only PATTERN   only run cases whose name contains PATTERN
#   --work DIR       scratch directory (default tests/work)
#
# What is checked, per case (1-D, 2-D, 3-D and spherical geometries, with and without
# Stokes, CRD only, two active atoms, ...):
#   * candidate at 1 thread runs cleanly and matches the baseline
#   * candidate at each thread count runs cleanly and matches its own
#     1-thread result (up to floating-point summation order)
#   * with --asan: no sanitizer reports at 1 thread and the largest count
# plus a memory-scaling case and a Kurucz line list with a leading blank line.
#
# The code reads atomic data through hard-coded ../../Atoms paths, so runs
# are made in WORK/runs/<name>, with WORK/Atoms etc. symlinked to the
# candidate tree.

HERE=$(cd "$(dirname "$0")" && pwd)
BASE=; CAND=; ASAN=; EXACT=0; TOL=; THREADS="2 3 8"; MEMTOL=1.05
ONLY=; WORK="$HERE/work"

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
    --work)      WORK=$2; shift 2;;
    -h|--help)   sed -n '2,38p' "$0"; exit 0;;
    *) echo "unknown option: $1"; exit 1;;
  esac
done
[ -n "$BASE" ] && [ -n "$CAND" ] || { sed -n '2,6p' "$0"; exit 1; }
BASE=$(cd "$BASE" && pwd); CAND=$(cd "$CAND" && pwd)
[ -z "$ASAN" ] || ASAN=$(cd "$ASAN" && pwd)
mkdir -p "$WORK/runs"; WORK=$(cd "$WORK" && pwd)
for d in Atoms Atmos Molecules; do ln -sfn "$CAND/$d" "$WORK/$d"; done
THREADS_MAX=$(echo $THREADS | tr ' ' '\n' | sort -n | tail -1)

for t in "$BASE" "$CAND" ${ASAN:+"$ASAN"}; do
  for g in rhf1d rhsc2d rhsc3d rhsphere; do
    [ -x "$t/$g/$g" ] || { echo "missing binary: $t/$g/$g"; exit 1; }
  done
done

if [ "$(uname)" = Darwin ]; then TIMEFLAG=-l; else TIMEFLAG=-v; fi
[ -x /usr/bin/time ] || echo "warning: /usr/bin/time missing, memory case skipped"

NPASS=0; NFAIL=0; NSKIP=0; FAILED=()
record() {   # PASS|FAIL|SKIP  name  detail
  case "$1" in
    PASS) NPASS=$((NPASS+1));;
    FAIL) NFAIL=$((NFAIL+1)); FAILED+=("$2");;
    SKIP) NSKIP=$((NSKIP+1));;
  esac
  printf '[%s] %-34s %s\n' "$1" "$2" "$3"
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
      @H6_ACTIVE) sed -E -i.bak 's/(H_6\.atom[[:space:]]+)PASSIVE/\1ACTIVE/' "$RUNDIR/atoms.input"
                  rm -f "$RUNDIR/atoms.input.bak";;
      ?*=*)       setkw "$RUNDIR/keyword.input" "${m%%=*}" "${m#*=}";;
    esac
  done
}

# execute TREE GEOM -> runs in RUNDIR, returns the program's exit status
execute() {
  local bin="$1/$(geomdir $2)/$(geomdir $2)"
  if [ -x /usr/bin/time ]; then
    ( cd "$RUNDIR" && ASAN_OPTIONS=detect_leaks=0 \
        /usr/bin/time $TIMEFLAG -o time.txt "$bin" > log.log 2>&1 )
  else
    ( cd "$RUNDIR" && ASAN_OPTIONS=detect_leaks=0 "$bin" > log.log 2>&1 )
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
  "sph:small:"
)

run_case() {
  local geom=$1 name=$2 mods=$3 id="$1/$2" n rc why
  prepare $geom $name base 1 "$mods"; local B1=$RUNDIR
  execute "$BASE" $geom; rc=$?
  if [ $rc -ne 0 ]; then
    record SKIP "$id" "baseline itself failed (rc=$rc), nothing to compare against"; return
  fi

  prepare $geom $name cand 1 "$mods"; local C1=$RUNDIR
  execute "$CAND" $geom; rc=$?; why=$(log_problem)
  if [ $rc -ne 0 ] || [ -n "$why" ]; then
    record FAIL "$id candidate N=1" "rc=$rc $why"; return
  fi
  if [ $EXACT -eq 1 ]; then compare "$B1" "$C1" $geom exact; else compare "$B1" "$C1" $geom; fi
  [ $? -eq 0 ] && record PASS "$id candidate N=1 vs baseline" "$CMPOUT" \
               || record FAIL "$id candidate N=1 vs baseline" "$CMPOUT"

  for n in $THREADS; do
    prepare $geom $name cand $n "$mods"; local CN=$RUNDIR
    execute "$CAND" $geom; rc=$?; why=$(log_problem)
    if [ $rc -ne 0 ] || [ -n "$why" ]; then
      record FAIL "$id candidate N=$n" "rc=$rc $why"; continue
    fi
    compare "$C1" "$CN" $geom
    [ $? -eq 0 ] && record PASS "$id N=$n vs N=1" "$CMPOUT" \
                 || record FAIL "$id N=$n vs N=1" "$CMPOUT"
  done

  [ -n "$ASAN" ] || return
  for n in 1 $THREADS_MAX; do
    prepare $geom $name asan $n "$mods"
    execute "$ASAN" $geom; rc=$?; why=$(log_problem)
    [ $rc -eq 0 ] && [ -z "$why" ] && record PASS "$id ASan N=$n" "clean" \
                                   || record FAIL "$id ASan N=$n" "rc=$rc $why"
  done
}

# Peak memory must not grow with the number of wavelengths on a large grid.
# (Per-wavelength caching of opacity buffers once added Nspect*Nspace*40 bytes.)
run_mem_case() {
  [ -x /usr/bin/time ] || { record SKIP "2d/mem" "no /usr/bin/time"; return; }
  python3 -I "$HERE/make_tiled_2d_atmos.py" "$CAND/Atmos/FALC_03x82.atmos" \
          "$WORK/tiled_2d.atmos" 192 >/dev/null || { record FAIL "2d/mem" "tiling failed"; return; }
  local mods="ATMOS_FILE=$WORK/tiled_2d.atmos;N_MAX_ITER=1" rc why bm cm
  prepare 2d mem base 1 "$mods"; local B=$RUNDIR; execute "$BASE" 2d; rc=$?
  [ $rc -eq 0 ] || { record SKIP "2d/mem" "baseline failed (rc=$rc)"; return; }
  prepare 2d mem cand 1 "$mods"; local C=$RUNDIR; execute "$CAND" 2d; rc=$?; why=$(log_problem)
  [ $rc -eq 0 ] && [ -z "$why" ] || { record FAIL "2d/mem candidate" "rc=$rc $why"; return; }
  bm=$(rss_mb "$B"); cm=$(rss_mb "$C")
  if awk -v c="$cm" -v b="$bm" -v t="$MEMTOL" 'BEGIN{exit !(c <= b*t)}'; then
    record PASS "2d/mem peak memory" "candidate ${cm} MB vs baseline ${bm} MB (limit x$MEMTOL)"
  else
    record FAIL "2d/mem peak memory" "candidate ${cm} MB vs baseline ${bm} MB (limit x$MEMTOL)"
  fi
}

# A Kurucz list with a leading blank line (as downloaded from the current
# kurucz.cfa.harvard.edu) must give the same result as the file without it.
run_kurucz_blank_case() {
  local mods="STOKES_MODE=NO_STOKES;STOKES_INPUT=none" rc why
  prepare 1d kurucz cand 1 "$mods"; local A=$RUNDIR
  echo "../../Atoms/Kurucz/gf0400.10" > "$A/kurucz.input"
  execute "$CAND" 1d; rc=$?; why=$(log_problem)
  [ $rc -eq 0 ] && [ -z "$why" ] || { record FAIL "1d/kurucz_blank reference" "rc=$rc $why"; return; }
  prepare 1d kurucz_blank cand 1 "$mods"; local B=$RUNDIR
  { echo; cat "$CAND/Atoms/Kurucz/gf0400.10"; } > "$B/kurucz_blank.dat"
  echo "./kurucz_blank.dat" > "$B/kurucz.input"
  execute "$CAND" 1d; rc=$?; why=$(log_problem)
  [ $rc -eq 0 ] && [ -z "$why" ] || { record FAIL "1d/kurucz_blank" "rc=$rc $why"; return; }
  compare "$A" "$B" 1d exact
  [ $? -eq 0 ] && record PASS "1d/kurucz_blank" "same as without the blank line" \
               || record FAIL "1d/kurucz_blank" "$CMPOUT"
}

selected() { [ -z "$ONLY" ] || [[ "$1" == *"$ONLY"* ]]; }

echo "baseline : $BASE"
echo "candidate: $CAND${ASAN:+   (ASan: $ASAN)}"
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
