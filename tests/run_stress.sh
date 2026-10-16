#!/bin/bash
# Stress, convergence and performance cases for RH: NLTE problems run to
# convergence (or to a fixed number of iterations for the timing cases) with
# one or two builds, reporting iterations, final population change, wall time,
# peak memory, and the difference between the builds.
#
# usage: run_stress.sh --tree NAME=DIR [--tree NAME=DIR ...] [options]
#
#   --tree NAME=DIR  a build of RH: a source tree with rhf1d, rhsc2d, rhsc3d
#                    and rhsphere built in place or in build/, or a CMake
#                    build directory. The first tree is the reference for
#                    comparisons.
#   --only PATTERN   only run cases whose name contains PATTERN (grep -E)
#   --set conv|perf|all   which group of cases to run (default conv)
#   --jobs N         run up to N cases at the same time (default 4). Timing
#                    cases always run one at a time.
#   --threads N      override N_THREADS of the multi-D cases
#   --work DIR       scratch directory (default tests/work_stress)
#   --list           print the cases and exit
#
# Results go to WORK/results_<date-time>.tsv (one file per invocation, one line
# per case and tree) and are printed as a table. Atomic data, atmospheres and
# molecules are taken from the first tree; generated multi-D atmospheres are
# written to WORK.

HERE=$(cd "$(dirname "$0")" && pwd)
TREES=(); ONLY=; SET=conv; JOBS=4; THREADS=; WORK="$HERE/work_stress"; LIST=0

while [ $# -gt 0 ]; do
  case "$1" in
    --tree)    TREES+=("$2"); shift 2;;
    --only)    ONLY=$2; shift 2;;
    --set)     SET=$2; shift 2;;
    --jobs)    JOBS=$2; shift 2;;
    --threads) THREADS=$2; shift 2;;
    --work)    WORK=$2; shift 2;;
    --list)    LIST=1; shift;;
    -h|--help) sed -n '2,22p' "$0"; exit 0;;
    *) echo "unknown option: $1"; exit 1;;
  esac
done

# ---------------------------------------------------------------------------
# Cases: NAME | SET | GEOM | ATMOS | ACTIVE ATOMS | KEYWORDS
#
# ATMOS is a file in Atmos/ or tests/atmos/, or gen:2d:NX / gen:3d:NX:NY for
# a perturbed FALC (make_perturbed_atmos.py). ACTIVE ATOMS is a comma-separated
# list of atom files in Atoms/ to make active; each replaces the line of the
# same element in the geometry's example atoms.input. KEYWORDS are
# ';'-separated keyword.input settings, plus @NOMOL=NAME to drop a molecule
# (e.g. MgH, which needs neutral Mg, when Mg II is the Mg model).
#
# 1-D cases run on 1 thread unless they set N_THREADS, so that two builds can
# be compared without summation-order noise (some cases, e.g. He_large, end at
# visibly different "converged" solutions for different thread counts).
#
# conv: run to convergence (N_MAX_ITER large), compare converged results.
# perf: fixed number of iterations, for wall time and peak memory.
# The comment after each group says what the case is known to expose.
L='ITER_LIMIT=1e-3;N_MAX_ITER=200'
CASES=(
# Well-behaved references (converge in 15-30 iterations)
"1d_falc_h6_ca_mg  |conv|1d|FALC_82.atmos  |H_6,CaII,MgII     |$L"
"1d_fala_h_ca      |conv|1d|FALA_80.atmos  |H,CaII            |$L"
"1d_falc_h20       |conv|1d|FALC_82.atmos  |H_20              |$L"
"1d_falc_mgi66     |conv|1d|FALC_82.atmos  |MgI_66            |$L"
"1d_falc_fe33      |conv|1d|FALC_82.atmos  |Fe33              |$L"
"1d_falc_nong      |conv|1d|FALC_82.atmos  |H_6,CaII,MgII     |$L;NG_ORDER=0"
"1d_falxco_h6crd   |conv|1d|FALXCO_80.atmos|H_6_CRD,CaII      |$L"
# Slow / non-convergence: PRD sub-iteration (H Lyman PRD), Ng without safeguards
"1d_falxco_h6_ca   |conv|1d|FALXCO_80.atmos|H_6,CaII          |$L"
"1d_t5750_h6       |conv|1d|T5750_g4.5_p00_ext.atmos|H_6     |$L"
"1d_falc_helarge   |conv|1d|FALC_82.atmos  |He_large          |$L"
"1d_falc_h6_tight  |conv|1d|FALC_82.atmos  |H_6               |ITER_LIMIT=1e-8;N_MAX_ITER=300;NG_ORDER=0"
"1d_falc_civi      |conv|1d|FALC_82.atmos  |CI-VI_Judge       |$L"
# Moving atmospheres: linear vs Bezier (they differ by 5-16 % in the
# chromosphere on this 82-point grid: depth resolution, not the scheme)
"1d_v10_ca_lin     |conv|1d|FALC_v10.atmos |CaII              |$L;S_INTERPOLATION=S_LINEAR"
"1d_v10_ca_bez     |conv|1d|FALC_v10.atmos |CaII              |$L;S_INTERPOLATION=S_BEZIER3"
"1d_shock_h6       |conv|1d|FALC_shock.atmos|H_6              |$L"
# Crashes / wrong answers from options
"1d_falc_cail      |conv|1d|FALC_82.atmos  |CaII_large        |$L"
"1d_falc_solvene   |conv|1d|FALC_82.atmos  |H_6,CaII          |$L;SOLVE_NE=ITERATION"
"1d_falc_prdng     |conv|1d|FALC_82.atmos  |CaII              |$L;PRD_NG_ORDER=2"
"1d_falxco_collrel |conv|1d|FALXCO_80.atmos|H_6               |$L;COLL_RELAX_NSTEP=20;COLL_RELAX_FACTOR=0.01"
# Multi-D
"2d_pert_h6_ca_mg  |conv|2d|gen:2d:32      |H_6,CaII,MgII     |$L;N_THREADS=8;S_INTERPOLATION=S_BEZIER3;@NOMOL=MgH"
"3d_pert_h6_ca     |conv|3d|gen:3d:8:8     |H_6,CaII          |$L;N_THREADS=8"
"2d_perf_many      |perf|2d|gen:2d:128     |H_6,CaII,MgII,Fe33|N_MAX_ITER=5;N_THREADS=16;S_INTERPOLATION=S_BEZIER3;@NOMOL=MgH"
"3d_perf_many      |perf|3d|gen:3d:24:24   |H_6,CaII,MgII,Fe33|N_MAX_ITER=3;N_THREADS=16;@NOMOL=MgH"
"3d_perf_mgi66     |perf|3d|gen:3d:16:16   |MgI_66            |N_MAX_ITER=3;N_THREADS=16"
)

trim() { local s=$1; s=${s#"${s%%[![:space:]]*}"}; echo "${s%"${s##*[![:space:]]}"}"; }
selected() { [ "$SET" = all ] || [ "$2" = "$SET" ] || return 1
             [ -z "$ONLY" ] || echo "$1" | grep -Eq "$ONLY"; }

if [ $LIST -eq 1 ]; then
  for c in "${CASES[@]}"; do IFS='|' read -r n s g a at kw <<<"$c"
    printf '%-18s %-4s %-3s %-15s %-28s %s\n' "$(trim $n)" "$(trim $s)" "$(trim $g)" \
           "$(trim $a)" "$(trim $at)" "$(trim "$kw")"; done
  exit 0
fi

[ ${#TREES[@]} -ge 1 ] || { sed -n '2,10p' "$0"; exit 1; }
. "$HERE/tree.sh"
NAMES=(); DIRS=(); SRCS=()      # DIRS: programs, SRCS: source trees
for t in "${TREES[@]}"; do
  NAMES+=("${t%%=*}"); d=$(cd "${t#*=}" && pwd) || exit 1
  SRCS+=("$(tree_src "$d")"); d=$(tree_bin "$d"); DIRS+=("$d")
  for g in rhf1d rhsc2d rhsc3d rhsphere; do
    [ -x "$d/$g/$g" ] || { echo "missing binary: $d/$g/$g"; exit 1; }
  done
done
REF=${SRCS[0]}
mkdir -p "$WORK/runs"; WORK=$(cd "$WORK" && pwd)
for d in Atoms Atmos Molecules; do ln -sfn "$REF/$d" "$WORK/$d"; done
[ -x /usr/bin/time ] || { echo "need GNU /usr/bin/time"; exit 1; }

geomdir() { case $1 in 1d) echo rhf1d;; 2d) echo rhsc2d;; 3d) echo rhsc3d;; sph) echo rhsphere;; esac; }

setkw() {
  if grep -Eq "^[[:space:]]*$2[[:space:]]*=" "$1"; then
    sed -E -i "s|^[[:space:]]*$2[[:space:]]*=.*|  $2 = $3|" "$1"
  else
    printf '  %s = %s\n' "$2" "$3" >> "$1"
  fi
}

atmos_path() {   # spec -> path of the atmosphere file (generates if needed)
  case $1 in
    gen:2d:*) local f="$WORK/pert_2d_${1#gen:2d:}.atmos"
              [ -f "$f" ] || python3 -I "$HERE/make_perturbed_atmos.py" \
                  "$REF/Atmos/FALC_03x82.atmos" "$f" ${1#gen:2d:} >/dev/null
              echo "$f";;
    gen:3d:*) local s=${1#gen:3d:} f="$WORK/pert_3d_${s/:/x}.atmos"
              [ -f "$f" ] || python3 -I "$HERE/make_perturbed_atmos.py" \
                  "$REF/Atmos/FALC_03x82.atmos" "$f" ${s%%:*} ${s##*:} >/dev/null
              echo "$f";;
    *)        if [ -f "$HERE/atmos/$1" ]; then echo "$HERE/atmos/$1"
              else echo "../../Atmos/$1"; fi;;
  esac
}

prepare() {   # NAME GEOM ATMOS ATOMS KEYWORDS TREENAME -> RUNDIR
  local name=$1 geom=$2 atm=$3 atoms=$4 kws=$5 tree=$6 a el m
  RUNDIR="$WORK/runs/${name}__$tree"
  rm -rf "$RUNDIR"; mkdir -p "$RUNDIR"
  cp "$REF/$(geomdir $geom)/run/"* "$RUNDIR/"
  # Plain intensity output unless a case asks otherwise
  setkw "$RUNDIR/keyword.input" ATMOS_FILE "$(atmos_path $atm)"
  setkw "$RUNDIR/keyword.input" STOKES_MODE NO_STOKES
  setkw "$RUNDIR/keyword.input" STOKES_INPUT none
  setkw "$RUNDIR/keyword.input" PRINT_CPU FALSE
  [ "$geom" = 1d ] && setkw "$RUNDIR/keyword.input" N_THREADS 1
  local IFS=';'
  for m in $kws; do
    case "$m" in
      @NOMOL=*) awk -v m="/${m#@NOMOL=}.molecule" '
                  index($0, m) { drop = 1; next }
                  { print }' "$RUNDIR/molecules.input" > "$RUNDIR/mol.tmp"
                # decrement the molecule count (first non-comment line)
                awk 'done || /^[[:space:]]*#/ || NF == 0 { print; next }
                     { print $1 - 1; done = 1 }' "$RUNDIR/mol.tmp" > "$RUNDIR/molecules.input"
                rm -f "$RUNDIR/mol.tmp";;
      ?*=*)     setkw "$RUNDIR/keyword.input" "${m%%=*}" "${m#*=}";;
    esac
  done
  [ -n "$THREADS" ] && [ "$geom" != 1d ] && setkw "$RUNDIR/keyword.input" N_THREADS "$THREADS"
  # All atoms passive, then the requested ones active, each replacing the
  # entry of the same element
  python3 -I - "$RUNDIR/atoms.input" "$REF/Atoms" "$atoms" <<'EOF' \
    2> "$RUNDIR/prepare_error" && rm -f "$RUNDIR/prepare_error"
import os, sys
path, atomdir, wanted = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
def element(f):
    for l in open(f, errors="replace"):
        if l.strip() and not l.lstrip().startswith("#"):
            return l.split()[0].upper()
lines = open(path).read().splitlines()
for i, l in enumerate(lines):
    lines[i] = l.replace(" ACTIVE ", " PASSIVE ").replace("\tACTIVE\t", "\tPASSIVE\t")
for a in filter(None, wanted):
    el = element(f"{atomdir}/{a}.atom")
    for i, l in enumerate(lines):
        w = l.split()
        if w and "Atoms/" in w[0] and element(f"{atomdir}/{os.path.basename(w[0])}") == el:
            lines[i] = f"  ../../Atoms/{a}.atom   ACTIVE   ZERO_RADIATION"
            break
    else:
        sys.exit(f"no {el} entry in atoms.input to replace with {a}")
open(path, "w").write("\n".join(lines) + "\n")
EOF
}

execute() {   # TREE_DIR GEOM  (in RUNDIR)
  local bin="$1/$(geomdir $2)/$(geomdir $2)"
  ( cd "$RUNDIR" && /usr/bin/time -v -o time.txt "$bin" > log.log 2>&1 )
}

RESULTS="$WORK/results_$(date +%Y%m%d-%H%M%S).tsv"
: > "$RESULTS"

run_one() {   # case-line tree-index
  local n s g a at kw lim
  IFS='|' read -r n s g a at kw <<<"$1"
  n=$(trim "$n"); g=$(trim "$g"); a=$(trim "$a"); at=$(trim "$at"); kw=$(trim "$kw")
  local tn=${NAMES[$2]} td=${DIRS[$2]}
  prepare "$n" "$g" "$a" "$at" "$kw" "$tn"
  if [ -f "$RUNDIR/prepare_error" ]; then
    printf '%s\t%s\t%s\n' "$n" "$tn" "{\"status\": \"$(cat "$RUNDIR/prepare_error")\"}" >> "$RESULTS"; return
  fi
  execute "$td" "$g"
  lim=$(awk -F= '/^[[:space:]]*ITER_LIMIT/{gsub(/ /,"",$2); print $2}' "$RUNDIR/keyword.input")
  printf '%s\t%s\t%s\n' "$n" "$tn" "$(python3 -I "$HERE/analyze_run.py" summary "$RUNDIR" "${lim:-1e-2}")" >> "$RESULTS"
}

echo "trees: $(for i in "${!NAMES[@]}"; do printf '%s=%s ' "${NAMES[$i]}" "${DIRS[$i]}"; done)"
echo "work : $WORK"
# Pre-generate atmospheres serially (avoids races between parallel jobs)
for c in "${CASES[@]}"; do IFS='|' read -r n s g a rest <<<"$c"
  selected "$(trim "$n")" "$(trim "$s")" && atmos_path "$(trim "$a")" >/dev/null; done

for c in "${CASES[@]}"; do
  IFS='|' read -r n s rest <<<"$c"; n=$(trim "$n"); s=$(trim "$s")
  selected "$n" "$s" || continue
  for i in "${!NAMES[@]}"; do
    if [ "$s" = perf ]; then
      run_one "$c" $i
    else
      while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
      run_one "$c" $i &
    fi
  done
done
wait

# --- Report ----------------------------------------------------------------
python3 -I - "$RESULTS" "$WORK/runs" "$HERE" "${NAMES[@]}" <<'EOF'
import json, subprocess, sys
res, runs, here, names = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
rows = {}
for line in open(res):
    case, tree, js = line.rstrip("\n").split("\t", 2)
    rows.setdefault(case, {})[tree] = json.loads(js)
print(f"\n{'case':18s} {'tree':10s} {'iters':>5s} {'delta':>9s} {'conv':>4s} "
      f"{'trend':11s} {'wall[s]':>8s} {'RSS[MB]':>8s}  status / diff vs {names[0]}")
for case, by in rows.items():
    for t in names:
        r = by.get(t)
        if r is None: continue
        extra = r.get("status", "?")
        if t != names[0] and names[0] in by:
            a, b = f"{runs}/{case}__{names[0]}", f"{runs}/{case}__{t}"
            p = subprocess.run([sys.executable, "-I", f"{here}/analyze_run.py", "pops", a, b],
                               capture_output=True, text=True).stdout.strip()
            s = subprocess.run([sys.executable, "-I", f"{here}/compare_spectra.py",
                                f"{a}/spectrum.out", f"{b}/spectrum.out", "--tol", "1"],
                               capture_output=True, text=True).stdout.strip()
            extra += f"; pops {p}; spec {s.split('(')[-1].split(')')[-1].strip() or s}"
        fmt = lambda v, f: (f % v) if isinstance(v, (int, float)) and v is not None else "-"
        print(f"{case:18s} {t:10s} {fmt(r.get('iters'), '%d'):>5s} {fmt(r.get('delta'), '%.2e'):>9s} "
              f"{'yes' if r.get('converged') else 'no':>4s} {r.get('trend', '-'):11s} "
              f"{fmt(r.get('wall'), '%.1f'):>8s} {fmt(r.get('rss'), '%.0f'):>8s}  {extra}")
EOF
