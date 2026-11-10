#!/bin/bash
# Reference solutions for the solver comparison (S-1 in docs/SOLVER.md):
# runs cases of run_stress.sh with the reference recipe and keeps their
# outputs, with a MANIFEST of the build and the inputs, in OUT/CASE/.
# run_stress.sh --ref OUT then reports every run's errors against them.
#
# usage: make_references.sh --tree DIR --out OUT [--only PATTERN] [--jobs N]
#                           [--bifrost DIR] [--kw 'K=V;...'] [--commit HASH]
#
#   --tree DIR     a build of RH (as for run_stress.sh). Its source tree's
#                  git commit is recorded; for a tree without .git (e.g. a
#                  git archive export) give it with --commit
#   --out OUT      where the references go (OUT/work holds the full runs)
#   --only, --jobs, --bifrost   passed to run_stress.sh (--set is conv)
#   --kw           extra keyword settings after the recipe (e.g. collisional
#                  relaxation for a case that needs it)
#
# Recipe: 1 thread in 1-D (run_stress.sh), the NG solver (SOLVER = NG, pinned
# so the recipe stays plain MALI when the default changes, S-10) with no
# acceleration (NG_ORDER = 0, PRD_NG_ORDER = 0), a fixed number of PRD
# sub-iterations (PRD_ITER_LIMIT = 0, so every main iteration does
# PRD_N_MAX_ITER of them, as S-2's G),
# ITER_LIMIT = 1e-10 (main iterations and the final background-scattering
# iteration, up to N_MAX_SCATTER = 50 sweeps), N_MAX_ITER = 3000. A case
# that does not reach ITER_LIMIT is kept but marked "reference: no" in its
# MANIFEST, and run_stress.sh --ref still compares with it: judge such
# cases on the SE residual alone.

HERE=$(cd "$(dirname "$0")" && pwd)
RECIPE='SOLVER=NG;NG_ORDER=0;PRD_NG_ORDER=0;PRD_ITER_LIMIT=0;ITER_LIMIT=1e-10;N_MAX_ITER=3000;N_MAX_SCATTER=50'
TREE=; OUT=; PASS=(); KW=; COMMIT=
while [ $# -gt 0 ]; do
  case "$1" in
    --tree)   TREE=$2; shift 2;;
    --out)    OUT=$2; shift 2;;
    --only|--jobs|--bifrost) PASS+=("$1" "$2"); shift 2;;
    --kw)     KW=$2; shift 2;;
    --commit) COMMIT=$2; shift 2;;
    -h|--help) sed -n '2,25p' "$0"; exit 0;;
    *) echo "unknown option: $1"; exit 1;;
  esac
done
[ -n "$TREE" ] && [ -n "$OUT" ] || { sed -n '2,25p' "$0"; exit 1; }
mkdir -p "$OUT" && OUT=$(cd "$OUT" && pwd) || exit 1
. "$HERE/tree.sh"
TREE=$(cd "$TREE" && pwd) || exit 1
BIN=$(tree_bin "$TREE"); SRC=$(tree_src "$BIN")
if [ -z "$COMMIT" ]; then
  COMMIT=$(git -C "$SRC" rev-parse HEAD 2>/dev/null) || { echo "no git commit for $SRC: give --commit"; exit 1; }
  [ -z "$(git -C "$SRC" status --porcelain --untracked-files=no)" ] || COMMIT="$COMMIT (modified)"
fi

"$HERE/run_stress.sh" --tree ref="$TREE" --set conv --work "$OUT/work" \
  --kw "$RECIPE${KW:+;$KW}" "${PASS[@]}" || exit 1
RESULTS=$(ls -t "$OUT"/work/results_*.tsv | head -1)

python3 -I - "$RESULTS" "$OUT" "$BIN" "$SRC" "$COMMIT" "$RECIPE${KW:+;$KW}" <<'EOF'
import hashlib, json, os, re, shutil, subprocess, sys, time
results, out, bindir, src, commit, recipe = sys.argv[1:]
KEEP = re.compile(r"^(pops|atom|radrate|collrate)\..+\.out$|^(spectrum|atmos|geometry|input|"
                  r"molecules|metals|flux)\.out$|^(keyword|atoms|molecules|ray|kurucz)\.input$|"
                  r"^log\.log$|^time\.txt$")
cache = {}
for l in open(os.path.join(bindir, "CMakeCache.txt")):
    m = re.match(r"(CMAKE_C_COMPILER|CMAKE_BUILD_TYPE|CMAKE_C_FLAGS|RH_\w+):\w+=(.*)", l)
    if m:
        cache[m.group(1)] = m.group(2)
cc = cache.get("CMAKE_C_COMPILER", "cc")
version = subprocess.run([cc, "--version"], capture_output=True, text=True).stdout.splitlines()[:1]
md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
for line in open(results):
    case, tree, js = line.rstrip("\n").split("\t", 2)
    r = json.loads(js)
    run = os.path.join(out, "work", "runs", f"{case}__{tree}")
    dest = os.path.join(out, case)
    shutil.rmtree(dest, ignore_errors=True)
    os.makedirs(dest)
    for f in os.listdir(run):
        if KEEP.match(f):
            shutil.copy2(os.path.join(run, f), dest)
    kw = {}
    for l in open(os.path.join(run, "keyword.input")):
        f = l.split("#")[0].split("=", 1)
        if len(f) == 2:
            kw[f[0].strip()] = f[1].strip()
    inputs = {}
    atmos = kw.get("ATMOS_FILE", "")
    for name in [atmos, kw.get("STOKES_INPUT", "none")] + [
            w.split()[0] for w in open(os.path.join(run, "atoms.input"))
            if "ACTIVE" in w.split()[1:2] and not w.lstrip().startswith("#")]:
        p = os.path.join(run, name)
        if name != "none" and os.path.exists(p):
            inputs[os.path.normpath(name)] = md5(p)
    ok = r.get("status") == "ok" and r.get("converged")
    manifest = dict(case=case, created=time.strftime("%Y-%m-%d %H:%M"), commit=commit,
                    build=cache, compiler=version[0] if version else cc,
                    recipe=recipe, reference="yes" if ok else "no",
                    iterations=r.get("iters"), last_change=r.get("delta"),
                    status=r.get("status"), fse=r.get("fse"), work=r.get("work"),
                    input_md5=inputs)
    json.dump(manifest, open(os.path.join(dest, "MANIFEST.json"), "w"), indent=1)
    print(f"{case:20s} reference {'yes' if ok else 'NO '}  {r.get('iters')} iterations, "
          f"last change {r.get('delta'):.1e}, F_SE max@1e-6 "
          f"{r['fse']['max_rel@1e-06']:.1e}" if isinstance(r.get("fse"), dict) else
          f"{case:20s} reference NO   status {r.get('status')}")
EOF
