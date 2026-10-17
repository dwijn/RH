#!/usr/bin/env python3
"""Sync the issue tables of docs/MODERNIZATION.md with GitHub.

usage: make_issues.py issues  [--dry-run]   create missing issues, link IDs in the log
       make_issues.py project [--dry-run]   create/fill the "RH modernization" Project

Needs the GitHub CLI (`gh`), logged in with the `repo` and `project` scopes.

issues: every row of a table in section 3 whose first cell is an ID such as
B-1 or N-10 becomes one issue titled "[ID] <text>", with labels for area,
severity, plan step and whether the problem also exists upstream. Issues
that already exist (matched by the "[ID]" title prefix) are left alone, so the
command can be re-run after new rows are added. The ID cells in the log are
then turned into links to their issues.

project: creates the user-owned Project (if missing), its custom fields
(Area, Severity, Plan step) and adds every "[ID]" issue to it with those
fields set. Views (board, grouping) cannot be created from the CLI; set them
up once in the web UI.
"""
import json
import os
import re
import subprocess
import sys
import time

REPO = "dwijn/RH"
OWNER = "dwijn"
PROJECT = "RH modernization"
LOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "MODERNIZATION.md")
LOG_URL = f"https://github.com/{REPO}/blob/modernize/docs/MODERNIZATION.md"

AREA = {"B": "bug", "N": "numerics", "P": "speed", "M": "memory",
        "T": "threading", "Q": "quality", "F": "feature", "S": "solver"}
AREA_COLOR = {"bug": "b60205", "numerics": "5319e7", "speed": "0e8a16",
              "memory": "1d76db", "threading": "fbca04", "quality": "c5def5",
              "feature": "a2eeef", "solver": "8a2be2"}
SEV_COLOR = {"high": "d93f0b", "medium": "f9d0c4", "low": "ededed"}

# Plan step (section 4 of the log) for every ID
STEP = {}
for ids, step in [
    ("B-1 B-2 B-3 B-4 B-8 B-9 B-13 B-15 B-16 B-17 B-18 B-19 N-1 N-5 N-11", 1),
    ("Q-1 Q-2 Q-3 Q-10 Q-12 Q-13 Q-15 Q-17", 2),
    ("B-20 B-21 B-22 B-23 N-17 N-18 B-24 N-20 N-21 N-22 N-23 N-24 Q-11 Q-14 Q-18", 3),
    ("B-5 B-6 B-7 N-2 N-3 N-4 N-6 N-7 N-8 N-9 N-10 N-12 N-13 N-14 N-15 N-16 N-19", 4),
    ("S-1 S-2 S-3 S-4 S-5 S-6 S-7 S-8 S-9", 4),
    ("P-1 P-2 P-3 P-4 P-5 P-6 P-7 P-8 P-9 P-10", 5),
    ("B-10 B-11 M-1 M-2 M-3 M-4 M-5 M-6 M-7 M-8 M-9", 6),
    ("B-12 B-14 T-1 T-2 T-3 Q-4 Q-5 Q-6 Q-7 Q-8 Q-9 Q-16 P-11", 7),
]:
    for i in ids.split():
        STEP[i] = step
# Problems introduced by the uncommitted threading work, not present upstream
NOT_UPSTREAM = {"M-3", "M-9", "P-5", "Q-2", "Q-12"}


def gh(*args, json_out=False, check=True):
    r = subprocess.run(["gh", *args], capture_output=True, text=True)
    if check and r.returncode != 0:
        sys.exit(f"gh {' '.join(args)}\n{r.stderr}")
    return json.loads(r.stdout) if json_out else r.stdout.strip()


def plain(md):
    """Markdown cell -> plain text for titles."""
    md = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", md)
    return md.replace("**", "").replace("`", "").strip()


def parse_rows():
    """Yield dicts {ID, Sev, ..., _section, _line} for every issue row."""
    lines = open(LOG).read().split("\n")
    section, header = None, None
    for n, line in enumerate(lines):
        if line.startswith("### "):
            section, header = line[4:].strip(), None
        if not line.startswith("|"):
            header = None if not line.strip() else header
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells[0] == "ID":
            header = cells
            continue
        if header is None or set(cells[0]) <= set("-: "):
            continue
        m = re.match(r"\[?([BNPMTQFS]-\d+)\]?", cells[0])
        if not m:
            continue
        row = dict(zip(header, cells))
        row.update(ID=m.group(1), _section=section, _line=n)
        yield row


def title_of(row):
    text = plain(row.get("Issue") or row.get("Item") or "")
    first = re.split(r"(?<=[a-z0-9)]):\s|;\s", text)[0]    # first clause
    if len(first) >= 30 and first.count("(") == first.count(")"):
        text = first
    if re.search(r"\d", row.get("3-D size", "")):
        text += f": {plain(row['3-D size'])} in a 256²×200 3-D run"
    if len(text) > 90:
        text = text[:87].rsplit(" ", 1)[0] + "…"
    return f"[{row['ID']}] {text}"


def labels_of(row):
    sev = row.get("Sev", "").lower()
    out = [AREA[row["ID"][0]]]
    if sev in SEV_COLOR:
        out.append(f"sev:{sev}")
    if row["ID"] in STEP:
        out.append(f"step-{STEP[row['ID']]}")
    if row["ID"] not in NOT_UPSTREAM and not row["ID"].startswith(("F-", "S-")):
        out.append("upstream")
    return out


def body_of(row):
    skip = {"ID", "Status", "_section", "_line"}
    fields = "\n".join(f"**{k}:** {v}" for k, v in row.items()
                       if k not in skip and v and v != "—")
    return (f"From the modernization log, section *{row['_section']}*.\n\n"
            f"{fields}\n\n"
            f"Plan step: {STEP.get(row['ID'], '—')}. "
            f"Full context, baseline results and reproducer cases: {LOG_URL}\n")


def existing_issues():
    found = {}
    for it in gh("issue", "list", "-R", REPO, "--state", "all", "--limit", "1000",
                 "--json", "number,title", json_out=True):
        m = re.match(r"\[([BNPMTQFS]-\d+)\]", it["title"])
        if m:
            found[m.group(1)] = it["number"]
    return found


def ensure_labels(dry):
    have = {l["name"] for l in gh("label", "list", "-R", REPO, "--limit", "200",
                                  "--json", "name", json_out=True)}
    want = {**{a: AREA_COLOR[a] for a in AREA.values()},
            **{f"sev:{s}": c for s, c in SEV_COLOR.items()},
            **{f"step-{i}": "bfdadc" for i in range(1, 8)},
            "upstream": "cccccc"}
    for name, color in want.items():
        if name not in have:
            print(f"label {name}")
            if not dry:
                gh("label", "create", name, "-R", REPO, "--color", color)


def link_ids(numbers):
    """Turn the ID cell of each row into a link to its issue."""
    text = open(LOG).read()
    def sub(m):
        i = m.group(1)
        if i not in numbers:
            return m.group(0)
        return f"| [{i}](https://github.com/{REPO}/issues/{numbers[i]}) |"
    new = re.sub(r"^\| ([BNPMTQFS]-\d+) \|", sub, text, flags=re.M)
    if new != text:
        open(LOG, "w").write(new)
        print("linked IDs in", LOG)


def cmd_issues(dry):
    rows = list(parse_rows())
    have = existing_issues()
    ensure_labels(dry)
    for row in rows:
        if row["ID"] in have:
            continue
        t, labels = title_of(row), labels_of(row)
        print(f"create {t}  [{', '.join(labels)}]")
        if dry:
            continue
        url = gh("issue", "create", "-R", REPO, "--title", t,
                 "--body", body_of(row), "--label", ",".join(labels))
        have[row["ID"]] = int(url.rstrip("/").rsplit("/", 1)[1])
        time.sleep(1.5)            # stay under GitHub's content-creation limits
    if not dry:
        link_ids(have)
    print(f"{len(rows)} IDs in the log, {len(have)} with issues")


def cmd_project(dry):
    projects = gh("project", "list", "--owner", OWNER, "--format", "json",
                  json_out=True)["projects"]
    proj = next((p for p in projects if p["title"] == PROJECT), None)
    if proj is None:
        print(f"create project {PROJECT}")
        if dry:
            return
        proj = gh("project", "create", "--owner", OWNER, "--title", PROJECT,
                  "--format", "json", json_out=True)
        gh("project", "link", str(proj["number"]), "--owner", OWNER, "--repo", REPO)
    num, pid = str(proj["number"]), proj["id"]

    want = {"Area": list(AREA.values()), "Severity": ["high", "medium", "low"],
            "Plan step": [str(i) for i in range(1, 8)]}
    fields = {f["name"]: f for f in gh("project", "field-list", num, "--owner", OWNER,
                                       "--format", "json", json_out=True)["fields"]}
    for name, opts in want.items():
        if name not in fields:
            print(f"field {name}")
            if not dry:
                gh("project", "field-create", num, "--owner", OWNER, "--name", name,
                   "--data-type", "SINGLE_SELECT", "--single-select-options", ",".join(opts))
    if dry:
        return
    fields = {f["name"]: f for f in gh("project", "field-list", num, "--owner", OWNER,
                                       "--format", "json", json_out=True)["fields"]}
    opt = {n: {o["name"]: o["id"] for o in fields[n]["options"]} for n in want}

    items = gh("project", "item-list", num, "--owner", OWNER, "--limit", "1000",
               "--format", "json", json_out=True)["items"]
    on_board = {it.get("content", {}).get("number") for it in items}
    rows = {r["ID"]: r for r in parse_rows()}
    for i, n in sorted(existing_issues().items(), key=lambda kv: kv[1]):
        if n in on_board:
            continue
        print(f"add #{n} [{i}]")
        item = gh("project", "item-add", num, "--owner", OWNER, "--url",
                  f"https://github.com/{REPO}/issues/{n}", "--format", "json", json_out=True)
        sev = rows.get(i, {}).get("Sev", "").lower()
        values = {"Area": AREA[i[0]], "Severity": sev if sev in opt["Severity"] else None,
                  "Plan step": str(STEP[i]) if i in STEP else None}
        for f, v in values.items():
            if v and v not in opt[f]:
                # gh cannot add options to an existing single-select field
                print(f"  no option '{v}' in field {f}: add it in the web UI, "
                      f"then set it on #{n} by hand")
            elif v:
                gh("project", "item-edit", "--id", item["id"], "--project-id", pid,
                   "--field-id", fields[f]["id"], "--single-select-option-id", opt[f][v])
    print(f"project {PROJECT}: https://github.com/users/{OWNER}/projects/{num}")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("issues", "project"):
        sys.exit(__doc__)
    {"issues": cmd_issues, "project": cmd_project}[sys.argv[1]]("--dry-run" in sys.argv)
