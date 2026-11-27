#!/usr/bin/env python3
"""Render docs/MODERNIZATION.md as a self-contained HTML page.

usage: make_log_page.py [--sync-md | --check] [--no-github] [IN.md] [OUT.html]

The Markdown file is the source of truth; the HTML is what gets published as
the shared "RH modernization log" page. Severity (high/medium/low) and status
(open/fixed/wontfix/in progress) cells in the issue tables become chips.

Changelog fragments: each change adds docs/changes/<ID>.md instead of editing
MODERNIZATION.md, so that pull requests do not conflict. A fragment starts
with "key: value" header lines and a blank line, then the entry text:

    ids: B-2            IDs the change fixes (comma- or space-separated;
                        empty for a change that closes no issue)
    date: 2026-10-15
    title: One-line summary

The page shows every fragment at the top of the changelog and marks its IDs
fixed. --sync-md folds the fragments into MODERNIZATION.md (changelog entry,
status cells) and deletes them. Code pull requests do not edit
MODERNIZATION.md at all: the row is logged beforehand (docs pull request),
the code pull request adds the fragment and closes the issue, and after
merging a docs pull request runs --sync-md. --check writes nothing and exits
1 when a sync is due (fragments to fold, status cells that GitHub changes,
or At a glance counts that differ).

A fragment does not mark an ID fixed while its issue is open on GitHub: the
fragment then fixes a part of a multi-part item, or the issue has not been
closed yet; the row stays open and a warning names it.

At a glance: the sentence "N issues are logged; N are fixed and N closed
without a fix. Of the N open, N are on the wishlist." is recounted from the
issue tables of section 3 (3.1 to 3.8; 3.9 repeats IDs), the wishlist being
the open IDs whose plan step in make_issues.py is Wishlist; so is "N of the
problems found also exist upstream (§3.9), and N of them are fixed here",
from the IDs of 3.9 and their status. A status cell
other than open, fixed or wontfix ("superseded by B-64") counts as the state
of its issue. --sync-md also sets the date of "*Status on YYYY-MM-DD.*" when
it changes anything. The rest of At a glance is prose, updated by hand.

Status from GitHub: every log ID is an issue on dwijn/RH titled "[ID] ...",
and the issue state is the authority for the status column. Unless
--no-github is given, the states are fetched with `gh` and each status cell
is set to fixed (closed as completed), wontfix (closed as not planned) or
open (open issue; an "in progress" cell is kept). A cell that becomes fixed
names the merged pull requests that closed the issue, "fixed (#455)"; a cell
already fixed, with or without numbers, is kept. GitHub links an issue to the
pull request that closed it only for pull requests into the default branch;
for the library track's branch (chrysalis, docs/LIBRARY.md §6.1), whose
issues are closed by hand, the numbers come from the merged pull requests
into that branch that say "Closes #N". In the table of bugs also
present upstream (a "Fixed here" column), an empty cell ("—") of a closed
issue gets links to those pull requests. Cells that change are listed.
Without a working `gh` the cells are left as they are.
"""
import datetime
import json
import subprocess
import glob
import html
import os
import re
import sys

import markdown

here = os.path.dirname(os.path.abspath(__file__))
flags = {"--sync-md", "--check", "--no-github"}
args = [a for a in sys.argv[1:] if a not in flags]
sync_md = "--sync-md" in sys.argv
check = "--check" in sys.argv
if sync_md and check:
    sys.exit("make_log_page.py: --sync-md and --check exclude each other")
use_github = "--no-github" not in sys.argv
src = args[0] if len(args) > 0 else os.path.join(here, "MODERNIZATION.md")
dst = args[1] if len(args) > 1 else os.path.join(here, "MODERNIZATION.html")
frag_dir = os.path.join(os.path.dirname(os.path.abspath(src)), "changes")


def read_fragments():
    """[(path, {ids, date, title}, body)], newest first."""
    out = []
    for f in sorted(glob.glob(os.path.join(frag_dir, "*.md"))):
        head, _, body = open(f).read().partition("\n\n")
        meta = dict(l.split(":", 1) for l in head.splitlines() if ":" in l)
        meta = {k.strip(): v.strip() for k, v in meta.items()}
        meta["ids"] = [i for i in re.split(r"[,\s]+", meta.get("ids", "").strip())
                       if i]       # no ids: a change that closes no issue
        out.append((f, meta, body.strip()))
    return sorted(out, key=lambda x: x[1].get("date", ""), reverse=True)


def apply_fragments(text, frags, states=None):
    """Changelog entries for the fragments, and their IDs marked fixed. An ID
    whose issue is still open on GitHub stays open: the fragment fixes a part
    of it (a multi-part item), or the issue has yet to be closed."""
    entries = []
    for path, m, body in frags:
        ids = ", ".join(i for i in m["ids"] if i)
        heading = f"### {m.get('date', '')} — {m.get('title', ids)}" + (f" ({ids})" if ids else "")
        if heading not in text:
            entries.append(f"{heading}\n\n{body}\n")
        for i in m["ids"]:
            if states and states.get(i, (None,))[0] == "open":
                print(f"warning: {os.path.basename(path)} names {i}, whose issue is open: "
                      f"its row stays open (a part of a multi-part item, or close the issue "
                      f"and sync again)", file=sys.stderr)
                continue
            text = re.sub(rf"^(\| \[?{re.escape(i)}(?![0-9])\]?[^\n]*\| )open( \|)$",
                          r"\1fixed\2", text, flags=re.M)
    if entries:
        text = text.replace("## 5. Changelog\n\n",
                            "## 5. Changelog\n\n" + "\n".join(entries) + "\n", 1)
    return text


# The integration branch of the library track (docs/LIBRARY.md §6.1)
SIDE_BRANCH = "chrysalis"


def github_states(repo="dwijn/RH"):
    """{ID: (status, [merged PR numbers that closed it])} from the issue
    states, or None if gh is unavailable."""
    def gh_json(*a):
        return json.loads(subprocess.run(["gh", *a, "-R", repo, "--limit", "1000"],
                                         capture_output=True, text=True, check=True).stdout)
    try:
        issues = gh_json("issue", "list", "--state", "all", "--json",
                         "number,title,state,stateReason,closedByPullRequestsReferences")
        merged = {p["number"] for p in gh_json("pr", "list", "--state", "merged",
                                                "--json", "number")}
        # Pull requests merged into the library track's branch: GitHub does
        # not close or link issues for them, so read their "Closes #N"
        side = {}
        for p in gh_json("pr", "list", "--state", "merged", "--base", SIDE_BRANCH,
                         "--json", "number,body"):
            for n in re.findall(r"(?i)\b(?:close[sd]?|fix(?:e[sd])?|resolve[sd]?)\s+#(\d+)",
                                p.get("body") or ""):
                side.setdefault(int(n), set()).add(p["number"])
    except (OSError, subprocess.CalledProcessError) as e:
        print(f"warning: no issue states from GitHub ({e}); status cells unchanged",
              file=sys.stderr)
        return None
    states = {}
    for i in issues:
        m = re.match(r"\[([A-Z]+-[0-9]+)\]", i["title"])
        if not m:
            continue
        if i["state"] == "OPEN":
            states[m.group(1)] = ("open", [])
        elif i.get("stateReason") == "NOT_PLANNED":
            states[m.group(1)] = ("wontfix", [])
        else:
            prs = sorted(p["number"] for p in i.get("closedByPullRequestsReferences") or []
                         if p["number"] in merged
                         and p.get("repository", {}).get("name", "RH") == repo.split("/")[1])
            if not prs:
                prs = sorted(side.get(i.get("number"), ()))
            states[m.group(1)] = ("fixed", prs)
    return states


def cells(line):
    """Split a table row at the pipes that are not escaped (\\|)."""
    return re.split(r"(?<!\\)\|", line)


def apply_github(text, states):
    """Set the status cell of every issue-table row from its issue state, and
    fill the empty "Fixed here" cells of the upstream table."""
    def fix(m):
        id_, old = m.group(2), m.group(4)
        new, prs = states.get(id_, (None, []))
        kind = old.split(" (")[0]
        bare_fixed = new == kind == "fixed" and old == "fixed" and prs
        if new is None or (new == kind and not bare_fixed) or \
           (new == "open" and kind == "in progress"):
            return m.group(0)
        if new == "fixed" and prs:
            new = "fixed (" + ", ".join(f"#{n}" for n in prs) + ")"
        print(f"status {id_}: {old} -> {new}", file=sys.stderr)
        return f"{m.group(1)}{m.group(3)}{new} |"
    text = re.sub(r"^(\| \[?([A-Z]+-[0-9]+)\]?)([^\n]*\| )"
                  r"(open|fixed(?: \([^)|\n]*\))?|wontfix|in progress) \|$", fix, text,
                  flags=re.M)

    lines, col = text.split("\n"), None
    for n, line in enumerate(lines):
        if not line.startswith("|"):
            col = None
            continue
        c = cells(line)
        names = [x.strip() for x in c]
        if names[1] == "ID":
            col = names.index("Fixed here") if "Fixed here" in names else None
            continue
        m = re.match(r"\[?([A-Z]+-[0-9]+)\]?", names[1])
        if col is None or not m or names[col] != "—":
            continue
        new, prs = states.get(m.group(1), (None, []))
        if new != "fixed" or not prs:
            continue
        c[col] = " " + ", ".join(f"[#{p}](https://github.com/dwijn/RH/pull/{p})"
                                 for p in prs) + " "
        print(f"fixed here {m.group(1)}: — -> " + ", ".join(f"#{p}" for p in prs),
              file=sys.stderr)
        lines[n] = "|".join(c)
    return "\n".join(lines)


def wishlist_ids():
    """The IDs whose plan step is the wishlist, from make_issues.py (the
    authority for the plan step of every ID), or None."""
    try:
        sys.path.insert(0, here)
        sys.dont_write_bytecode = True  # no docs/__pycache__
        import make_issues
        return {i for i, s in make_issues.STEP.items() if s == "Wishlist"}
    except Exception as e:  # noqa: BLE001 - any failure leaves the counts alone
        print(f"warning: no plan steps from make_issues.py ({e})", file=sys.stderr)
        return None


def tally(text, states):
    """Counts for "At a glance" from the issue tables of section 3 (3.1 to
    3.8; 3.9 repeats IDs): {logged, fixed, wontfix, open, wishlist,
    upstream, upstream_fixed}, or None when a status cannot be classified.
    A status cell is fixed ("fixed", "fixed (#n)"), wontfix, or open ("open",
    "open — ...", "in progress"); any other cell ("superseded by B-64") takes
    the state of its issue. upstream counts the IDs of 3.9, upstream_fixed
    those of them that are fixed (by the same classification); both None
    when 3.9 names an ID that 3.1 to 3.8 do not."""
    m = re.search(r"^## 3\. Issues$(.*?)^### 3\.9 ", text, flags=re.M | re.S)
    wish = wishlist_ids()
    if not m or wish is None:
        print("warning: At a glance counts not computed (section 3 or plan steps not found)",
              file=sys.stderr)
        return None
    kind, unknown = {}, []
    for line in m.group(1).split("\n"):
        if not line.startswith("| "):
            continue
        c = cells(line)
        mm = re.match(r"\[?([A-Z]+-[0-9]+)\]?", c[1].strip())
        if not mm or mm.group(1) in kind:
            continue
        id_, s = mm.group(1), c[-2].strip()
        if s.startswith("fixed"):
            k = "fixed"
        elif s == "wontfix":
            k = "wontfix"
        elif s.startswith("open") or s == "in progress":
            k = "open"
        else:
            k = (states or {}).get(id_, (None,))[0]
            if k is None:
                unknown.append(f"{id_} ({s})")
        kind[id_] = k
    if unknown:
        print("warning: At a glance counts not computed: status unknown for " +
              ", ".join(unknown), file=sys.stderr)
        return None
    n = {k: sum(1 for v in kind.values() if v == k) for k in ("fixed", "wontfix", "open")}
    n["logged"] = len(kind)
    n["wishlist"] = sum(1 for i, v in kind.items() if v == "open" and i in wish)
    up = []
    m9 = re.search(r"^### 3\.9 .*?(?=^## |\Z)", text, flags=re.M | re.S)
    for line in (m9.group(0) if m9 else "").split("\n"):
        if line.startswith("| "):
            mm = re.match(r"\[?([A-Z]+-[0-9]+)\]?", cells(line)[1].strip())
            if mm and mm.group(1) not in up:
                up.append(mm.group(1))
    missing = [i for i in up if i not in kind]
    if missing or not up:
        print("warning: upstream counts not computed: " +
              ("3.9 not found" if not up else "not in 3.1 to 3.8: " + ", ".join(missing)),
              file=sys.stderr)
        n["upstream"] = n["upstream_fixed"] = None
    else:
        n["upstream"] = len(up)
        n["upstream_fixed"] = sum(1 for i in up if kind[i] == "fixed")
    return n


GLANCE = re.compile(r"(\d+)( issues are logged;\s+)(\d+)( are fixed and\s+)(\d+)"
                    r"(\s+closed without a fix\.\s+Of the\s+)(\d+)( open,\s+)(\d+)"
                    r"( are on the wishlist\.)")
GLANCE_KEYS = ("logged", "fixed", "wontfix", "open", "wishlist")
UPSTREAM = re.compile(r"(\d+)( of the problems found also exist upstream \(§3\.9\), and\s+)"
                      r"(\d+)( of them\s+are fixed here\.)")
UPSTREAM_KEYS = ("upstream", "upstream_fixed")


def apply_counts(text, n):
    """Write the counts into the two "At a glance" sentences that have them;
    list what changed."""
    if n is None:
        return text
    for regex, keys, what in ((GLANCE, GLANCE_KEYS, "the issue counts"),
                              (UPSTREAM, UPSTREAM_KEYS, "the upstream counts")):
        m = regex.search(text)
        if n[keys[0]] is None:
            continue
        if not m:
            print(f"warning: the At a glance sentence with {what} was not found",
                  file=sys.stderr)
            continue
        g = list(m.groups())
        for k, i in zip(keys, range(0, 2 * len(keys), 2)):
            if int(g[i]) != n[k]:
                print(f"at a glance {k}: {g[i]} -> {n[k]}", file=sys.stderr)
            g[i] = str(n[k])
        text = text[:m.start()] + "".join(g) + text[m.end():]
    return text


frags = read_fragments()
text = orig = open(src).read()
states = None
if use_github:
    states = github_states()
    if states is not None:
        text = apply_github(text, states)
github_changed = text != orig
text = apply_fragments(text, frags, states)
before_counts = text
text = apply_counts(text, tally(text, states))
counts_changed = text != before_counts
if check:
    due = [f"{len(frags)} fragment(s) to fold"] if frags else []
    if github_changed:
        due.append("cells to update from GitHub (listed above)")
    if counts_changed:
        due.append("At a glance counts (listed above)")
    print("log sync due: " + "; ".join(due) if due else "log in sync")
    sys.exit(1 if due else 0)
if sync_md:
    if text != orig:
        text = re.sub(r"\*Status on \d{4}-\d{2}-\d{2}\.\*",
                      f"*Status on {datetime.date.today().isoformat()}.*", text, count=1)
    open(src, "w").write(text)
    for f, _, _ in frags:
        os.remove(f)
    print(f"folded {len(frags)} fragment(s) into {src}")
title_md, body_md = text.split("\n", 1)
body = markdown.markdown(body_md, extensions=["tables", "sane_lists", "toc"],
                         extension_configs={"toc": {"permalink": False}})

# Chips for severity and status cells
chip = {
    "high": "sev-high", "medium": "sev-med", "low": "sev-low",
    "open": "st-open", "fixed": "st-fixed", "wontfix": "st-wont",
    "in progress": "st-prog",
}
def chipify(m):
    word = m.group(1)
    key = word.lower()
    cls = chip.get(key) or (chip["open"] if key.startswith("open") else None)
    if not cls:
        return m.group(0)
    return f'<td><span class="chip {cls}">{word}</span></td>'
body = re.sub(r"<td>(high|medium|low|open|fixed|wontfix|in progress)</td>", chipify, body)
body = re.sub(r"<td>fixed (\([^<]*\))</td>",
              lambda m: f'<td><span class="chip st-fixed">fixed</span> {m.group(1)}</td>', body)
body = re.sub(r"<td>(open — [^<]*)</td>",
              lambda m: f'<td><span class="chip st-open">open</span> {m.group(1)[7:]}</td>', body)

# Wrap tables so wide ones scroll inside their own box
body = body.replace("<table>", '<div class="tbl"><table>').replace("</table>", "</table></div>")
# Issue IDs as anchors-to-be: make B-1 etc. monospace in the first column
body = re.sub(r'<td>((?:<a href="[^"]+">)?[BNPMTQ]-\d+(?:</a>)?)</td>', r'<td class="id">\1</td>', body)

# Section navigation from h2 headings
navs = re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', body)
nav = "".join(f'<a href="#{i}">{re.sub("<[^>]+>", "", t)}</a>' for i, t in navs)

page = f"""<title>RH Modernization Log</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,600;12..96,700&family=JetBrains+Mono:wght@400;500&family=Source+Sans+3:ital,wght@0,400;0,600;1,400&display=swap">
<style>
/* Layout: an engineering logbook that pairs with the "RH code map" page.
   Same faces and teal/amber family; one reading column, wide tables scroll in
   their own boxes; severity and status as chips so open problems stand out. */
:root {{
  --bg: #eef2f5; --surface: #ffffff; --fg: #13202b; --muted: #4b5b69; --line: #c3ced7;
  --accent: #0b7a86; --accent-bg: #d8eff1; --warm: #a84a0a; --warm-bg: #fae3d0;
  --bad: #b3261e; --bad-bg: #fbe1de; --good: #2d6a2f; --good-bg: #dcefd9; --link: #0b5fa5;
  --f-display: "Bricolage Grotesque", "Segoe UI", system-ui, sans-serif;
  --f-body: "Source Sans 3", "Segoe UI", system-ui, sans-serif;
  --f-mono: "JetBrains Mono", ui-monospace, "SF Mono", Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --bg: #0d1419; --surface: #141d24; --fg: #e4ecf2; --muted: #9db0bf; --line: #2c3b47;
    --accent: #4cc3cf; --accent-bg: #10353a; --warm: #f2a766; --warm-bg: #3b2814;
    --bad: #f28b82; --bad-bg: #3d1a17; --good: #8fd18b; --good-bg: #183019; --link: #7ab8f0;
    color-scheme: dark;
  }}
}}
:root[data-theme="dark"] {{
  --bg: #0d1419; --surface: #141d24; --fg: #e4ecf2; --muted: #9db0bf; --line: #2c3b47;
  --accent: #4cc3cf; --accent-bg: #10353a; --warm: #f2a766; --warm-bg: #3b2814;
  --bad: #f28b82; --bad-bg: #3d1a17; --good: #8fd18b; --good-bg: #183019; --link: #7ab8f0;
  color-scheme: dark;
}}
* {{ box-sizing: border-box; }}
body {{ background: var(--bg); color: var(--fg); font: 400 16px/1.55 var(--f-body); }}
.page {{ max-width: 1080px; margin: 0 auto; padding-block: 28px 72px;
  padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px)); display: grid; gap: 14px; }}
.page > * {{ min-width: 0; }}
h1, h2, h3 {{ font-family: var(--f-display); text-wrap: balance; line-height: 1.15; margin: 0; }}
h1 {{ font-size: 2.4rem; font-weight: 700; letter-spacing: -0.01em; }}
h2 {{ font-size: 1.6rem; font-weight: 700; margin-top: 30px; padding-top: 14px; border-top: 2px solid var(--line); scroll-margin-top: 12px; }}
h3 {{ font-size: 1.12rem; font-weight: 600; margin-top: 14px; }}
p, li {{ max-width: 76ch; }}
p, ul, ol {{ margin: 0; }}
ul, ol {{ padding-left: 1.2em; display: grid; gap: 5px; }}
a {{ color: var(--link); }}
code {{ font-family: var(--f-mono); font-size: 0.84em; background: var(--accent-bg); padding: 0.05em 0.3em; border-radius: 4px; overflow-wrap: anywhere; }}
strong {{ font-weight: 600; }}
hr {{ display: none; }}
nav {{ display: flex; flex-wrap: wrap; gap: 6px 18px; padding-block: 10px; border-block: 1px solid var(--line); font-size: 0.95rem; }}
.sub {{ color: var(--muted); font-size: 1.05rem; }}
#at-a-glance + ul {{ background: var(--surface); border: 1px solid var(--line); border-left: 4px solid var(--accent);
  border-radius: 8px; padding: 14px 18px 14px 34px; }}
.tbl {{ overflow-x: auto; background: var(--surface); border: 1px solid var(--line); border-radius: 8px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.92rem; font-variant-numeric: tabular-nums; }}
th, td {{ text-align: left; vertical-align: top; padding: 8px 11px; border-bottom: 1px solid var(--line); }}
th {{ font-size: 0.78rem; letter-spacing: 0.05em; text-transform: uppercase; color: var(--muted); font-weight: 600; white-space: nowrap; }}
tr:last-child td {{ border-bottom: 0; }}
td.id {{ font-family: var(--f-mono); font-size: 0.85em; white-space: nowrap; font-weight: 500; }}
td code {{ white-space: normal; }}
.chip {{ display: inline-block; font-size: 0.75rem; font-weight: 600; letter-spacing: 0.03em; padding: 1px 8px;
  border-radius: 999px; border: 1.5px solid; white-space: nowrap; }}
.sev-high {{ color: var(--bad); background: var(--bad-bg); border-color: var(--bad); }}
.sev-med  {{ color: var(--warm); background: var(--warm-bg); border-color: var(--warm); }}
.sev-low  {{ color: var(--muted); background: transparent; border-color: var(--line); }}
.st-open  {{ color: var(--muted); background: transparent; border-color: var(--muted); }}
.st-prog  {{ color: var(--accent); background: var(--accent-bg); border-color: var(--accent); }}
.st-fixed {{ color: var(--good); background: var(--good-bg); border-color: var(--good); }}
.st-wont  {{ color: var(--muted); background: transparent; border-color: var(--line); text-decoration: line-through; }}
:focus-visible {{ outline: 2px solid var(--link); outline-offset: 2px; }}
</style>
<div class="page">
<header style="display:grid;gap:10px">
<h1>{html.escape(title_md.lstrip('# ').strip())}</h1>
<p class="sub">Issues found in RH, baseline measurements against upstream, and every change made, with how it was verified. Generated from <code>docs/MODERNIZATION.md</code>.</p>
<nav aria-label="Sections">{nav}</nav>
</header>
{body}
</div>
"""
open(dst, "w").write(page)
print(f"wrote {dst}")
