#!/usr/bin/env python3
"""Render docs/MODERNIZATION.md as a self-contained HTML page.

usage: make_log_page.py [--sync-md] [--no-github] [IN.md] [OUT.html]

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
status cells) and deletes them; run it on the integration branch now and then.

Status from GitHub: every log ID is an issue on dwijn/RH titled "[ID] ...",
and the issue state is the authority for the status column. Unless
--no-github is given, the states are fetched with `gh` and each status cell
is set to fixed (closed as completed), wontfix (closed as not planned) or
open (open issue; an "in progress" cell is kept). Cells that change are
listed. Without a working `gh` the cells are left as they are.
"""
import json
import subprocess
import glob
import html
import os
import re
import sys

import markdown

here = os.path.dirname(os.path.abspath(__file__))
flags = {"--sync-md", "--no-github"}
args = [a for a in sys.argv[1:] if a not in flags]
sync_md = "--sync-md" in sys.argv
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


def apply_fragments(text, frags):
    """Changelog entries for the fragments, and their IDs marked fixed."""
    entries = []
    for _, m, body in frags:
        ids = ", ".join(i for i in m["ids"] if i)
        heading = f"### {m.get('date', '')} — {m.get('title', ids)}" + (f" ({ids})" if ids else "")
        if heading not in text:
            entries.append(f"{heading}\n\n{body}\n")
        for i in m["ids"]:
            text = re.sub(rf"^(\| \[?{re.escape(i)}(?![0-9])\]?[^\n]*\| )open( \|)$",
                          r"\1fixed\2", text, flags=re.M)
    if entries:
        text = text.replace("## 5. Changelog\n\n",
                            "## 5. Changelog\n\n" + "\n".join(entries) + "\n", 1)
    return text


def github_states(repo="dwijn/RH"):
    """{ID: status} from the issue states, or None if gh is unavailable."""
    try:
        out = subprocess.run(
            ["gh", "issue", "list", "-R", repo, "--state", "all", "--limit", "1000",
             "--json", "title,state,stateReason"],
            capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError) as e:
        print(f"warning: no issue states from GitHub ({e}); status cells unchanged",
              file=sys.stderr)
        return None
    states = {}
    for i in json.loads(out):
        m = re.match(r"\[([A-Z]+-[0-9]+)\]", i["title"])
        if not m:
            continue
        if i["state"] == "OPEN":
            states[m.group(1)] = "open"
        else:
            states[m.group(1)] = ("wontfix" if i.get("stateReason") == "NOT_PLANNED"
                                  else "fixed")
    return states


def apply_github(text, states):
    """Set the status cell of every issue-table row from its issue state."""
    def fix(m):
        id_, old = m.group(2), m.group(4)
        new = states.get(id_)
        if new is None or new == old or (new == "open" and old == "in progress"):
            return m.group(0)
        print(f"status {id_}: {old} -> {new}", file=sys.stderr)
        return f"{m.group(1)}{m.group(3)}{new} |"
    return re.sub(r"^(\| \[?([A-Z]+-[0-9]+)\]?)([^\n]*\| )"
                  r"(open|fixed|wontfix|in progress) \|$", fix, text, flags=re.M)


frags = read_fragments()
text = open(src).read()
if use_github:
    states = github_states()
    if states is not None:
        text = apply_github(text, states)
text = apply_fragments(text, frags)
if sync_md:
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
