ids:
date: 2026-11-19
title: CI fails a code pull request that edits the log; status cells come from the log sync

Code pull requests no longer edit `docs/MODERNIZATION.md`. Each one marked its own row fixed, and those edits conflicted with each other (new rows sit next to each other at the end of a section, and git does not merge edits to adjacent lines) and with docs pull requests that rewrote the same rows; resolving the conflicts of #445, #455, #456 and #458 lost five status cells (restored in #462). A code pull request now adds its change note `docs/changes/<ID>.md` and closes its issue, and a docs pull request runs `docs/make_log_page.py --sync-md` after the merge: it folds the notes into the changelog and writes the status cells from GitHub, as `fixed (#PR)` with the merged pull requests that closed the issue (#464). The new CI job `log` fails a code pull request that changes `docs/MODERNIZATION.md`.
