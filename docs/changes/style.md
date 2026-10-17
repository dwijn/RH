ids:
date: 2026-10-16
title: One code style, enforced by clang-format

The code followed Emacs's C mode (2-space indentation with leading runs
of 8 spaces stored as tabs), which looked broken in any editor not set
to 8-column tabs. All C sources are now formatted by clang-format
(`.clang-format`): K&R braces, one tab per level shown as 4 columns,
spaces for alignment, aligned declarations and assignments, no line
limit yet (existing line breaks kept). `.editorconfig` sets tab width 4
for editors and GitHub, a CI job rejects unformatted code, and
`.git-blame-ignore-revs` lets `git blame` skip the reformat commit. The
reformat changes no machine code (all 27 programs compile identically).
With 4-column tabs, 960 lines exceed 80 columns (56 before), mostly in
deeply nested loops; they are left for the structural work (step 7).
Function and variable naming remains inconsistent (camelCase,
Capitalized, snake_case); also for step 7.
