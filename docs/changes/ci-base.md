ids: Q-14
date: 2026-10-17
title: CI compares each change with its base commit; drift from upstream is reported

- The regression job builds the base commit of the change (the target
  branch before it) and runs `run_tests.sh` against it in tolerance mode.
  Results must match, unless the pull request carries the label
  `changes-results` (an intended change of results): then the differences
  from the base are reported in the job summary instead of failing, while
  sanitizer, thread-count and memory checks still fail.
- The comparison with upstream RH (c95bdf6) remains, as a report that never
  fails. Upstream is the reference we trust but check: every difference
  from it needs an explanation in this log (a known numerical change or a
  fixed bug). The analytic unit tests guard against slow drift.
- Why: the first fixes that change results on purpose (Bezier weights,
  N-17: up to 7e-7) would otherwise fail every later comparison with
  upstream, or force a tolerance too loose to catch accidents.
