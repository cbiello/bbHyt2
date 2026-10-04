---
name: bbhytplots-repo
description: Where bbH yt^2 M1M0 plots go - bbhytplots repo, oneloop dir, mirrors tthplots/oneloop with GMVFNS-style stat error bars
metadata:
  type: reference
---

Plots for the bbHyt2 M1M0 study go in git@github.com:cbiello/bbhytplots.git (clone at /scratch/cbiello/plots-repos/bbhytplots), subdir `oneloop/`, same organisation as git@github.com:cbiello/tthplots.git `oneloop/` (topfiles converters -> datfiles/<tag>-run/distributions -> myscprit script, MATRIX-style gnuplot).
Statistical errors must be error bars at bin centres as in the bbH GMVFNS project on this host (yueh): bbh_GMFVNS_main/integrator/runs-alphas2/*/plots/make_ebars.py + extra_main/extra_ratio hooks in myscript-eulerrun.py.
Chain: `oneloop/make_all.sh <TAG> <exact dir> <mass dir>` (tags HT4, HT4lam, mH2; RUNS env var for the run location). The same oneloop/ is also committed inside bbHyt2 as `bbhyt2plots/oneloop/` (user asked for it, 2026-10-04); keep both in sync.
Related: [[m1m0-massification-status]].
