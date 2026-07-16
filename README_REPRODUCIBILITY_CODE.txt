Causal Synergy Curve reproducibility archive
2026-07-16

This archive accompanies the Statistics in Medicine submission:
"The Causal Synergy Curve: A Bayesian Closure-Capacity Framework for
Intervention Architecture."

Contents
--------
code/
  R scripts used to run the full Bayesian simulation study, fit the NHANES
  applied model, merge simulation chunks, and generate manuscript tables
  and figures.

stan/
  Stan source files for the regularized coalition-surface models. Compiled
  Windows executables are intentionally omitted.

data/
  nhanes_8sdoh_analytic.csv: analysis-ready public-use NHANES 8-SDOH
  illustration file used by the applied scripts.
  nhanes_8sdoh_codebook.md: variable definitions for the analytic file.

*.md
  Run logs and validation notes retained for provenance.

Recommended run order
---------------------
1. Install R packages used by the scripts, including cmdstanr, posterior,
   data.table, dplyr, tidyr, purrr, ggplot2, splines, and survival.
2. Install CmdStan through cmdstanr if it is not already installed.
3. For the applied NHANES analysis:
   Rscript code/run_fullbayes_nhanes_applied.R
   Rscript code/make_fullbayes_applied_tables_figures.R
4. For the full Bayesian simulation study:
   Rscript code/run_fullbayes_simulation.R
   Rscript code/merge_fullbayes_simulation_chunks.R
   Rscript code/make_fullbayes_simulation_figures.R

Notes
-----
- The production simulation is computationally intensive. The manuscript
  reports the completed full Bayesian run; rerunning the full grid can take
  substantial wall time depending on chains, iterations, and CPU resources.
- Intermediate posterior draw CSVs, logs, worker scripts, and compiled Stan
  executables are not included because they are platform-specific or large
  transient artifacts. They can be regenerated from the R and Stan source.
- The manuscript submission archive is intentionally manuscript-only. This
  reproducibility archive is the code/data companion.
