# Causal Synergy Curve Reproducibility Archive

This repository contains the reproducibility materials for the Statistics in Medicine submission:

**The Causal Synergy Curve: Bayesian Estimation of Gap-Closure Architecture Across Modifiable Domains**

The archive contains the Stan models, R scripts, analytic NHANES illustration file, codebook, and run notes used to reproduce the simulation and applied analyses reported in the manuscript.

## Contents

```text
code/
  run_fullbayes_simulation.R
  merge_fullbayes_simulation_chunks.R
  make_fullbayes_simulation_figures.R
  run_fullbayes_nhanes_applied.R
  make_fullbayes_applied_tables_figures.R
  csc_design_readouts.R
  compare_vb_mcmc_agreement.R
  run_vb_simulation.R
  merge_vb_simulation_chunks.R

stan/
  csc_coalition_horseshoe_binary.stan
  csc_coalition_horseshoe_weibull.stan

data/
  nhanes_8sdoh_analytic.csv
  nhanes_8sdoh_codebook.md
```

## Requirements

The production analyses were run in R with CmdStan through `cmdstanr`.

Required R packages include:

- `cmdstanr`
- `posterior`
- `data.table`
- `dplyr`
- `tidyr`
- `purrr`
- `ggplot2`
- `splines`
- `survival`

Install CmdStan through `cmdstanr` before running the full Bayesian scripts.

## Recommended Run Order

### Applied NHANES analysis

```r
Rscript code/run_fullbayes_nhanes_applied.R
Rscript code/make_fullbayes_applied_tables_figures.R
```

### Full Bayesian simulation study

```r
Rscript code/run_fullbayes_simulation.R
Rscript code/merge_fullbayes_simulation_chunks.R
Rscript code/make_fullbayes_simulation_figures.R
```

The full simulation grid is computationally intensive. Runtime depends on the number of chains, iterations, and available CPU cores.

## Data Notes

The applied illustration uses an analysis-ready public-use NHANES 8-SDOH file included in `data/`. The file is provided for manuscript reproducibility. Source NHANES public-use files are available from the National Center for Health Statistics. Users should follow NCHS/CDC terms and documentation when reconstructing or extending the analysis.

## Files Omitted Intentionally

The repository excludes compiled Stan executables, CmdStan draw CSV files, intermediate logs, and platform-specific worker scripts. These are transient artifacts that can be regenerated from the supplied R and Stan source files.

## Citation

If using these materials, please cite the archived release DOI once minted by Zenodo. Citation metadata are provided in `CITATION.cff` and `.zenodo.json`.
