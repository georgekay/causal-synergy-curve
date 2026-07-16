# SiM revision: fully Bayesian Causal Synergy Curve (2026-07-11)

This package updates the Causal Synergy Curve estimator so the simulation and applied example use
the same fully Bayesian Stan MCMC workflow.

## What changed vs `ready_upload_SiM_EIntegrated_20260703`

1. One estimator everywhere: a regularized-horseshoe coalition surface. The outcome model is an
   ordinary GLM, Cox, or hazard model with a structured linear predictor: covariates, main effects,
   burden-depth, anchored gates, and a coalition block under a regularized horseshoe prior.
2. The simulation harness is fully Bayesian Stan MCMC, not a Gaussian or Laplace posterior
   approximation. It uses the same posterior g-computation readouts as the application.
3. The simulation reports posterior regime-probability calibration, recovery, and credible-interval
   coverage.
4. The applied example has been refit with the same horseshoe coalition surface using a Weibull PH
   model, full MCMC, and weighted 10-year g-computation.

## Files

- `stan/csc_coalition_horseshoe_binary.stan`: cloglog binary model for simulation.
- `stan/csc_coalition_horseshoe_weibull.stan`: Weibull PH model for applied survival analysis.
- `code/csc_design_readouts.R`: shared design construction and posterior g-computation of closure,
  SC(k), ESR, Shapley values, gate B/E/G, and second differences.
- `code/run_fullbayes_simulation.R`: ADEMP simulation harness; writes recovery, coverage, and
  calibration CSV and LaTeX files.
- `code/run_fullbayes_nhanes_applied.R`: applied refit; writes architecture, gate, pairwise,
  summary, and diagnostic CSVs.
- `code/make_fullbayes_applied_tables_figures.R`: builds applied paper tables and figures from the
  applied MCMC outputs.
- `manuscript/methods_estimating_model_revised.tex`: drop-in estimating-model text.
- `manuscript/simulation_fullbayes_revised.tex`: drop-in simulation-study text.

## How to run locally

Run from the `code/` directory.

```powershell
cd D:\AMC_Papers\Causal_Synergy_Curve_Methods_Paper\ready_upload_SiM_FullBayes_20260711\code

# Production simulation, single serial R process. This is reproducible but slow.
$env:CSC_REPS = "500"
$env:CSC_N = "8000"
$env:CSC_CHAINS = "4"
$env:CSC_ITER_WARMUP = "500"
$env:CSC_ITER_SAMPLING = "500"
Rscript run_fullbayes_simulation.R

# Production simulation, machine-scale launch. On the 80-logical-processor Xeon machine this
# defaults to 20 concurrent workers x 4 chains per fit = about 80 chains in flight.
.\launch_fullbayes_simulation_max.ps1 -Reps 500 -N 8000 -ChainsPerFit 4

# After the launched worker processes finish, merge the chunk files.
# Use the CSC_OUTDIR path printed by launch_fullbayes_simulation_max.ps1.
$env:CSC_OUTDIR = "D:\AMC_Papers\Causal_Synergy_Curve_Methods_Paper\ready_upload_SiM_FullBayes_20260711\simulation_outputs\production_fullbayes_YYYYMMDD_HHMMSS"
Rscript merge_fullbayes_simulation_chunks.R

# Applied cohorts
$env:CSC_NHANES = "..\data\nhanes_8sdoh_analytic.csv"
$env:CSC_COHORT = "pooled"
Rscript run_fullbayes_nhanes_applied.R
$env:CSC_COHORT = "black"
Rscript run_fullbayes_nhanes_applied.R
$env:CSC_COHORT = "white"
Rscript run_fullbayes_nhanes_applied.R

# Applied tables and figures
Rscript make_fullbayes_applied_tables_figures.R
```

## Run status after 2026-07-14 audit

The NHANES applied analysis has been run with full Stan MCMC for pooled, Black, and White cohorts:

- 4 chains;
- 1000 warmup iterations per chain;
- 1000 post-warmup iterations per chain;
- adapt delta 0.99;
- maximum tree depth 12;
- full-cohort posterior g-computation readouts using 250 retained posterior draws.

CmdStan CSVs are preserved under `applied_outputs/cmdstan_csv_*_20260714/`. Summarized diagnostics,
architecture/gate tables, and figures are in `applied_outputs/`, `tables/`, and `figures/`.

The top-level simulation outputs required by `manuscript/simulation_fullbayes_revised.tex`
(`simulation_outputs/fullbayes_recovery_coverage.tex` and
`simulation_outputs/fullbayes_regime_calibration.tex`) have not yet been generated from a production
replication run. A one-replication smoke check was run only to verify compilation, sign convention,
and output plumbing. Those files are isolated in `simulation_outputs/smoke_check_20260714/` and
should not be cited in the manuscript.

## Notes

- Sign convention: a positive cloglog interaction gives redundant closure; complementary regimes are
  planted with negative link interactions.
- Regime is read on the link scale; magnitude is read on the response scale.
- Horseshoe hyperparameters: `scale_global = p0/((Pc-p0)*sqrt(N))` with `p0` equal to the expected
  number of active coalitions, `slab_scale = 2`, `slab_df = 4`, and `beta_scale = 2.5`.
- The full MCMC simulation is compute-heavy. The scripted 500-replication run requires 2000 Stan
  fits.
