# Full-Bayes package audit and run log, 2026-07-14

Package audited:

`D:\AMC_Papers\Causal_Synergy_Curve_Methods_Paper\ready_upload_SiM_FullBayes_20260711`

## Package status

This package is a Full-Bayes revision module for the SiM manuscript. It now contains:

- Stan models for the binary simulation and Weibull PH NHANES applied analysis;
- R scripts for simulation, applied MCMC, shared CSC readouts, and applied table/figure generation;
- the NHANES 8-SDOH analytic CSV copied into `data/`;
- completed applied MCMC outputs for pooled, Black, and White cohorts;
- applied tables and figures generated from the MCMC readouts;
- a quarantined simulation smoke check.

It is not yet a complete final submission package because the production Full-Bayes simulation tables
have not been generated and the manuscript directory contains drop-in fragments rather than the full
compiled manuscript source.

## Code fixes made during the audit

1. Stabilized the binary cloglog Stan likelihood against overflow by evaluating the Bernoulli
   likelihood on a clipped cumulative-hazard scale.
2. Vectorized and stabilized the Weibull PH Stan likelihood.
3. Added `init = 0` to Stan calls to avoid unstable random initial values.
4. Made the applied harness resumable from preserved CmdStan CSVs via `CSC_CMDSTAN_CSV_DIR`.
5. Added environment controls for chains, iterations, tree depth, g-computation sample size, output
   directory, draw-retention seed, and horizon.
6. Added MCMC diagnostic output files for applied analyses.
7. Added simulation truth-gate computation and gate-coverage output.
8. Updated simulation sign strength so planted complementary/redundant regimes clear the ROPE under
   the default domain prevalence.
9. Fixed the shared design builder so `order = 1` correctly creates no coalition block.
10. Added `code/make_fullbayes_applied_tables_figures.R` to generate paper tables and figures from
    the applied MCMC output files.
11. Added `code/launch_fullbayes_simulation_max.ps1` and
    `code/merge_fullbayes_simulation_chunks.R` so the production simulation can be sharded across
    the workstation. On this machine the launcher defaults to 20 concurrent workers x 4 chains per
    fit, using the 80 logical processors.

## Data audit

NHANES analytic dataset:

`data/nhanes_8sdoh_analytic.csv`

Required analysis columns were present. Final cohort sizes used by the applied MCMC:

| Cohort | n | deaths |
|---|---:|---:|
| Pooled | 26353 | 3022 |
| Black | 8968 | 1030 |
| White | 17385 | 1992 |

## Applied MCMC settings

For each cohort:

- model: Weibull PH regularized-horseshoe coalition surface;
- domains: 8 SDOH;
- coalition order: pairwise (`CSC_ORDER=2`);
- horizon: 10 years;
- chains: 4;
- warmup iterations per chain: 1000;
- post-warmup iterations per chain: 1000;
- adapt delta: 0.99;
- maximum tree depth: 12;
- posterior readout: 250 retained posterior draws, deterministic seed `20260711`;
- standardization sample: full cohort for each cohort-specific readout.

CmdStan CSVs are preserved in:

- `applied_outputs/cmdstan_csv_pooled_20260714/`
- `applied_outputs/cmdstan_csv_black_20260714/`
- `applied_outputs/cmdstan_csv_white_20260714/`

## MCMC diagnostics

| Cohort | max Rhat | min bulk ESS | min tail ESS | divergences | max treedepth |
|---|---:|---:|---:|---:|---:|
| Pooled | 1.007 | 537 | 303 | 0 | 12 |
| Black | 1.007 | 531 | 894 | 0 | 12 |
| White | 1.007 | 562 | 990 | 0 | 11 |

Stan reported 11/4000 transitions at max tree depth for pooled and 1/4000 for Black. White did not
hit the maximum tree depth. No divergent transitions were detected.

## Applied MCMC results

Closure is reported on the response scale in percentage points. ESR is reported on both response and
link scales; regime classification is based on the link-scale posterior second-difference family.

| Cohort | Full closure pp (95% CrI) | Response ESR (95% CrI) | Link ESR (95% CrI) | Regime posterior |
|---|---:|---:|---:|---|
| Pooled | 3.99 (3.65, 4.31) | -0.053 (-0.208, 0.065) | 0.176 (0.059, 0.267) | additive 1.000 |
| Black | 4.87 (4.12, 5.59) | 0.018 (-0.334, 0.345) | 0.229 (-0.080, 0.493) | additive 0.844, complementary 0.132, redundant 0.012, mixed 0.012 |
| White | 3.95 (3.60, 4.32) | -0.069 (-0.235, 0.072) | 0.182 (0.067, 0.286) | additive 1.000 |

Observed all-favorable support:

- pooled: 25.5%;
- Black: 9.1%;
- White: 28.2%.

## Domain-role results

Pooled domain closures and gate indices:

| Domain | Closure pp (95% CrI) | Gate index G (95% CrI) |
|---|---:|---:|
| Employment | 0.74 (0.57, 0.89) | -0.016 (-0.100, 0.076) |
| Income | 1.12 (0.88, 1.41) | 0.059 (-0.081, 0.167) |
| Food security | 0.23 (0.08, 0.40) | 0.088 (-0.028, 0.206) |
| Education | 0.48 (0.34, 0.63) | 0.371 (0.218, 0.541) |
| Access to care | 0.07 (-0.03, 0.18) | 0.038 (-0.057, 0.139) |
| Health insurance | 0.60 (0.40, 0.80) | 0.046 (-0.051, 0.135) |
| Housing instability | 0.13 (-0.02, 0.28) | 0.035 (-0.045, 0.122) |
| Partnership | 0.62 (0.42, 0.83) | 0.083 (0.003, 0.172) |

Interpretation for the applied illustration:

- income is the largest single-domain closure lever in the pooled fit;
- education is the clearest pooled gate, with the other-domain closure larger when education is
  favorable;
- the global surface is predominantly additive by the posterior regime classifier, with a positive
  link-scale ESR but no non-additive regime label under the pairwise ROPE rule.

## Generated applied artifacts

Tables:

- `tables/fullbayes_applied_summary.csv`
- `tables/fullbayes_applied_summary.tex`
- `tables/fullbayes_domain_roles.csv`
- `tables/fullbayes_domain_roles.tex`

Figures:

- `figures/fig_fullbayes_closure_by_cohort.pdf`
- `figures/fig_fullbayes_closure_by_cohort.png`
- `figures/fig_fullbayes_domain_architecture.pdf`
- `figures/fig_fullbayes_domain_architecture.png`
- `figures/fig_fullbayes_gate_indices.pdf`
- `figures/fig_fullbayes_gate_indices.png`

Visual QA: all three PNG figures were opened after generation. The plots render correctly, labels are
legible, intervals are visible, and no plot is blank or malformed.

Applied raw and summary outputs:

- `applied_outputs/applied_summary_*.csv`
- `applied_outputs/applied_architecture_gates_*.csv`
- `applied_outputs/applied_pairwise_*.csv`
- `applied_outputs/applied_sampler_diagnostics_*.csv`
- `applied_outputs/applied_parameter_summary_*.csv`

## Remaining before SiM submission

1. Run the production Full-Bayes simulation. The required top-level simulation files are not present:
   `simulation_outputs/fullbayes_recovery_coverage.tex` and
   `simulation_outputs/fullbayes_regime_calibration.tex`.
2. Integrate the drop-in Full-Bayes methods/simulation text and applied tables/figures into the main
   manuscript source.
3. Recompile the manuscript and supplement, then audit references, figure/table numbering, and
   cross-references.
4. Decide whether the production simulation replication count is 500 as scripted or a smaller
   explicitly labeled run. The current 500-replication plan entails 2000 Stan fits.

The max-resource launcher for item 1 is:

```powershell
cd D:\AMC_Papers\Causal_Synergy_Curve_Methods_Paper\ready_upload_SiM_FullBayes_20260711\code
.\launch_fullbayes_simulation_max.ps1 -Reps 500 -N 8000 -ChainsPerFit 4
```

After all workers finish, set `CSC_OUTDIR` to the production output directory printed by the
launcher and run:

```powershell
Rscript merge_fullbayes_simulation_chunks.R
```
