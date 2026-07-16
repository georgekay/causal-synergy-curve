# Validation Run Report, 2026-07-15

This report records the validation runs executed in the `ready_upload_SiM_FullBayes_20260711` package on 2026-07-15.

## Code Changes Made

The simulation and bridge scripts were updated and parse-checked:

- `code/run_vb_simulation.R`
  - Adds PSIS diagnostics for ADVI fits when enabled.
  - Writes posterior means and interval endpoints for ESR, full closure, and gate summaries.
  - Writes `rep` identifiers for matching VB and MCMC rows.
- `code/run_fullbayes_simulation.R`
  - Writes `rep` identifiers, fit-failure rows, and posterior means and interval endpoints.
- `code/merge_vb_simulation_chunks.R`
  - Merges launcher chunk files with names such as `vb_simulation_raw_vb_chunk01_...csv`.
  - Handles no-PSIS agreement panels without failing on empty PSIS summaries.
- `code/merge_fullbayes_simulation_chunks.R`
  - Merges launcher chunk files with names such as `fullbayes_simulation_raw_mcmc_chunk01_...csv`.
- `code/compare_vb_mcmc_agreement.R`
  - Matches VB and MCMC rows by `scenario`, `rep`, and `truth`.
  - Reports modal-regime agreement and mean absolute differences for ESR, full closure, and gate summaries.
- `code/launch_vb_simulation_max.ps1`
  - Runs PSIS-enabled VB by default, with a `-NoPSIS` switch.
- `code/launch_vb_mcmc_agreement_panel.ps1`
  - Runs matched VB and MCMC panels and compares them.

All edited R scripts were parse-checked successfully.

## Production VB Simulation With PSIS

Command:

```powershell
& code\launch_vb_simulation_max.ps1 `
  -Reps 500 -N 8000 -Workers 60 -DrawsKeep 500 `
  -Algorithm meanfield -VBIter 20000 -ElboSamples 100 -EvalElbo 100 -GradSamples 1
```

Output directory:

```text
simulation_outputs\production_vb_20260715_175349
```

Merged output:

- 56 chunk files
- 2,000 scenario-replicate rows

Main operating characteristics:

| Scenario | Recovery | ESR coverage | Full closure coverage | Gate coverage | Failure rate |
|---|---:|---:|---:|---:|---:|
| Additive | 0.983 | 0.928 | 0.940 | 0.932 | 0.032 |
| Complementary | 0.752 | 0.975 | 0.959 | 0.973 | 0.026 |
| Mixed | 0.271 | 0.929 | 0.947 | 0.976 | 0.012 |
| Redundant | 0.946 | 0.901 | 0.944 | 0.938 | 0.032 |

Regime-probability calibration:

| Confidence bin | Mean confidence | Accuracy | n |
|---|---:|---:|---:|
| [0.4,0.5] | 0.451 | 0.093 | 151 |
| (0.5,0.6] | 0.548 | 0.458 | 236 |
| (0.6,0.7] | 0.653 | 0.642 | 176 |
| (0.7,0.8] | 0.749 | 0.902 | 143 |
| (0.8,0.9] | 0.861 | 0.798 | 198 |
| (0.9,1] | 0.986 | 0.920 | 985 |

PSIS diagnostics:

- Rows: 2,000
- Fit failures: 51
- PSIS available: 1,949
- PSIS passing `k < 0.7`: 0
- Finite Pareto-k values: 1,894
- Infinite or missing Pareto-k values: 106
- Finite Pareto-k summary: min 0.916, median 7.115, mean 6.680, 95th percentile 11.415, max 120.764

Interpretation for manuscript use:

- The VB simulation is useful as a fast operating-characteristics engine for closure, ESR, and gate summaries.
- Mean-field VB is not PSIS-certified as a close approximation to the full posterior for this Stan model. The package should not claim PSIS validation of VB.
- The mixed regime remains hard to recover under this design; this should be reported as a detectability limitation of the regime classifier under weak mixed structure.

## VB/MCMC Agreement Panel

An initial `N=8000`, 20-replicate-per-scenario agreement panel was started but stopped because the MCMC component was days-scale for simulation use. This is consistent with the decision to use VB for large simulation operating characteristics and reserve full MCMC for the applied NHANES analysis.

A compact full-MCMC agreement panel was then run:

```powershell
& code\launch_vb_mcmc_agreement_panel.ps1 `
  -RepsPerScenario 8 -N 2000 `
  -VBWorkers 8 -MCMCWorkers 8 -ChainsPerFit 4 `
  -VBAlgorithm meanfield -VBIter 12000 -DrawsKeep 400 `
  -IterWarmup 300 -IterSampling 300 -AdaptDelta 0.95 -MaxTreeDepth 12 `
  -NoPSIS
```

Output directory:

```text
simulation_outputs\vb_mcmc_agreement_20260715_203839
```

Matched rows:

- 32 total rows
- 8 rows per scenario

Agreement summary:

| Scenario | n | Modal agreement | Mean abs ESR diff | Mean abs full-closure diff | Mean abs gate diff | VB recovery | MCMC recovery |
|---|---:|---:|---:|---:|---:|---:|---:|
| Additive | 8 | 1.000 | 0.598 | 0.009 | 0.156 | 1.000 | 1.000 |
| Complementary | 8 | 0.500 | 1.550 | 0.045 | 0.175 | 0.500 | 1.000 |
| Mixed | 8 | 0.250 | 0.274 | 0.004 | 0.129 | 0.000 | 0.625 |
| Redundant | 8 | 0.750 | 0.080 | 0.050 | 0.055 | 0.625 | 0.875 |
| Overall | 32 | 0.625 | 0.626 | 0.027 | 0.129 | 0.531 | 0.875 |

Interpretation for manuscript use:

- MCMC recovers the simulated regimes more often than mean-field VB in this compact panel.
- Full-closure agreement is good in absolute terms; mean absolute differences are small.
- ESR and modal-regime agreement are not close enough to treat mean-field VB as a substitute for MCMC for regime probabilities.
- The safest reporting structure is:
  - Use VB for large-scale simulation operating characteristics and speed.
  - Use full MCMC for the applied NHANES analysis and any final posterior regime probabilities.
  - State that mean-field VB is a computational approximation for simulation screening, not the inferential gold standard for the final applied posterior.

## Applied NHANES Full MCMC

Applied full-Bayes outputs already exist in:

```text
applied_outputs
```

Key scalar summaries:

| Cohort | n | 10-year full closure, pp | 95% CrI | Response ESR | 95% CrI | Link ESR | 95% CrI | All-favorable support |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Pooled | 26,353 | 3.989 | [3.651, 4.313] | -0.053 | [-0.208, 0.065] | 0.176 | [0.059, 0.267] | 0.255 |
| Black | 8,968 | 4.867 | [4.124, 5.593] | 0.018 | [-0.334, 0.345] | 0.229 | [-0.080, 0.493] | 0.091 |
| White | 17,385 | 3.948 | [3.602, 4.323] | -0.069 | [-0.235, 0.072] | 0.182 | [0.067, 0.286] | 0.282 |

Sampler diagnostics:

| Cohort | Chains | Sampling iterations | Max Rhat | Min bulk ESS | Min tail ESS | Divergences | Max treedepth |
|---|---:|---:|---:|---:|---:|---:|---:|
| Pooled | 4 | 1,000 | 1.007 | 537 | 303 | 0 | 12 |
| Black | 4 | 1,000 | 1.007 | 531 | 894 | 0 | 12 |
| White | 4 | 1,000 | 1.007 | 562 | 990 | 0 | 11 |

Selected gate results:

- Pooled education gate: `G_link = 0.371`, CrI `[0.218, 0.541]`, `P_gate = 1.000`.
- White education gate: `G_link = 0.310`, CrI `[0.181, 0.497]`, `P_gate = 1.000`.
- Black education gate: `G_link = 0.158`, CrI `[-0.043, 0.421]`, `P_gate = 0.928`.
- Pooled partnership gate: `G_link = 0.083`, CrI `[0.003, 0.172]`, `P_gate = 0.980`.
- Income is the largest single-domain closure contributor, but not the primary gate: pooled `G_link = 0.059`, CrI `[-0.081, 0.167]`.

Interpretation for manuscript use:

- Applied MCMC diagnostics are acceptable.
- The final applied results should rely on the full MCMC outputs, not VB.
- The response-scale ESR is near additive and crosses zero; the link-scale ESR is mildly positive in pooled and White cohorts and uncertain in Black.
- The most distinctive applied role separation remains: income is the dominant closure contributor, while education is the clearest gate in pooled and White analyses.

## Current Status

Completed:

- PSIS-enabled production VB simulation.
- Compact VB/MCMC agreement panel.
- Primary full-MCMC simulation panel with five known structures, including a planted gate.
- Applied NHANES full-MCMC output inspection.
- Script parse checks after edits.

Not recommended:

- Claiming PSIS validates mean-field VB. The diagnostic rejects that.
- Claiming VB and MCMC give interchangeable regime probabilities. The compact agreement panel does not support that.
- Running full MCMC for every simulation replicate. The attempted `N=8000`, 20-rep-per-scenario panel demonstrated that this is not practical for the simulation study.

Recommended manuscript framing:

1. Full MCMC is the primary inferential engine for the applied NHANES analysis.
2. VB is a fast simulation/screening engine used for broad operating characteristics, with PSIS diagnostics reported.
3. The method has strong closure/gate interval behavior in the VB simulation, but mixed-regime classification is difficult under weak mixed truth.
4. The final applied claims should be tied to the full-MCMC NHANES tables and sampler diagnostics.

## Primary Full-MCMC Simulation Panel

After the VB/MCMC bridge, a manuscript-primary full-MCMC simulation panel was run with five known
data-generating structures:

- additive
- complementary
- redundant
- mixed
- domain-anchored gate

Final strict sampler command:

```powershell
& code\launch_fullbayes_simulation_max.ps1 `
  -Reps 25 -N 2000 -Workers 20 -ChainsPerFit 4 `
  -IterWarmup 500 -IterSampling 600 -DrawsKeep 500 `
  -AdaptDelta 0.99 -MaxTreeDepth 13
```

Output directory:

```text
simulation_outputs\production_fullbayes_20260715_213021
```

The final strict run completed 125/125 scenario-replication fits and was promoted to the standard
manuscript input paths:

```text
simulation_outputs\fullbayes_simulation_raw.csv
simulation_outputs\fullbayes_recovery_coverage.csv
simulation_outputs\fullbayes_recovery_coverage.tex
simulation_outputs\fullbayes_regime_calibration.csv
simulation_outputs\fullbayes_regime_calibration.tex
```

Main operating characteristics:

| Scenario | Regime recovery | ESR coverage | Full closure coverage | Gate-index coverage | Gate-domain recovery | Mean divergences |
|---|---:|---:|---:|---:|---:|---:|
| Additive | 1.00 | 1.00 | 0.96 | 0.92 | 0.08 | 0.00 |
| Complementary | 0.92 | 1.00 | 0.96 | 0.92 | 0.32 | 0.00 |
| Gate | 0.92 | 1.00 | 0.92 | 0.96 | 1.00 | 0.00 |
| Mixed | 0.76 | 0.92 | 0.88 | 0.88 | 0.56 | 0.00 |
| Redundant | 0.96 | 0.88 | 0.92 | 0.96 | 0.00 | 0.00 |

Diagnostics:

- Fitted models: 125
- Fit failures: 0
- Total divergent transitions: 0
- Fits with divergent transitions: 0
- Maximum model-level Rhat: 1.023
- Minimum bulk ESS: 256.9
- Minimum tail ESS: 227.5

Generated simulation figures:

```text
figures\fig_fullbayes_sim_recovery_coverage.pdf/.png
figures\fig_fullbayes_sim_regime_confusion.pdf/.png
figures\fig_fullbayes_sim_calibration.pdf/.png
```

The simulation manuscript fragment was updated:

```text
manuscript\simulation_fullbayes_revised.tex
```

This strict full-MCMC panel is the recommended primary simulation evidence for the SiM manuscript.
