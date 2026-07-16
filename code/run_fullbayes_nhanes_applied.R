## =====================================================================================
## Full-Bayes applied illustration: NHANES 8-SDOH -> 10-year all-cause mortality.
## Weibull PH with the SAME regularized-horseshoe coalition surface as the simulation
## (baseline + main effects + burden-depth + anchored gates + horseshoe coalition block),
## fit by full Stan MCMC. Posterior g-computation at horizon t = 10 with weighted (pooled MEC)
## standardization gives response-scale closure (absolute risk) and link-scale structure.
##
## Run:  Rscript run_fullbayes_nhanes_applied.R
## Needs: cmdstanr (+ CmdStan), csc_design_readouts.R, and the NHANES analytic CSV (path below).
## =====================================================================================
suppressPackageStartupMessages({ library(cmdstanr) })
source("csc_design_readouts.R")

DATA   <- Sys.getenv("CSC_NHANES", "../data/nhanes_8sdoh_analytic.csv")
COHORT <- Sys.getenv("CSC_COHORT", "pooled")     # pooled | black | white
HORIZON<- as.numeric(Sys.getenv("CSC_HORIZON", "10"))
ORDER <- as.integer(Sys.getenv("CSC_ORDER", "2"))
DRAWS_KEEP <- as.integer(Sys.getenv("CSC_DRAWS_KEEP", "1000"))
CHAINS <- as.integer(Sys.getenv("CSC_CHAINS", "4"))
ITER_WARM <- as.integer(Sys.getenv("CSC_ITER_WARMUP", "750"))
ITER_SAMP <- as.integer(Sys.getenv("CSC_ITER_SAMPLING", "750"))
ADAPT_DELTA <- as.numeric(Sys.getenv("CSC_ADAPT_DELTA", "0.99"))
MAX_TREEDEPTH <- as.integer(Sys.getenv("CSC_MAX_TREEDEPTH", "12"))
GCOMP_N <- as.integer(Sys.getenv("CSC_GCOMP_N", "0"))  # 0 = exact/full standardization sample
CMDSTAN_CSV_DIR <- Sys.getenv("CSC_CMDSTAN_CSV_DIR", "")
CMDSTAN_OUTDIR <- Sys.getenv("CSC_CMDSTAN_OUTDIR", "")
STAN   <- "../stan/csc_coalition_horseshoe_weibull.stan"
OUT    <- "../applied_outputs"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DOM <- c("D_employment","D_income","D_food_security","D_education",
         "D_access_to_care","D_health_insurance","D_housing_instability","D_partnership")

d <- read.csv(DATA)
if (COHORT == "black") d <- d[grepl("Black|black|NHB|^2$", as.character(d$analysis_race)), ]
if (COHORT == "white") d <- d[grepl("White|white|NHW|^1$", as.character(d$analysis_race)), ]
d <- d[complete.cases(d[, c(DOM, "analysis_age", "sex_female", "followup_years", "MORTSTAT")]), ]
d <- d[d$followup_years > 0, ]
n <- nrow(d)
D <- as.matrix((d[, DOM] > 0) * 1)                        # 1 = unfavorable
W <- cbind(age = as.numeric(scale(d$analysis_age)), sex = d$sex_female)
time <- d$followup_years; event <- as.integer(d$MORTSTAT > 0)
wt <- if ("pooled_mec_weight" %in% names(d)) d$pooled_mec_weight else rep(1, n)
wt <- wt / mean(wt)
u_ref <- sapply(seq_len(ncol(D)), function(g) { x <- D[D[, g] > 0, g]; if (length(x)) mean(x) else 0 })

meta <- csc_build_design(D, W, ORDER); Pc <- ncol(meta$Xc); p <- ncol(D); p0 <- 3
sdat <- list(N = n, Pb = ncol(meta$Xb), Pc = Pc, Xb = meta$Xb, Xc = meta$Xc,
             time = time, event = event,
             scale_global = p0 / ((Pc - p0) * sqrt(n)), slab_scale = 2, slab_df = 4, beta_scale = 2.5)
if (nzchar(CMDSTAN_CSV_DIR)) {
  csv_files <- list.files(CMDSTAN_CSV_DIR, pattern = "\\.csv$", full.names = TRUE)
  if (!length(csv_files)) stop("CSC_CMDSTAN_CSV_DIR contains no CmdStan CSV files: ", CMDSTAN_CSV_DIR)
  fit <- as_cmdstan_fit(files = csv_files)
} else {
  mod <- cmdstan_model(STAN)
  if (nzchar(CMDSTAN_OUTDIR)) dir.create(CMDSTAN_OUTDIR, showWarnings = FALSE, recursive = TRUE)
  sample_args <- list(data = sdat, chains = CHAINS, parallel_chains = CHAINS,
                      iter_warmup = ITER_WARM, iter_sampling = ITER_SAMP,
                      adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH,
                      init = 0, refresh = 50)
  if (nzchar(CMDSTAN_OUTDIR)) sample_args$output_dir <- CMDSTAN_OUTDIR
  fit <- do.call(mod$sample, sample_args)
}
sm <- fit$summary()
sd <- fit$sampler_diagnostics(format = "df")
dr <- fit$draws(format = "df")
chains_recorded <- if (nzchar(CMDSTAN_CSV_DIR)) length(csv_files) else CHAINS
iter_sampling_recorded <- if (nzchar(CMDSTAN_CSV_DIR)) floor(nrow(dr) / chains_recorded) else ITER_SAMP
diag_tbl <- data.frame(
  cohort = COHORT,
  n = n,
  chains = chains_recorded,
  iter_warmup = if (nzchar(CMDSTAN_CSV_DIR)) NA_integer_ else ITER_WARM,
  iter_sampling = iter_sampling_recorded,
  max_rhat = suppressWarnings(max(sm$rhat, na.rm = TRUE)),
  min_ess_bulk = suppressWarnings(min(sm$ess_bulk, na.rm = TRUE)),
  min_ess_tail = suppressWarnings(min(sm$ess_tail, na.rm = TRUE)),
  divergences = if ("divergent__" %in% names(sd)) sum(sd$divergent__, na.rm = TRUE) else NA_integer_,
  max_treedepth = if ("treedepth__" %in% names(sd)) max(sd$treedepth__, na.rm = TRUE) else NA_integer_
)
set.seed(as.integer(Sys.getenv("CSC_READOUT_SEED", "20260711")))
idx <- sort(sample(nrow(dr), min(DRAWS_KEEP, nrow(dr))))
gb <- function(pat) as.matrix(dr[idx, grep(pat, colnames(dr)), drop = FALSE])
post <- list(alpha0 = dr[[grep("^alpha0$", colnames(dr))]][idx],
             beta_b = gb("^beta_b\\["), theta = gb("^theta\\["))
shp  <- dr[[grep("^shape$", colnames(dr))]][idx]

## response = absolute risk at HORIZON; link = linear predictor (structure)
resp_fun <- function(E) 1 - exp(-sweep(exp(pmin(E, 30)), 2, HORIZON^shp, "*"))
if (GCOMP_N > 0 && GCOMP_N < n) {
  set.seed(20260711)
  gidx <- sort(sample(seq_len(n), GCOMP_N))
} else {
  gidx <- seq_len(n)
}
Dg <- D[gidx, , drop = FALSE]
Wg <- W[gidx, , drop = FALSE]
wtg <- wt[gidx]
wtg <- wtg / mean(wtg)
clo <- csc_closures(Dg, Wg, meta, post, alpha = 1, resp_fun = resp_fun, w = wtg)
M <- clo$M

## per-draw functionals
ESR_link <- ESR_resp <- full_resp <- numeric(M)
phi_resp <- matrix(0, p, M); pair_link <- NULL; regime <- character(M)
for (m in seq_len(M)) {
  fl <- csc_functionals(clo$link[, m], p); fr <- csc_functionals(clo$resp[, m], p)
  ESR_link[m] <- fl$ESR; ESR_resp[m] <- fr$ESR; full_resp[m] <- fr$full
  phi_resp[, m] <- fr$phi; regime[m] <- csc_regime(fl$pair_delta, 0.10)
  pair_link <- if (is.null(pair_link)) matrix(0, length(fl$pair_delta), M) else pair_link
  pair_link[, m] <- fl$pair_delta
}
gt <- csc_gates(Dg, Wg, meta, post, u_ref, w = wtg)   # link-scale gate B/E/G

qs <- function(x) c(mean = mean(x), lo = quantile(x, .025), hi = quantile(x, .975))
arch <- data.frame(domain = DOM,
                   closure_pp = rowMeans(phi_resp) * 100,
                   lo = apply(phi_resp, 1, quantile, .025) * 100,
                   hi = apply(phi_resp, 1, quantile, .975) * 100,
                   G_link = rowMeans(gt$G), G_lo = apply(gt$G, 1, quantile, .025),
                   G_hi = apply(gt$G, 1, quantile, .975), P_gate = rowMeans(gt$G > 0))
regime_probs <- prop.table(table(factor(regime, levels = c("additive","complementary","redundant","mixed"))))
fav_all <- weighted.mean(rowSums(D) == 0, wt)
summary_tbl <- data.frame(
  cohort = COHORT, n = n, horizon = HORIZON,
  n_gcomp = length(gidx),
  full_closure_pp = mean(full_resp) * 100,
  full_lo = quantile(full_resp, .025) * 100, full_hi = quantile(full_resp, .975) * 100,
  ESR_resp = mean(ESR_resp), ESR_resp_lo = quantile(ESR_resp, .025), ESR_resp_hi = quantile(ESR_resp, .975),
  ESR_link = mean(ESR_link), ESR_link_lo = quantile(ESR_link, .025), ESR_link_hi = quantile(ESR_link, .975),
  P_additive = regime_probs["additive"], P_complementary = regime_probs["complementary"],
  P_redundant = regime_probs["redundant"], P_mixed = regime_probs["mixed"], all_favorable = fav_all)

pr <- paste0("_", COHORT)
write.csv(arch, file.path(OUT, paste0("applied_architecture_gates", pr, ".csv")), row.names = FALSE)
write.csv(summary_tbl, file.path(OUT, paste0("applied_summary", pr, ".csv")), row.names = FALSE)
write.csv(diag_tbl, file.path(OUT, paste0("applied_sampler_diagnostics", pr, ".csv")), row.names = FALSE)
write.csv(sm, file.path(OUT, paste0("applied_parameter_summary", pr, ".csv")), row.names = FALSE)
saveRDS(list(pair_link_mean = rowMeans(pair_link), post_shape_mean = mean(shp)),
        file.path(OUT, paste0("applied_pairwise", pr, ".rds")))
cat("\n== cohort:", COHORT, " n =", n, "==\n")
cat(sprintf("Full favorable shift closes %.1f pp (95%% CrI %.1f, %.1f) of %d-year risk.\n",
            mean(full_resp)*100, quantile(full_resp,.025)*100, quantile(full_resp,.975)*100, HORIZON))
cat(sprintf("Regime: additive %.0f%%, complementary %.0f%%, redundant %.0f%%, mixed %.0f%%.\n",
            100*regime_probs["additive"],100*regime_probs["complementary"],
            100*regime_probs["redundant"],100*regime_probs["mixed"]))
print(arch, digits = 3)
print(diag_tbl, digits = 3)
cat("Wrote outputs to", OUT, "\n")
