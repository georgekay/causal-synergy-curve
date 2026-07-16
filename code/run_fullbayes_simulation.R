## =====================================================================================
## Full-Bayes operating-characteristic simulation for the Causal Synergy Curve.
## FULL posterior sampling (Stan MCMC) of the regularized-horseshoe coalition surface --
## no Laplace/normal approximation -- matched to the applied estimator.
## Reports: (1) regime recovery, (2) posterior credible-interval COVERAGE of closure/ESR/gate,
##          (3) posterior REGIME-PROBABILITY CALIBRATION (confidence vs accuracy).
## Truth is g-computed from the true cloglog surface with the same definitions the estimator uses.
##
## Run:  Rscript run_fullbayes_simulation.R          (edit REPS/CHAINS/etc. below)
## Needs: cmdstanr (+ CmdStan), and csc_design_readouts.R in the same folder.
## =====================================================================================
suppressPackageStartupMessages({ library(cmdstanr) })
source("csc_design_readouts.R")
set.seed(20260711)
options(warn = 1)

## ---- configuration (reduce REPS for a fast first pass; 500 for the paper) -------------
ALL_SCEN  <- c("additive", "complementary", "redundant", "mixed", "gate")
SCEN_ENV  <- Sys.getenv("CSC_SCENARIOS", paste(ALL_SCEN, collapse = ","))
SCEN      <- trimws(strsplit(SCEN_ENV, ",", fixed = TRUE)[[1]])
SCEN      <- SCEN[nzchar(SCEN)]
if (!all(SCEN %in% ALL_SCEN)) stop("CSC_SCENARIOS must be drawn from: ", paste(ALL_SCEN, collapse = ", "))
P         <- 6
N         <- as.integer(Sys.getenv("CSC_N", "8000"))
REPS      <- as.integer(Sys.getenv("CSC_REPS", "500"))
REP_START <- as.integer(Sys.getenv("CSC_REP_START", "1"))
REP_END   <- as.integer(Sys.getenv("CSC_REP_END", as.character(REPS)))
REP_SEQ   <- seq.int(REP_START, REP_END)
ORDER     <- as.integer(Sys.getenv("CSC_ORDER", "2"))
DRAWS_KEEP<- as.integer(Sys.getenv("CSC_DRAWS_KEEP", "500"))
CHAINS    <- as.integer(Sys.getenv("CSC_CHAINS", "4"))
ITER_WARM <- as.integer(Sys.getenv("CSC_ITER_WARMUP", "500"))
ITER_SAMP <- as.integer(Sys.getenv("CSC_ITER_SAMPLING", "500"))
ADAPT_DELTA <- as.numeric(Sys.getenv("CSC_ADAPT_DELTA", "0.95"))
MAX_TREEDEPTH <- as.integer(Sys.getenv("CSC_MAX_TREEDEPTH", "12"))
ROPE      <- as.numeric(Sys.getenv("CSC_ROPE", "0.10"))  # practical-equivalence region for second differences
GAMMA     <- as.numeric(Sys.getenv("CSC_GAMMA", "1.50")) # link interaction strength; 1.5 clears ROPE at Pr(D_jD_k) ~= .09
MAINS     <- seq(0.55, 0.20, length.out = P)
STAN_FILE <- "../stan/csc_coalition_horseshoe_binary.stan"
OUTDIR    <- Sys.getenv("CSC_OUTDIR", "../simulation_outputs")
JOB_ID    <- Sys.getenv("CSC_JOB_ID", "")
TAG       <- if (nzchar(JOB_ID)) paste0("_", JOB_ID) else ""
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
RAW_PATH  <- file.path(OUTDIR, paste0("fullbayes_simulation_raw", TAG, ".csv"))
if (file.exists(RAW_PATH)) file.remove(RAW_PATH)

mod <- cmdstan_model(STAN_FILE)

## ---- data-generating surfaces (SIGN CONVENTION: +link interaction => REDUNDANT closure) ----
## To plant a COMPLEMENTARY closure regime we use NEGATIVE link interactions, and vice versa.
eta_true_fun <- function(scen) function(D, age) {
  e <- -2.2 + 0.5 * age + as.numeric(D %*% MAINS)
  if (scen == "complementary") e <- e - GAMMA * (D[,1]*D[,2] + D[,3]*D[,4])
  if (scen == "redundant")     e <- e + GAMMA * (D[,1]*D[,2] + D[,3]*D[,4])
  if (scen == "mixed")         e <- e - GAMMA * D[,1]*D[,2] + GAMMA * D[,3]*D[,4]
  if (scen == "gate") {
    ## Domain 1 is the gate: deficits 2:P matter much more when domain 1 is favorable.
    ## Algebraically this is a domain-anchored set of negative pairwise link interactions:
    ##   D_j * (1-D_1) = D_j - D_1 D_j.
    e <- -2.2 + 0.5 * age + MAINS[1] * D[,1] +
      0.15 * rowSums(D[, 2:P, drop = FALSE]) +
      0.70 * (1 - D[,1]) * rowSums(D[, 2:P, drop = FALSE])
  }
  e
}

gen_data <- function(scen, seed) {
  set.seed(seed)
  D   <- matrix(rbinom(N * P, 1, 0.30), N, P)
  age <- rnorm(N)
  ef  <- eta_true_fun(scen)
  y   <- rbinom(N, 1, inv_cloglog(ef(D, age)))
  list(D = D, W = matrix(age, ncol = 1), y = y, ef = ef, age = age)
}

## true closure vector (2^P) on link + response, then true functionals
truth_functionals <- function(dat) {
  D <- dat$D; age <- dat$age; ef <- dat$ef; p <- P
  e0 <- ef(D, age); r0 <- inv_cloglog(e0)
  cl_link <- numeric(2^p); cl_resp <- numeric(2^p)
  for (mask in 1:(2^p - 1)) {
    S <- which(bitwAnd(mask, 2^(0:(p-1))) > 0)
    Dp <- D; Dp[, S] <- 0
    cl_link[mask + 1] <- mean(e0 - ef(Dp, age))
    cl_resp[mask + 1] <- mean(r0 - inv_cloglog(ef(Dp, age)))
  }
  fl <- csc_functionals(cl_link, p); fr <- csc_functionals(cl_resp, p)
  gt <- truth_gates(dat)
  gate_domain <- which.max(abs(gt$G))
  regime <- if (any(abs(fl$pair_delta) > ROPE)) csc_regime(fl$pair_delta, ROPE) else "additive"
  list(regime = regime, ESR_link = fl$ESR,
       full_resp = fr$full, phi_resp = fr$phi,
       gate_domain = gate_domain, gate_G = gt$G[gate_domain])
}

truth_gates <- function(dat) {
  D <- dat$D; W <- dat$W; age <- dat$age; ef <- dat$ef; p <- P
  Bg <- Eg <- numeric(p)
  for (g in seq_len(p)) {
    oth <- setdiff(seq_len(p), g)
    Bu <- D; Bu[, g] <- 1; Bu2 <- Bu; Bu2[, oth] <- 0
    Ef <- D; Ef[, g] <- 0; Ef2 <- Ef; Ef2[, oth] <- 0
    Bg[g] <- mean(ef(Bu, age) - ef(Bu2, age))
    Eg[g] <- mean(ef(Ef, age) - ef(Ef2, age))
  }
  G <- (Eg - Bg) / (abs(Eg) + abs(Bg) + 1e-9)
  list(B = Bg, E = Eg, G = G)
}

fit_one <- function(dat) {
  meta <- csc_build_design(dat$D, dat$W, ORDER)
  Pc <- ncol(meta$Xc); p0 <- 2
  sdat <- list(N = N, Pb = ncol(meta$Xb), Pc = Pc, Xb = meta$Xb, Xc = meta$Xc, y = as.integer(dat$y),
               scale_global = p0 / ((Pc - p0) * sqrt(N)), slab_scale = 2, slab_df = 4, beta_scale = 2.5)
  fit <- mod$sample(data = sdat, chains = CHAINS, parallel_chains = CHAINS,
                    iter_warmup = ITER_WARM, iter_sampling = ITER_SAMP, refresh = 0,
                    adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH,
                    init = 0, show_messages = FALSE)
  sm <- fit$summary()
  sd <- fit$sampler_diagnostics(format = "df")
  diag <- data.frame(
    max_rhat = suppressWarnings(max(sm$rhat, na.rm = TRUE)),
    min_ess_bulk = suppressWarnings(min(sm$ess_bulk, na.rm = TRUE)),
    min_ess_tail = suppressWarnings(min(sm$ess_tail, na.rm = TRUE)),
    divergences = if ("divergent__" %in% names(sd)) sum(sd$divergent__, na.rm = TRUE) else NA_integer_,
    max_treedepth = if ("treedepth__" %in% names(sd)) max(sd$treedepth__, na.rm = TRUE) else NA_integer_
  )
  dr <- fit$draws(format = "df")
  idx <- sort(sample(nrow(dr), min(DRAWS_KEEP, nrow(dr))))
  gb <- function(pat) as.matrix(dr[idx, grep(pat, colnames(dr)), drop = FALSE])
  post <- list(alpha0 = dr[[grep("^alpha0$", colnames(dr))]][idx],
               beta_b = gb("^beta_b\\["), theta = gb("^theta\\["))
  list(meta = meta, post = post, diag = diag)
}

## posterior readouts -> per-draw ESR_link, regime label, full response closure, gate G(max)
posterior_readouts <- function(dat, meta, post) {
  clo <- csc_closures(dat$D, dat$W, meta, post, alpha = 1, resp_fun = inv_cloglog)
  M <- clo$M
  ESR_link <- numeric(M); regime <- character(M); full_resp <- numeric(M)
  for (m in seq_len(M)) {
    fl <- csc_functionals(clo$link[, m], P); fr <- csc_functionals(clo$resp[, m], P)
    ESR_link[m] <- fl$ESR; regime[m] <- csc_regime(fl$pair_delta, ROPE); full_resp[m] <- fr$full
  }
  gt <- csc_gates(dat$D, dat$W, meta, post, rep(1, P))
  list(ESR_link = ESR_link, regime = regime, full_resp = full_resp, gate_G = gt$G)
}

## ---- run ------------------------------------------------------------------------------
rows <- list(); calib <- list(); ri <- 0
TOTAL <- length(SCEN) * length(REP_SEQ)
for (scen in SCEN) for (rep in REP_SEQ) {
  ri <- ri + 1
  dat <- gen_data(scen, seed = 1000 * match(scen, ALL_SCEN) + rep)
  tru <- truth_functionals(dat)
  ff  <- tryCatch(fit_one(dat), error = function(e) NULL)
  if (is.null(ff)) {
    rows[[ri]] <- data.frame(
      scenario = scen, rep = rep, truth = tru$regime, modal = "failed", correct = FALSE,
      p_true = NA_real_, p_modal = NA_real_,
      esr_mean = NA_real_, esr_lo = NA_real_, esr_hi = NA_real_, esr_cov = NA,
      full_mean = NA_real_, full_lo = NA_real_, full_hi = NA_real_, full_cov = NA,
      gate_domain = tru$gate_domain,
      gate_G_truth = tru$gate_G,
      gate_G_mean = NA_real_, gate_G_lo = NA_real_, gate_G_hi = NA_real_, gate_cov = NA,
      gate_modal_domain = NA_integer_, gate_domain_correct = NA,
      max_rhat = NA_real_, min_ess_bulk = NA_real_, min_ess_tail = NA_real_,
      divergences = NA_integer_, max_treedepth = NA_integer_,
      fit_failed = TRUE
    )
    cat(sprintf("[%d/%d] rep=%d %-13s truth=%-13s modal=failed fit error\n",
                ri, TOTAL, rep, scen, tru$regime))
    write.table(rows[[ri]], RAW_PATH, sep = ",", row.names = FALSE,
                col.names = !file.exists(RAW_PATH), append = file.exists(RAW_PATH))
    flush.console()
    next
  }
  pr  <- posterior_readouts(dat, ff$meta, ff$post)
  probs <- prop.table(table(factor(pr$regime, levels = c("additive","complementary","redundant","mixed"))))
  modal <- names(which.max(probs))
  gdraw <- pr$gate_G[tru$gate_domain, ]
  esr_q <- as.numeric(quantile(pr$ESR_link, c(.025, .975), na.rm = TRUE, names = FALSE))
  full_q <- as.numeric(quantile(pr$full_resp, c(.025, .975), na.rm = TRUE, names = FALSE))
  gate_q <- as.numeric(quantile(gdraw, c(.025, .975), na.rm = TRUE, names = FALSE))
  gate_mean_by_domain <- rowMeans(pr$gate_G, na.rm = TRUE)
  gate_modal_domain <- which.max(gate_mean_by_domain)
  rows[[ri]] <- data.frame(
    scenario = scen, rep = rep, truth = tru$regime, modal = modal, correct = (modal == tru$regime),
    p_true = as.numeric(probs[tru$regime]), p_modal = as.numeric(max(probs)),
    esr_mean = mean(pr$ESR_link), esr_lo = esr_q[1], esr_hi = esr_q[2],
    esr_cov  = (esr_q[1] <= tru$ESR_link && tru$ESR_link <= esr_q[2]),
    full_mean = mean(pr$full_resp), full_lo = full_q[1], full_hi = full_q[2],
    full_cov = (full_q[1] <= tru$full_resp && tru$full_resp <= full_q[2]),
    gate_domain = tru$gate_domain,
    gate_G_truth = tru$gate_G,
    gate_G_mean = mean(gdraw), gate_G_lo = gate_q[1], gate_G_hi = gate_q[2],
    gate_cov = (gate_q[1] <= tru$gate_G && tru$gate_G <= gate_q[2]),
    gate_modal_domain = gate_modal_domain,
    gate_domain_correct = (gate_modal_domain == tru$gate_domain),
    ff$diag, fit_failed = FALSE)
  cat(sprintf("[%d/%d] rep=%d %-13s truth=%-13s modal=%-13s p_true=%.2f esr_cov=%d\n",
              ri, TOTAL, rep, scen, tru$regime, modal, as.numeric(probs[tru$regime]), rows[[ri]]$esr_cov))
  write.table(rows[[ri]], RAW_PATH, sep = ",", row.names = FALSE,
              col.names = !file.exists(RAW_PATH), append = file.exists(RAW_PATH))
  flush.console()
}
res <- do.call(rbind, rows)
write.csv(res, RAW_PATH, row.names = FALSE)

## ---- summaries ------------------------------------------------------------------------
res_fit <- subset(res, !fit_failed)
recov <- aggregate(cbind(recovery = correct, esr_coverage = esr_cov, full_coverage = full_cov,
                         gate_coverage = gate_cov, mean_divergences = divergences,
                         gate_domain_recovery = gate_domain_correct,
                         max_rhat = max_rhat) ~ scenario,
                   data = res_fit, FUN = function(x) mean(x, na.rm = TRUE))
fail <- aggregate(fit_failed ~ scenario, data = res, FUN = function(x) mean(x, na.rm = TRUE))
names(fail)[names(fail) == "fit_failed"] <- "failure_rate"
recov <- merge(recov, fail, by = "scenario", all = TRUE)
write.csv(recov, file.path(OUTDIR, paste0("fullbayes_recovery_coverage", TAG, ".csv")), row.names = FALSE)

## regime-probability calibration: bin modal-label confidence, compare to empirical accuracy
res_cal <- subset(res, !fit_failed & is.finite(p_modal))
if (nrow(res_cal)) {
  res_cal$bin <- cut(res_cal$p_modal, breaks = seq(0.4, 1.0, by = 0.1), include.lowest = TRUE)
  cal <- do.call(rbind, lapply(split(res_cal, droplevels(res_cal$bin)), function(d)
    data.frame(bin = d$bin[1], mean_confidence = mean(d$p_modal), accuracy = mean(d$correct), n = nrow(d))))
} else {
  cal <- data.frame(bin = character(), mean_confidence = numeric(), accuracy = numeric(), n = integer())
}
write.csv(cal, file.path(OUTDIR, paste0("fullbayes_regime_calibration", TAG, ".csv")), row.names = FALSE)

cat("\n== recovery + coverage ==\n"); print(recov)
cat("\n== regime-probability calibration ==\n"); print(cal)
cat("\nWrote CSVs to", OUTDIR, "\n")

## write booktabs LaTeX tables the manuscript \input's
write_tex <- function(df, file, caption, label) {
  df <- as.data.frame(lapply(df, function(x) if (is.numeric(x)) sprintf("%.3f", x) else as.character(x)),
                      stringsAsFactors = FALSE)
  al <- paste(rep("l", ncol(df)), collapse = "")
  con <- file(file, "w")
  writeLines(c("\\begin{table}[!htbp]", "\\centering",
               sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
               "\\resizebox{\\textwidth}{!}{%",
               sprintf("\\begin{tabular}{%s}", al), "\\toprule",
               paste(paste(gsub("_", "\\\\_", names(df)), collapse = " & "), "\\\\"),
               "\\midrule"), con)
  for (i in seq_len(nrow(df))) writeLines(paste(paste(df[i, ], collapse = " & "), "\\\\"), con)
  writeLines(c("\\bottomrule", "\\end{tabular}", "}", "\\end{table}"), con); close(con)
}
write_tex(recov, file.path(OUTDIR, paste0("fullbayes_recovery_coverage", TAG, ".tex")),
          "Full-Bayes structural recovery and posterior interval coverage by regime.", "tab:fb-recovery")
cal_tex <- cal
if (nrow(cal_tex)) {
  cal_tex$bin <- paste0("$", cal_tex$bin, "$")
  cal_tex$n <- as.character(as.integer(cal_tex$n))
  names(cal_tex) <- c("Posterior probability bin", "Mean confidence", "Accuracy", "n")
}
write_tex(cal_tex, file.path(OUTDIR, paste0("fullbayes_regime_calibration", TAG, ".tex")),
          "Posterior regime-probability calibration: declared-regime confidence versus empirical accuracy.",
          "tab:fb-calibration")
