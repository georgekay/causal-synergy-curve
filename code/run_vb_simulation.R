## =====================================================================================
## Variational-Bayes operating-characteristic simulation for the Causal Synergy Curve.
## Stan ADVI approximation to the same regularized-horseshoe coalition surface used by
## the MCMC reference estimator. This is the scalable simulation engine; NHANES applied
## results remain full MCMC.
## Reports: (1) regime recovery, (2) posterior credible-interval COVERAGE of closure/ESR/gate,
##          (3) posterior REGIME-PROBABILITY CALIBRATION (confidence vs accuracy).
## Truth is g-computed from the true cloglog surface with the same definitions the estimator uses.
##
## Run:  Rscript run_vb_simulation.R          (edit REPS/VB settings/etc. below)
## Needs: cmdstanr (+ CmdStan), and csc_design_readouts.R in the same folder.
## =====================================================================================
suppressPackageStartupMessages({ library(cmdstanr); library(loo) })
source("csc_design_readouts.R")
set.seed(20260711)
options(warn = 1)

## ---- configuration (reduce REPS for a fast first pass; 500 for the paper) -------------
ALL_SCEN  <- c("additive", "complementary", "redundant", "mixed")
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
VB_ALGORITHM <- Sys.getenv("CSC_VB_ALGORITHM", "fullrank")
VB_ITER <- as.integer(Sys.getenv("CSC_VB_ITER", "20000"))
VB_GRAD_SAMPLES <- as.integer(Sys.getenv("CSC_VB_GRAD_SAMPLES", "1"))
VB_ELBO_SAMPLES <- as.integer(Sys.getenv("CSC_VB_ELBO_SAMPLES", "100"))
VB_ADAPT_ITER <- as.integer(Sys.getenv("CSC_VB_ADAPT_ITER", "50"))
VB_TOL_REL_OBJ <- as.numeric(Sys.getenv("CSC_VB_TOL_REL_OBJ", "0.005"))
VB_EVAL_ELBO <- as.integer(Sys.getenv("CSC_VB_EVAL_ELBO", "100"))
USE_PSIS  <- tolower(Sys.getenv("CSC_USE_PSIS", "TRUE")) %in% c("true", "t", "1", "yes", "y")
PSIS_K_OK <- as.numeric(Sys.getenv("CSC_PSIS_K_OK", "0.70"))
ROPE      <- as.numeric(Sys.getenv("CSC_ROPE", "0.10"))  # practical-equivalence region for second differences
GAMMA     <- as.numeric(Sys.getenv("CSC_GAMMA", "1.50")) # link interaction strength; 1.5 clears ROPE at Pr(D_jD_k) ~= .09
MAINS     <- seq(0.55, 0.20, length.out = P)
STAN_FILE <- "../stan/csc_coalition_horseshoe_binary.stan"
OUTDIR    <- Sys.getenv("CSC_OUTDIR", "../simulation_outputs")
JOB_ID    <- Sys.getenv("CSC_JOB_ID", "")
TAG       <- if (nzchar(JOB_ID)) paste0("_", JOB_ID) else ""
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
RAW_PATH  <- file.path(OUTDIR, paste0("vb_simulation_raw", TAG, ".csv"))
if (file.exists(RAW_PATH)) file.remove(RAW_PATH)

mod <- cmdstan_model(STAN_FILE)

## ---- data-generating surfaces (SIGN CONVENTION: +link interaction => REDUNDANT closure) ----
## To plant a COMPLEMENTARY closure regime we use NEGATIVE link interactions, and vice versa.
eta_true_fun <- function(scen) function(D, age) {
  e <- -2.2 + 0.5 * age + as.numeric(D %*% MAINS)
  if (scen == "complementary") e <- e - GAMMA * (D[,1]*D[,2] + D[,3]*D[,4])
  if (scen == "redundant")     e <- e + GAMMA * (D[,1]*D[,2] + D[,3]*D[,4])
  if (scen == "mixed")         e <- e - GAMMA * D[,1]*D[,2] + GAMMA * D[,3]*D[,4]
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
  list(regime = csc_regime(fl$pair_delta, ROPE), ESR_link = fl$ESR,
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

## ---- PSIS correction for ADVI draws --------------------------------------------------
## CmdStan ADVI exposes lp__ (target log density) and lp_approx__ (variational log density).
## Their difference is a valid importance log-ratio up to an irrelevant additive constant.
psis_from_vb_draws <- function(dr) {
  empty <- list(available = FALSE, k = NA_real_, ess = NA_real_, ess_frac = NA_real_,
                weights = rep(1 / nrow(dr), nrow(dr)), message = "")
  if (!USE_PSIS) {
    empty$message <- "PSIS disabled by CSC_USE_PSIS"
    return(empty)
  }
  if (!all(c("lp__", "lp_approx__") %in% names(dr))) {
    empty$message <- "ADVI output lacks lp__ and/or lp_approx__"
    return(empty)
  }
  log_ratios <- as.numeric(dr[["lp__"]] - dr[["lp_approx__"]])
  ok <- is.finite(log_ratios)
  if (sum(ok) < 20) {
    empty$message <- "Too few finite ADVI log-ratios for PSIS"
    return(empty)
  }
  ps <- tryCatch(loo::psis(log_ratios[ok]), error = function(e) e)
  if (inherits(ps, "error")) {
    empty$message <- conditionMessage(ps)
    return(empty)
  }
  w_ok <- as.numeric(weights(ps, log = FALSE, normalize = TRUE))
  w <- rep(0, length(log_ratios))
  w[ok] <- w_ok
  ess <- 1 / sum(w^2)
  list(available = TRUE,
       k = as.numeric(loo::pareto_k_values(ps))[1],
       ess = ess,
       ess_frac = ess / length(log_ratios),
       weights = w / sum(w),
       message = "")
}

wquantile <- function(x, probs, w = NULL) {
  x <- as.numeric(x)
  if (is.null(w)) return(as.numeric(stats::quantile(x, probs, na.rm = TRUE, names = FALSE)))
  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]; w <- w[ok]
  if (!length(x) || sum(w) <= 0) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o] / sum(w)
  cw <- cumsum(w)
  as.numeric(approx(cw, x, xout = probs, method = "linear", ties = "ordered", rule = 2)$y)
}

weighted_regime_probs <- function(regime, w = NULL) {
  lev <- c("additive", "complementary", "redundant", "mixed")
  if (is.null(w)) w <- rep(1 / length(regime), length(regime))
  w <- as.numeric(w); w <- w / sum(w)
  pr <- sapply(lev, function(z) sum(w[regime == z], na.rm = TRUE))
  pr / sum(pr)
}

summarize_readout <- function(pr, tru, w = NULL) {
  probs <- weighted_regime_probs(pr$regime, w)
  modal <- names(which.max(probs))
  gdraw <- pr$gate_G[tru$gate_domain, ]
  esr_q <- wquantile(pr$ESR_link, c(.025, .975), w)
  full_q <- wquantile(pr$full_resp, c(.025, .975), w)
  gate_q <- wquantile(gdraw, c(.025, .975), w)
  data.frame(
    modal = modal,
    correct = (modal == tru$regime),
    p_true = as.numeric(probs[tru$regime]),
    p_modal = as.numeric(max(probs)),
    esr_mean = if (is.null(w)) mean(pr$ESR_link) else sum(as.numeric(w) * pr$ESR_link),
    esr_lo = esr_q[1], esr_hi = esr_q[2],
    esr_cov = (esr_q[1] <= tru$ESR_link && tru$ESR_link <= esr_q[2]),
    full_mean = if (is.null(w)) mean(pr$full_resp) else sum(as.numeric(w) * pr$full_resp),
    full_lo = full_q[1], full_hi = full_q[2],
    full_cov = (full_q[1] <= tru$full_resp && tru$full_resp <= full_q[2]),
    gate_G_mean = if (is.null(w)) mean(gdraw) else sum(as.numeric(w) * gdraw),
    gate_G_lo = gate_q[1], gate_G_hi = gate_q[2],
    gate_cov = (gate_q[1] <= tru$gate_G && tru$gate_G <= gate_q[2])
  )
}

fit_one <- function(dat, seed) {
  meta <- csc_build_design(dat$D, dat$W, ORDER)
  Pc <- ncol(meta$Xc); p0 <- 2
  sdat <- list(N = N, Pb = ncol(meta$Xb), Pc = Pc, Xb = meta$Xb, Xc = meta$Xc, y = as.integer(dat$y),
               scale_global = p0 / ((Pc - p0) * sqrt(N)), slab_scale = 2, slab_df = 4, beta_scale = 2.5)
  fit <- mod$variational(data = sdat, seed = seed, refresh = 0, init = 0,
                         algorithm = VB_ALGORITHM, iter = VB_ITER,
                         grad_samples = VB_GRAD_SAMPLES, elbo_samples = VB_ELBO_SAMPLES,
                         eta = 1, adapt_engaged = TRUE, adapt_iter = VB_ADAPT_ITER,
                         tol_rel_obj = VB_TOL_REL_OBJ, eval_elbo = VB_EVAL_ELBO,
                         output_samples = DRAWS_KEEP,
                         show_messages = FALSE, show_exceptions = FALSE)
  diag <- data.frame(
    vb_algorithm = VB_ALGORITHM,
    vb_iter = VB_ITER,
    vb_output_samples = DRAWS_KEEP
  )
  dr <- fit$draws(format = "df")
  idx <- sort(sample(nrow(dr), min(DRAWS_KEEP, nrow(dr))))
  psis <- psis_from_vb_draws(dr[idx, , drop = FALSE])
  gb <- function(pat) as.matrix(dr[idx, grep(pat, colnames(dr)), drop = FALSE])
  post <- list(alpha0 = dr[[grep("^alpha0$", colnames(dr))]][idx],
               beta_b = gb("^beta_b\\["), theta = gb("^theta\\["))
  list(meta = meta, post = post, diag = diag, psis = psis)
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
  ff  <- tryCatch(
    fit_one(dat, seed = 500000 + 1000 * match(scen, ALL_SCEN) + rep),
    error = function(e) structure(list(message = conditionMessage(e)), class = "vb_fit_failed")
  )
  if (inherits(ff, "vb_fit_failed")) {
    rows[[ri]] <- data.frame(
      scenario = scen, rep = rep, truth = tru$regime, modal = "failed", correct = FALSE,
      p_true = NA_real_, p_modal = NA_real_,
      esr_mean = NA_real_, esr_lo = NA_real_, esr_hi = NA_real_, esr_cov = NA,
      full_mean = NA_real_, full_lo = NA_real_, full_hi = NA_real_, full_cov = NA,
      gate_domain = tru$gate_domain, gate_G_truth = tru$gate_G,
      gate_G_mean = NA_real_, gate_G_lo = NA_real_, gate_G_hi = NA_real_, gate_cov = NA,
      psis_available = FALSE, psis_k = NA_real_, psis_ess_frac = NA_real_,
      psis_ok = FALSE, psis_message = "",
      modal_psis = "failed", correct_psis = FALSE,
      p_true_psis = NA_real_, p_modal_psis = NA_real_,
      esr_mean_psis = NA_real_, esr_lo_psis = NA_real_, esr_hi_psis = NA_real_,
      esr_cov_psis = NA,
      full_mean_psis = NA_real_, full_lo_psis = NA_real_, full_hi_psis = NA_real_,
      full_cov_psis = NA,
      gate_G_mean_psis = NA_real_, gate_G_lo_psis = NA_real_, gate_G_hi_psis = NA_real_,
      gate_cov_psis = NA,
      vb_algorithm = VB_ALGORITHM, vb_iter = VB_ITER, vb_output_samples = DRAWS_KEEP,
      fit_failed = TRUE, fit_message = ff$message
    )
    cat(sprintf("[%d/%d] rep=%d %-13s truth=%-13s modal=failed fit error: %s\n",
                ri, TOTAL, rep, scen, tru$regime, ff$message))
    write.table(rows[[ri]], RAW_PATH, sep = ",", row.names = FALSE,
                col.names = !file.exists(RAW_PATH), append = file.exists(RAW_PATH))
    flush.console()
    next
  }
  pr  <- posterior_readouts(dat, ff$meta, ff$post)
  sum_unw <- summarize_readout(pr, tru)
  sum_psis <- if (isTRUE(ff$psis$available)) summarize_readout(pr, tru, ff$psis$weights) else
    data.frame(modal = NA_character_, correct = NA, p_true = NA_real_, p_modal = NA_real_,
               esr_mean = NA_real_, esr_lo = NA_real_, esr_hi = NA_real_, esr_cov = NA,
               full_mean = NA_real_, full_lo = NA_real_, full_hi = NA_real_, full_cov = NA,
               gate_G_mean = NA_real_, gate_G_lo = NA_real_, gate_G_hi = NA_real_,
               gate_cov = NA)
  rows[[ri]] <- data.frame(
    scenario = scen, rep = rep, truth = tru$regime,
    modal = sum_unw$modal, correct = sum_unw$correct,
    p_true = sum_unw$p_true, p_modal = sum_unw$p_modal,
    esr_mean = sum_unw$esr_mean, esr_lo = sum_unw$esr_lo, esr_hi = sum_unw$esr_hi,
    esr_cov  = sum_unw$esr_cov,
    full_mean = sum_unw$full_mean, full_lo = sum_unw$full_lo, full_hi = sum_unw$full_hi,
    full_cov = sum_unw$full_cov,
    gate_domain = tru$gate_domain,
    gate_G_truth = tru$gate_G,
    gate_G_mean = sum_unw$gate_G_mean,
    gate_G_lo = sum_unw$gate_G_lo, gate_G_hi = sum_unw$gate_G_hi,
    gate_cov = sum_unw$gate_cov,
    psis_available = ff$psis$available,
    psis_k = ff$psis$k,
    psis_ess_frac = ff$psis$ess_frac,
    psis_ok = isTRUE(ff$psis$available) && is.finite(ff$psis$k) && ff$psis$k < PSIS_K_OK,
    psis_message = ff$psis$message,
    modal_psis = sum_psis$modal,
    correct_psis = sum_psis$correct,
    p_true_psis = sum_psis$p_true,
    p_modal_psis = sum_psis$p_modal,
    esr_mean_psis = sum_psis$esr_mean,
    esr_lo_psis = sum_psis$esr_lo,
    esr_hi_psis = sum_psis$esr_hi,
    esr_cov_psis = sum_psis$esr_cov,
    full_mean_psis = sum_psis$full_mean,
    full_lo_psis = sum_psis$full_lo,
    full_hi_psis = sum_psis$full_hi,
    full_cov_psis = sum_psis$full_cov,
    gate_G_mean_psis = sum_psis$gate_G_mean,
    gate_G_lo_psis = sum_psis$gate_G_lo,
    gate_G_hi_psis = sum_psis$gate_G_hi,
    gate_cov_psis = sum_psis$gate_cov,
    ff$diag, fit_failed = FALSE, fit_message = "")
  cat(sprintf("[%d/%d] rep=%d %-13s truth=%-13s modal=%-13s p_true=%.2f esr_cov=%d\n",
              ri, TOTAL, rep, scen, tru$regime, sum_unw$modal, sum_unw$p_true, rows[[ri]]$esr_cov))
  write.table(rows[[ri]], RAW_PATH, sep = ",", row.names = FALSE,
              col.names = !file.exists(RAW_PATH), append = file.exists(RAW_PATH))
  flush.console()
}
res <- do.call(rbind, rows)
write.csv(res, RAW_PATH, row.names = FALSE)

## ---- summaries ------------------------------------------------------------------------
res_fit <- subset(res, !fit_failed)
if (nrow(res_fit)) {
  recov <- aggregate(cbind(recovery = correct, esr_coverage = esr_cov, full_coverage = full_cov,
                           gate_coverage = gate_cov) ~ scenario,
                     data = res_fit, FUN = function(x) mean(x, na.rm = TRUE))
} else {
  recov <- data.frame(scenario = unique(res$scenario), recovery = NA_real_,
                      esr_coverage = NA_real_, full_coverage = NA_real_,
                      gate_coverage = NA_real_)
}
fail <- aggregate(fit_failed ~ scenario, data = res, FUN = function(x) mean(x, na.rm = TRUE))
names(fail)[names(fail) == "fit_failed"] <- "failure_rate"
if (nrow(res_fit)) {
  psis_sum <- aggregate(cbind(psis_available = psis_available, psis_ok = psis_ok,
                              psis_k = psis_k, psis_ess_frac = psis_ess_frac) ~ scenario,
                        data = res_fit, FUN = function(x) mean(x, na.rm = TRUE))
} else {
  psis_sum <- data.frame(scenario = unique(res$scenario), psis_available = NA_real_,
                         psis_ok = NA_real_, psis_k = NA_real_, psis_ess_frac = NA_real_)
}
recov <- Reduce(function(a, b) merge(a, b, by = "scenario", all = TRUE),
                list(recov, fail, psis_sum))
res_psis <- subset(res_fit, psis_available)
if (nrow(res_psis)) {
  recov_psis <- aggregate(cbind(recovery_psis = correct_psis, esr_coverage_psis = esr_cov_psis,
                                full_coverage_psis = full_cov_psis,
                                gate_coverage_psis = gate_cov_psis) ~ scenario,
                          data = res_psis, FUN = function(x) mean(x, na.rm = TRUE))
  recov <- merge(recov, recov_psis, by = "scenario", all = TRUE)
}
write.csv(recov, file.path(OUTDIR, paste0("vb_recovery_coverage", TAG, ".csv")), row.names = FALSE)

## regime-probability calibration: bin modal-label confidence, compare to empirical accuracy
res_ok <- subset(res, !fit_failed & is.finite(p_modal))
if (nrow(res_ok)) {
  res_ok$bin <- cut(res_ok$p_modal, breaks = seq(0.4, 1.0, by = 0.1), include.lowest = TRUE)
  cal <- do.call(rbind, lapply(split(res_ok, droplevels(res_ok$bin)), function(d)
    data.frame(bin = d$bin[1], mean_confidence = mean(d$p_modal), accuracy = mean(d$correct), n = nrow(d))))
} else {
  cal <- data.frame(bin = character(), mean_confidence = numeric(), accuracy = numeric(), n = integer())
}
write.csv(cal, file.path(OUTDIR, paste0("vb_regime_calibration", TAG, ".csv")), row.names = FALSE)

res_ok_psis <- subset(res, !fit_failed & psis_available & is.finite(p_modal_psis))
if (nrow(res_ok_psis)) {
  res_ok_psis$bin <- cut(res_ok_psis$p_modal_psis, breaks = seq(0.4, 1.0, by = 0.1), include.lowest = TRUE)
  cal_psis <- do.call(rbind, lapply(split(res_ok_psis, droplevels(res_ok_psis$bin)), function(d)
    data.frame(bin = d$bin[1], mean_confidence = mean(d$p_modal_psis),
               accuracy = mean(d$correct_psis), n = nrow(d))))
} else {
  cal_psis <- data.frame(bin = character(), mean_confidence = numeric(), accuracy = numeric(), n = integer())
}
write.csv(cal_psis, file.path(OUTDIR, paste0("vb_regime_calibration_psis", TAG, ".csv")), row.names = FALSE)

cat("\n== recovery + coverage ==\n"); print(recov)
cat("\n== regime-probability calibration ==\n"); print(cal)
cat("\n== PSIS-weighted regime-probability calibration ==\n"); print(cal_psis)
cat("\nWrote CSVs to", OUTDIR, "\n")

## write booktabs LaTeX tables the manuscript \input's
write_tex <- function(df, file, caption, label) {
  df <- as.data.frame(lapply(df, function(x) if (is.numeric(x)) sprintf("%.3f", x) else as.character(x)),
                      stringsAsFactors = FALSE)
  al <- paste(rep("l", ncol(df)), collapse = "")
  con <- file(file, "w")
  writeLines(c("\\begin{table}[!htbp]", "\\centering",
               sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
               sprintf("\\begin{tabular}{%s}", al), "\\toprule",
               paste(paste(gsub("_", "\\\\_", names(df)), collapse = " & "), "\\\\"),
               "\\midrule"), con)
  for (i in seq_len(nrow(df))) writeLines(paste(paste(df[i, ], collapse = " & "), "\\\\"), con)
  writeLines(c("\\bottomrule", "\\end{tabular}", "\\end{table}"), con); close(con)
}
write_tex(recov, file.path(OUTDIR, paste0("vb_recovery_coverage", TAG, ".tex")),
          "Variational-Bayes structural recovery and posterior interval coverage by regime.", "tab:vb-recovery")
write_tex(cal, file.path(OUTDIR, paste0("vb_regime_calibration", TAG, ".tex")),
          "Variational posterior regime-probability calibration: declared-regime confidence versus empirical accuracy.",
          "tab:vb-calibration")
write_tex(cal_psis, file.path(OUTDIR, paste0("vb_regime_calibration_psis", TAG, ".tex")),
          "PSIS-weighted variational regime-probability calibration: declared-regime confidence versus empirical accuracy.",
          "tab:vb-calibration-psis")
