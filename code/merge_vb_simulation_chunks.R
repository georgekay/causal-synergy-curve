## Merge chunked variational-Bayes simulation outputs into the production CSV and LaTeX tables
## expected by manuscript/simulation_fullbayes_revised.tex.
suppressPackageStartupMessages({ })

OUTDIR <- Sys.getenv("CSC_OUTDIR", "../simulation_outputs")
if (!dir.exists(OUTDIR)) stop("Output directory does not exist: ", OUTDIR)

raw_files <- list.files(
  OUTDIR,
  pattern = "^vb_simulation_raw_.*chunk.+\\.csv$",
  full.names = TRUE
)
if (!length(raw_files)) stop("No vb_simulation_raw*.csv files found in ", OUTDIR)

res <- do.call(rbind, lapply(raw_files, read.csv, stringsAsFactors = FALSE))
if (!nrow(res)) stop("Merged simulation output has zero rows.")

bool_cols <- c("correct", "esr_cov", "full_cov", "gate_cov", "fit_failed",
               "psis_available", "psis_ok", "correct_psis", "esr_cov_psis",
               "full_cov_psis", "gate_cov_psis")
as_bool_na <- function(x) {
  y <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1")] <- TRUE
  out[y %in% c("false", "f", "0")] <- FALSE
  out
}
for (cc in intersect(bool_cols, names(res))) res[[cc]] <- as_bool_na(res[[cc]])
res$fit_failed[is.na(res$fit_failed)] <- FALSE
numeric_cols <- c("p_true", "p_modal", "p_true_psis", "p_modal_psis",
                  "psis_k", "psis_ess_frac",
                  "esr_mean", "esr_lo", "esr_hi", "full_mean", "full_lo", "full_hi",
                  "esr_mean_psis", "esr_lo_psis", "esr_hi_psis",
                  "full_mean_psis", "full_lo_psis", "full_hi_psis",
                  "gate_mean", "gate_lo", "gate_hi",
                  "gate_G_mean", "gate_G_lo", "gate_G_hi",
                  "gate_G_mean_psis", "gate_G_lo_psis", "gate_G_hi_psis",
                  "true_esr", "true_full", "true_gate", "rep")
for (cc in intersect(numeric_cols, names(res))) res[[cc]] <- as.numeric(res[[cc]])

res_ok <- subset(res, !fit_failed)
recov <- aggregate(cbind(recovery = correct, esr_coverage = esr_cov, full_coverage = full_cov,
                         gate_coverage = gate_cov) ~ scenario,
                   data = res_ok, FUN = function(x) mean(x, na.rm = TRUE))
fail <- aggregate(fit_failed ~ scenario, data = res, FUN = function(x) mean(x, na.rm = TRUE))
names(fail)[names(fail) == "fit_failed"] <- "failure_rate"
recov <- merge(recov, fail, by = "scenario", all = TRUE)

if (all(c("psis_available", "psis_ok", "psis_k", "psis_ess_frac") %in% names(res_ok))) {
  psis_rate <- aggregate(cbind(psis_available = psis_available, psis_ok = psis_ok) ~ scenario,
                         data = res_ok, FUN = function(x) mean(x, na.rm = TRUE))
  if (any(is.finite(res_ok$psis_k)) || any(is.finite(res_ok$psis_ess_frac))) {
    psis_moments <- aggregate(cbind(psis_k = psis_k, psis_ess_frac = psis_ess_frac) ~ scenario,
                              data = subset(res_ok, is.finite(psis_k) | is.finite(psis_ess_frac)),
                              FUN = function(x) mean(x, na.rm = TRUE))
    psis_sum <- merge(psis_rate, psis_moments, by = "scenario", all = TRUE)
  } else {
    psis_sum <- psis_rate
    psis_sum$psis_k <- NA_real_
    psis_sum$psis_ess_frac <- NA_real_
  }
  recov <- merge(recov, psis_sum, by = "scenario", all = TRUE)
}
if (all(c("psis_available", "correct_psis", "esr_cov_psis", "full_cov_psis",
          "gate_cov_psis") %in% names(res_ok)) && any(res_ok$psis_available, na.rm = TRUE)) {
  res_psis <- subset(res_ok, psis_available)
  recov_psis <- aggregate(cbind(recovery_psis = correct_psis,
                                esr_coverage_psis = esr_cov_psis,
                                full_coverage_psis = full_cov_psis,
                                gate_coverage_psis = gate_cov_psis) ~ scenario,
                          data = res_psis, FUN = function(x) mean(x, na.rm = TRUE))
  recov <- merge(recov, recov_psis, by = "scenario", all = TRUE)
}

res_ok <- subset(res, !fit_failed & is.finite(p_modal))
if (nrow(res_ok)) {
  res_ok$bin <- cut(res_ok$p_modal, breaks = seq(0.4, 1.0, by = 0.1), include.lowest = TRUE)
  cal <- do.call(rbind, lapply(split(res_ok, droplevels(res_ok$bin)), function(d) {
    data.frame(bin = as.character(d$bin[1]),
               mean_confidence = mean(d$p_modal),
               accuracy = mean(d$correct),
               n = nrow(d))
  }))
} else {
  cal <- data.frame(bin = character(), mean_confidence = numeric(), accuracy = numeric(), n = integer())
}
rownames(cal) <- NULL

write.csv(res, file.path(OUTDIR, "vb_simulation_raw.csv"), row.names = FALSE)
write.csv(recov, file.path(OUTDIR, "vb_recovery_coverage.csv"), row.names = FALSE)
write.csv(cal, file.path(OUTDIR, "vb_regime_calibration.csv"), row.names = FALSE)

if (all(c("psis_available", "p_modal_psis", "correct_psis") %in% names(res))) {
  res_ok_psis <- subset(res, !fit_failed & psis_available & is.finite(p_modal_psis))
  if (nrow(res_ok_psis)) {
    res_ok_psis$bin <- cut(res_ok_psis$p_modal_psis, breaks = seq(0.4, 1.0, by = 0.1), include.lowest = TRUE)
    cal_psis <- do.call(rbind, lapply(split(res_ok_psis, droplevels(res_ok_psis$bin)), function(d) {
      data.frame(bin = as.character(d$bin[1]),
                 mean_confidence = mean(d$p_modal_psis),
                 accuracy = mean(d$correct_psis),
                 n = nrow(d))
    }))
  } else {
    cal_psis <- data.frame(bin = character(), mean_confidence = numeric(), accuracy = numeric(), n = integer())
  }
  rownames(cal_psis) <- NULL
  write.csv(cal_psis, file.path(OUTDIR, "vb_regime_calibration_psis.csv"), row.names = FALSE)
}

write_tex <- function(df, file, caption, label) {
  df <- as.data.frame(lapply(df, function(x) if (is.numeric(x)) sprintf("%.3f", x) else as.character(x)),
                      stringsAsFactors = FALSE)
  al <- paste(rep("l", ncol(df)), collapse = "")
  con <- file(file, "w")
  on.exit(close(con), add = TRUE)
  writeLines(c("\\begin{table}[!htbp]", "\\centering",
               sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
               sprintf("\\begin{tabular}{%s}", al), "\\toprule",
               paste(paste(gsub("_", "\\\\_", names(df)), collapse = " & "), "\\\\"),
               "\\midrule"), con)
  for (i in seq_len(nrow(df))) writeLines(paste(paste(df[i, ], collapse = " & "), "\\\\"), con)
  writeLines(c("\\bottomrule", "\\end{tabular}", "\\end{table}"), con)
}

write_tex(recov, file.path(OUTDIR, "vb_recovery_coverage.tex"),
          "Variational-Bayes structural recovery and posterior interval coverage by regime.", "tab:vb-recovery")
write_tex(cal, file.path(OUTDIR, "vb_regime_calibration.tex"),
          "Variational posterior regime-probability calibration: declared-regime confidence versus empirical accuracy.",
          "tab:vb-calibration")
if (exists("cal_psis")) {
  write_tex(cal_psis, file.path(OUTDIR, "vb_regime_calibration_psis.tex"),
            "PSIS-weighted variational regime-probability calibration: declared-regime confidence versus empirical accuracy.",
            "tab:vb-calibration-psis")
}

cat("Merged", length(raw_files), "chunk files and", nrow(res), "replication rows.\n")
cat("Wrote production outputs to:", normalizePath(OUTDIR, winslash = "\\", mustWork = FALSE), "\n")
