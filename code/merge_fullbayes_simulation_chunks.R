## Merge chunked Full-Bayes simulation outputs into the production CSV and LaTeX tables
## expected by manuscript/simulation_fullbayes_revised.tex.
suppressPackageStartupMessages({ })

OUTDIR <- Sys.getenv("CSC_OUTDIR", "../simulation_outputs")
if (!dir.exists(OUTDIR)) stop("Output directory does not exist: ", OUTDIR)

raw_files <- list.files(
  OUTDIR,
  pattern = "^fullbayes_simulation_raw_.*chunk.+\\.csv$",
  full.names = TRUE
)
if (!length(raw_files)) stop("No fullbayes_simulation_raw*.csv files found in ", OUTDIR)

res <- do.call(rbind, lapply(raw_files, read.csv, stringsAsFactors = FALSE))
if (!nrow(res)) stop("Merged simulation output has zero rows.")

bool_cols <- c("correct", "esr_cov", "full_cov", "gate_cov", "gate_domain_correct", "fit_failed")
as_bool_na <- function(x) {
  y <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1")] <- TRUE
  out[y %in% c("false", "f", "0")] <- FALSE
  out
}
for (cc in intersect(bool_cols, names(res))) res[[cc]] <- as_bool_na(res[[cc]])
if (!"fit_failed" %in% names(res)) res$fit_failed <- FALSE
res$fit_failed[is.na(res$fit_failed)] <- FALSE

res_fit <- subset(res, !fit_failed)
recov <- aggregate(cbind(recovery = correct, esr_coverage = esr_cov, full_coverage = full_cov,
                         gate_coverage = gate_cov, mean_divergences = divergences,
                         gate_domain_recovery = gate_domain_correct,
                         max_rhat = max_rhat) ~ scenario,
                   data = res_fit, FUN = function(x) mean(x, na.rm = TRUE))
fail <- aggregate(fit_failed ~ scenario, data = res, FUN = function(x) mean(x, na.rm = TRUE))
names(fail)[names(fail) == "fit_failed"] <- "failure_rate"
recov <- merge(recov, fail, by = "scenario", all = TRUE)

res_cal <- subset(res, !fit_failed & is.finite(p_modal))
if (nrow(res_cal)) {
  res_cal$bin <- cut(res_cal$p_modal, breaks = seq(0.4, 1.0, by = 0.1), include.lowest = TRUE)
  cal <- do.call(rbind, lapply(split(res_cal, droplevels(res_cal$bin)), function(d) {
    data.frame(bin = as.character(d$bin[1]),
               mean_confidence = mean(d$p_modal),
               accuracy = mean(d$correct),
               n = nrow(d))
  }))
} else {
  cal <- data.frame(bin = character(), mean_confidence = numeric(), accuracy = numeric(), n = integer())
}
rownames(cal) <- NULL

write.csv(res, file.path(OUTDIR, "fullbayes_simulation_raw.csv"), row.names = FALSE)
write.csv(recov, file.path(OUTDIR, "fullbayes_recovery_coverage.csv"), row.names = FALSE)
write.csv(cal, file.path(OUTDIR, "fullbayes_regime_calibration.csv"), row.names = FALSE)

write_tex <- function(df, file, caption, label) {
  df <- as.data.frame(lapply(df, function(x) if (is.numeric(x)) sprintf("%.3f", x) else as.character(x)),
                      stringsAsFactors = FALSE)
  al <- paste(rep("l", ncol(df)), collapse = "")
  con <- file(file, "w")
  on.exit(close(con), add = TRUE)
  writeLines(c("\\begin{table}[!htbp]", "\\centering",
               sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
               "\\resizebox{\\textwidth}{!}{%",
               sprintf("\\begin{tabular}{%s}", al), "\\toprule",
               paste(paste(gsub("_", "\\\\_", names(df)), collapse = " & "), "\\\\"),
               "\\midrule"), con)
  for (i in seq_len(nrow(df))) writeLines(paste(paste(df[i, ], collapse = " & "), "\\\\"), con)
  writeLines(c("\\bottomrule", "\\end{tabular}", "}", "\\end{table}"), con)
}

write_tex(recov, file.path(OUTDIR, "fullbayes_recovery_coverage.tex"),
          "Full-Bayes structural recovery and posterior interval coverage by regime.", "tab:fb-recovery")
cal_tex <- cal
if (nrow(cal_tex)) {
  cal_tex$bin <- paste0("$", cal_tex$bin, "$")
  cal_tex$n <- as.character(as.integer(cal_tex$n))
  names(cal_tex) <- c("Posterior probability bin", "Mean confidence", "Accuracy", "n")
}
write_tex(cal_tex, file.path(OUTDIR, "fullbayes_regime_calibration.tex"),
          "Posterior regime-probability calibration: declared-regime confidence versus empirical accuracy.",
          "tab:fb-calibration")

cat("Merged", length(raw_files), "chunk files and", nrow(res), "replication rows.\n")
cat("Wrote production outputs to:", normalizePath(OUTDIR, winslash = "\\", mustWork = FALSE), "\n")
