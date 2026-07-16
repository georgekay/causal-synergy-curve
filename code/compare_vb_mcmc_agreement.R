## Compare a stratified VB simulation panel with a matched MCMC panel.
##
## Expected inputs:
##   CSC_VB_DIR   = directory containing vb_simulation_raw.csv or vb_simulation_raw_chunk*.csv
##   CSC_MCMC_DIR = directory containing fullbayes_simulation_raw.csv or fullbayes_simulation_raw_chunk*.csv
##   CSC_AGREE_OUTDIR = optional output directory; defaults to dirname(CSC_VB_DIR)/agreement_compare

read_raw <- function(dir, stem) {
  merged <- file.path(dir, paste0(stem, ".csv"))
  if (file.exists(merged)) return(read.csv(merged, stringsAsFactors = FALSE))
  files <- list.files(dir, pattern = paste0("^", stem, "_chunk.+\\.csv$"), full.names = TRUE)
  if (!length(files)) stop("No raw files found for ", stem, " in ", dir)
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  y <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1")] <- TRUE
  out[y %in% c("false", "f", "0")] <- FALSE
  out
}

num <- function(x) suppressWarnings(as.numeric(x))

write_tex <- function(df, file, caption, label) {
  df <- as.data.frame(lapply(df, function(x) if (is.numeric(x)) sprintf("%.3f", x) else as.character(x)),
                      stringsAsFactors = FALSE)
  con <- file(file, "w")
  on.exit(close(con), add = TRUE)
  writeLines(c("\\begin{table}[!htbp]", "\\centering",
               sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
               sprintf("\\begin{tabular}{%s}", paste(rep("l", ncol(df)), collapse = "")),
               "\\toprule",
               paste(paste(gsub("_", "\\\\_", names(df)), collapse = " & "), "\\\\"),
               "\\midrule"), con)
  for (i in seq_len(nrow(df))) writeLines(paste(paste(df[i, ], collapse = " & "), "\\\\"), con)
  writeLines(c("\\bottomrule", "\\end{tabular}", "\\end{table}"), con)
}

vb_dir <- Sys.getenv("CSC_VB_DIR")
mcmc_dir <- Sys.getenv("CSC_MCMC_DIR")
if (!nzchar(vb_dir) || !dir.exists(vb_dir)) stop("Set CSC_VB_DIR to the VB panel output directory.")
if (!nzchar(mcmc_dir) || !dir.exists(mcmc_dir)) stop("Set CSC_MCMC_DIR to the MCMC panel output directory.")
outdir <- Sys.getenv("CSC_AGREE_OUTDIR")
if (!nzchar(outdir)) outdir <- file.path(dirname(vb_dir), "agreement_compare")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

vb <- read_raw(vb_dir, "vb_simulation_raw")
mc <- read_raw(mcmc_dir, "fullbayes_simulation_raw")
if (!"rep" %in% names(vb) || !"rep" %in% names(mc)) {
  stop("Both VB and MCMC raw files must contain a 'rep' column. Re-run with the updated scripts.")
}
if (!"fit_failed" %in% names(vb)) vb$fit_failed <- FALSE
if (!"fit_failed" %in% names(mc)) mc$fit_failed <- FALSE
vb$fit_failed <- as_bool(vb$fit_failed)
mc$fit_failed <- as_bool(mc$fit_failed)
vb <- subset(vb, !fit_failed)
mc <- subset(mc, !fit_failed)

keys <- c("scenario", "rep", "truth")
dat <- merge(vb, mc, by = keys, suffixes = c("_vb", "_mcmc"))
if (!nrow(dat)) stop("No matched VB/MCMC rows after joining on scenario, rep, truth.")

if (!"psis_available" %in% names(dat)) dat$psis_available <- FALSE
dat$psis_available <- as_bool(dat$psis_available)
if (!"psis_ok" %in% names(dat)) dat$psis_ok <- NA
dat$psis_ok <- as_bool(dat$psis_ok)

has_psis_means <- all(c("esr_mean_psis", "full_mean_psis", "gate_G_mean_psis", "modal_psis") %in% names(dat))
dat$modal_vb_used <- dat$modal_vb
dat$esr_mean_vb_used <- num(dat$esr_mean_vb)
dat$full_mean_vb_used <- num(dat$full_mean_vb)
dat$gate_G_mean_vb_used <- num(dat$gate_G_mean_vb)
if (has_psis_means) {
  use <- dat$psis_available & is.finite(num(dat$esr_mean_psis))
  dat$modal_vb_used[use] <- dat$modal_psis[use]
  dat$esr_mean_vb_used[use] <- num(dat$esr_mean_psis[use])
  dat$full_mean_vb_used[use] <- num(dat$full_mean_psis[use])
  dat$gate_G_mean_vb_used[use] <- num(dat$gate_G_mean_psis[use])
}

dat$esr_mean_mcmc <- num(dat$esr_mean_mcmc)
dat$full_mean_mcmc <- num(dat$full_mean_mcmc)
dat$gate_G_mean_mcmc <- num(dat$gate_G_mean_mcmc)
dat$same_modal <- dat$modal_vb == dat$modal_mcmc
dat$same_modal_used <- dat$modal_vb_used == dat$modal_mcmc
dat$abs_esr_diff <- abs(dat$esr_mean_vb - dat$esr_mean_mcmc)
dat$abs_full_diff <- abs(dat$full_mean_vb - dat$full_mean_mcmc)
dat$abs_gate_diff <- abs(dat$gate_G_mean_vb - dat$gate_G_mean_mcmc)
dat$abs_esr_diff_used <- abs(dat$esr_mean_vb_used - dat$esr_mean_mcmc)
dat$abs_full_diff_used <- abs(dat$full_mean_vb_used - dat$full_mean_mcmc)
dat$abs_gate_diff_used <- abs(dat$gate_G_mean_vb_used - dat$gate_G_mean_mcmc)

summ <- do.call(rbind, lapply(split(dat, dat$scenario), function(d) {
  data.frame(
    scenario = d$scenario[1],
    n_matched = nrow(d),
    modal_agreement = mean(d$same_modal, na.rm = TRUE),
    modal_agreement_psis_if_available = mean(d$same_modal_used, na.rm = TRUE),
    mean_abs_ESR_difference = mean(d$abs_esr_diff, na.rm = TRUE),
    mean_abs_ESR_difference_psis_if_available = mean(d$abs_esr_diff_used, na.rm = TRUE),
    mean_abs_full_closure_difference = mean(d$abs_full_diff, na.rm = TRUE),
    mean_abs_full_closure_difference_psis_if_available = mean(d$abs_full_diff_used, na.rm = TRUE),
    mean_abs_gate_difference = mean(d$abs_gate_diff, na.rm = TRUE),
    mean_abs_gate_difference_psis_if_available = mean(d$abs_gate_diff_used, na.rm = TRUE),
    VB_recovery = mean(as_bool(d$correct_vb), na.rm = TRUE),
    MCMC_recovery = mean(as_bool(d$correct_mcmc), na.rm = TRUE),
    PSIS_available = mean(d$psis_available, na.rm = TRUE),
    PSIS_ok_k_lt_0_7 = mean(d$psis_ok, na.rm = TRUE),
    mean_Pareto_k = mean(num(d$psis_k), na.rm = TRUE)
  )
}))
rownames(summ) <- NULL

overall <- data.frame(
  scenario = "overall",
  n_matched = nrow(dat),
  modal_agreement = mean(dat$same_modal, na.rm = TRUE),
  modal_agreement_psis_if_available = mean(dat$same_modal_used, na.rm = TRUE),
  mean_abs_ESR_difference = mean(dat$abs_esr_diff, na.rm = TRUE),
  mean_abs_ESR_difference_psis_if_available = mean(dat$abs_esr_diff_used, na.rm = TRUE),
  mean_abs_full_closure_difference = mean(dat$abs_full_diff, na.rm = TRUE),
  mean_abs_full_closure_difference_psis_if_available = mean(dat$abs_full_diff_used, na.rm = TRUE),
  mean_abs_gate_difference = mean(dat$abs_gate_diff, na.rm = TRUE),
  mean_abs_gate_difference_psis_if_available = mean(dat$abs_gate_diff_used, na.rm = TRUE),
  VB_recovery = mean(as_bool(dat$correct_vb), na.rm = TRUE),
  MCMC_recovery = mean(as_bool(dat$correct_mcmc), na.rm = TRUE),
  PSIS_available = mean(dat$psis_available, na.rm = TRUE),
  PSIS_ok_k_lt_0_7 = mean(dat$psis_ok, na.rm = TRUE),
  mean_Pareto_k = mean(num(dat$psis_k), na.rm = TRUE)
)
summ <- rbind(summ, overall)

write.csv(dat, file.path(outdir, "vb_mcmc_agreement_matched_rows.csv"), row.names = FALSE)
write.csv(summ, file.path(outdir, "vb_mcmc_agreement_summary.csv"), row.names = FALSE)
write_tex(summ, file.path(outdir, "vb_mcmc_agreement_summary.tex"),
          "Agreement between variational-Bayes and MCMC estimators on the stratified simulation panel.",
          "tab:vb-mcmc-agreement")

cat("Matched rows:", nrow(dat), "\n")
cat("Wrote agreement outputs to:", normalizePath(outdir, winslash = "\\", mustWork = FALSE), "\n")
print(summ)
