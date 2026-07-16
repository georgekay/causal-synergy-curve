## Build applied FullBayes tables and figures from generated posterior g-computation CSVs.
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

ROOT <- normalizePath("..", winslash = "/", mustWork = TRUE)
OUT <- file.path(ROOT, "applied_outputs")
TAB <- file.path(ROOT, "tables")
FIG <- file.path(ROOT, "figures")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

cohorts <- c("pooled", "black", "white")
domain_map <- c(
  D_employment = "Employment",
  D_income = "Income",
  D_food_security = "Food security",
  D_education = "Education",
  D_access_to_care = "Access to care",
  D_health_insurance = "Health insurance",
  D_housing_instability = "Housing stability",
  D_partnership = "Partnership"
)

read1 <- function(prefix, cohort) {
  f <- file.path(OUT, sprintf("%s_%s.csv", prefix, cohort))
  if (!file.exists(f)) stop("Missing required output: ", f)
  read.csv(f, check.names = FALSE)
}

summ <- do.call(rbind, lapply(cohorts, function(c) read1("applied_summary", c)))
diag <- do.call(rbind, lapply(cohorts, function(c) read1("applied_sampler_diagnostics", c)))
arch <- do.call(rbind, lapply(cohorts, function(c) {
  x <- read1("applied_architecture_gates", c)
  x$cohort <- c
  x
}))
summ <- merge(summ, diag[, c("cohort", "max_rhat", "min_ess_bulk", "min_ess_tail",
                             "divergences", "max_treedepth")], by = "cohort", all.x = TRUE)
summ$cohort_label <- factor(summ$cohort, cohorts, c("Pooled", "Black", "White"))
arch$cohort_label <- factor(arch$cohort, cohorts, c("Pooled", "Black", "White"))
arch$domain_label <- factor(unname(domain_map[arch$domain]), rev(unname(domain_map)))
arch$domain_order <- match(arch$domain, names(domain_map))
summ <- summ[order(summ$cohort_label), ]

fmt_num <- function(x, d = 2) formatC(x, digits = d, format = "f")
fmt_ci <- function(m, lo, hi, d = 2) sprintf("%s (%s, %s)", fmt_num(m, d), fmt_num(lo, d), fmt_num(hi, d))
tex_escape <- function(x) {
  x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  gsub("_", "\\\\_", x, fixed = TRUE)
}

write_booktabs <- function(df, file, caption, label, align = NULL, resize = FALSE) {
  if (is.null(align)) align <- paste(rep("l", ncol(df)), collapse = "")
  con <- file(file, "w")
  on.exit(close(con), add = TRUE)
  open_tabular <- sprintf("\\begin{tabular}{%s}", align)
  close_tabular <- "\\end{tabular}"
  if (resize) {
    open_tabular <- c("\\resizebox{\\textwidth}{!}{%", open_tabular)
    close_tabular <- c(close_tabular, "}")
  }
  writeLines(c("\\begin{table}[!htbp]",
               "\\centering",
               sprintf("\\caption{%s}", caption),
               sprintf("\\label{%s}", label),
               open_tabular,
               "\\toprule",
               paste(tex_escape(names(df)), collapse = " & "),
               "\\\\",
               "\\midrule"), con)
  for (i in seq_len(nrow(df))) {
    writeLines(paste(paste(tex_escape(as.character(df[i, ])), collapse = " & "), "\\\\"), con)
  }
  writeLines(c("\\bottomrule", close_tabular, "\\end{table}"), con)
}

summary_tex <- data.frame(
  Cohort = as.character(summ$cohort_label),
  `n` = summ$n,
  `G-computation n` = summ$n_gcomp,
  `Full closure, pp` = fmt_ci(summ$full_closure_pp, summ$full_lo, summ$full_hi, 2),
  `Link ESR` = fmt_ci(summ$ESR_link, summ$ESR_link_lo, summ$ESR_link_hi, 3),
  `P(additive)` = fmt_num(summ$P_additive, 2),
  `P(complementary)` = fmt_num(summ$P_complementary, 2),
  `All favorable` = percent(summ$all_favorable, accuracy = 0.1),
  `max Rhat` = fmt_num(summ$max_rhat, 2),
  `min bulk ESS` = fmt_num(summ$min_ess_bulk, 0),
  `divergences` = summ$divergences,
  check.names = FALSE
)
write_booktabs(summary_tex, file.path(TAB, "fullbayes_applied_summary.tex"),
               "Full-Bayes NHANES applied posterior summaries.",
               "tab:fullbayes-applied-summary", resize = TRUE)
write.csv(summary_tex, file.path(TAB, "fullbayes_applied_summary.csv"), row.names = FALSE)

gate_tex <- arch[order(arch$cohort_label, arch$domain_order),
                 c("cohort_label", "domain_label", "closure_pp", "lo", "hi",
                   "G_link", "G_lo", "G_hi", "P_gate")]
gate_tex <- data.frame(
  Cohort = as.character(gate_tex$cohort_label),
  Domain = as.character(gate_tex$domain_label),
  `Closure contribution, pp` = fmt_ci(gate_tex$closure_pp, gate_tex$lo, gate_tex$hi, 2),
  `Gate index G` = fmt_ci(gate_tex$G_link, gate_tex$G_lo, gate_tex$G_hi, 3),
  `P(G>0)` = fmt_num(gate_tex$P_gate, 2),
  check.names = FALSE
)
write_booktabs(gate_tex, file.path(TAB, "fullbayes_domain_roles.tex"),
               "Domain architecture and anchored gate posterior summaries.",
               "tab:fullbayes-domain-roles")
write.csv(gate_tex, file.path(TAB, "fullbayes_domain_roles.csv"), row.names = FALSE)

theme_csc <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(color = "#4b5563"),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
}

p1 <- ggplot(summ, aes(x = cohort_label, y = full_closure_pp, ymin = full_lo, ymax = full_hi,
                       color = cohort_label)) +
  geom_hline(yintercept = 0, color = "#9ca3af") +
  geom_pointrange(size = 0.9, fatten = 3) +
  scale_color_manual(values = c(Pooled = "#1f77b4", Black = "#7c3aed", White = "#059669"), guide = "none") +
  labs(title = "Full favorable SDOH shift: posterior closure by cohort",
       subtitle = "Weibull PH regularized-horseshoe coalition surface; 10-year absolute-risk closure",
       x = NULL, y = "Gap closure (percentage points)") +
  theme_csc()
ggsave(file.path(FIG, "fig_fullbayes_closure_by_cohort.pdf"), p1, width = 6.5, height = 4.2)
ggsave(file.path(FIG, "fig_fullbayes_closure_by_cohort.png"), p1, width = 6.5, height = 4.2, dpi = 320)

p2 <- ggplot(arch, aes(x = closure_pp, y = domain_label, xmin = lo, xmax = hi, color = cohort_label)) +
  geom_vline(xintercept = 0, color = "#9ca3af") +
  geom_pointrange(position = position_dodge(width = 0.55), fatten = 2, size = 0.6) +
  scale_color_manual(values = c(Pooled = "#1f77b4", Black = "#7c3aed", White = "#059669"), name = NULL) +
  labs(title = "Domain architecture: Shapley-style closure contributions",
       subtitle = "Posterior means and 95% credible intervals on the 10-year risk scale",
       x = "Contribution to closure (percentage points)", y = NULL) +
  theme_csc()
ggsave(file.path(FIG, "fig_fullbayes_domain_architecture.pdf"), p2, width = 8, height = 5.5)
ggsave(file.path(FIG, "fig_fullbayes_domain_architecture.png"), p2, width = 8, height = 5.5, dpi = 320)

p3 <- ggplot(arch, aes(x = G_link, y = domain_label, xmin = G_lo, xmax = G_hi, color = cohort_label)) +
  geom_vline(xintercept = 0, color = "#9ca3af") +
  geom_pointrange(position = position_dodge(width = 0.55), fatten = 2, size = 0.6) +
  scale_color_manual(values = c(Pooled = "#1f77b4", Black = "#7c3aed", White = "#059669"), name = NULL) +
  labs(title = "Anchored gate readout: which domain conditions the payoff of the rest?",
       subtitle = "Positive G means other domains close more gap when the anchor is favorable",
       x = "Gate index G", y = NULL) +
  theme_csc()
ggsave(file.path(FIG, "fig_fullbayes_gate_indices.pdf"), p3, width = 8, height = 5.5)
ggsave(file.path(FIG, "fig_fullbayes_gate_indices.png"), p3, width = 8, height = 5.5, dpi = 320)

cat("Wrote applied tables to", TAB, "and figures to", FIG, "\n")
