## Build manuscript-ready figures for the Full-Bayes simulation validation.
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

ROOT <- normalizePath("..", winslash = "/", mustWork = TRUE)
SIM <- file.path(ROOT, "simulation_outputs")
FIG <- file.path(ROOT, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

raw_path <- file.path(SIM, "fullbayes_simulation_raw.csv")
rec_path <- file.path(SIM, "fullbayes_recovery_coverage.csv")
cal_path <- file.path(SIM, "fullbayes_regime_calibration.csv")
if (!file.exists(raw_path)) stop("Missing simulation raw file: ", raw_path)
if (!file.exists(rec_path)) stop("Missing recovery file: ", rec_path)
if (!file.exists(cal_path)) stop("Missing calibration file: ", cal_path)

raw <- read.csv(raw_path, stringsAsFactors = FALSE)
rec <- read.csv(rec_path, stringsAsFactors = FALSE)
cal <- read.csv(cal_path, stringsAsFactors = FALSE)

scenario_levels <- c("additive", "complementary", "redundant", "mixed", "gate")
scenario_labels <- c(
  additive = "Additive",
  complementary = "Complementary",
  redundant = "Redundant",
  mixed = "Mixed",
  gate = "Gate"
)

rec$scenario <- factor(rec$scenario, scenario_levels, scenario_labels[scenario_levels])
raw$scenario <- factor(raw$scenario, scenario_levels, scenario_labels[scenario_levels])

metric_map <- c(
  recovery = "Regime recovery",
  esr_coverage = "ESR coverage",
  full_coverage = "Full-closure coverage",
  gate_coverage = "Gate-index coverage",
  gate_domain_recovery = "Gate-domain recovery"
)

long <- do.call(rbind, lapply(names(metric_map), function(m) {
  data.frame(
    scenario = rec$scenario,
    metric = metric_map[[m]],
    value = rec[[m]],
    stringsAsFactors = FALSE
  )
}))
long$metric <- factor(long$metric, unname(metric_map))

pal <- c(
  "Regime recovery" = "#2364AA",
  "ESR coverage" = "#3DA35D",
  "Full-closure coverage" = "#4C9F70",
  "Gate-index coverage" = "#7AC74F",
  "Gate-domain recovery" = "#8B5CF6"
)

p1 <- ggplot(long, aes(x = scenario, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68, color = "white", linewidth = 0.25) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "#3B4556", linewidth = 0.45) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Full-Bayes simulation operating characteristics",
       subtitle = "Primary MCMC validation panel: 25 replications per truth, N = 2,000, five known structures") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#536173"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

save_pair <- function(plot, stem, width = 8.5, height = 5.2) {
  ggsave(file.path(FIG, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(stem, ".png")), plot, width = width, height = height, dpi = 320)
}
save_pair(p1, "fig_fullbayes_sim_recovery_coverage")

tab <- as.data.frame(table(raw$scenario, raw$modal), stringsAsFactors = FALSE)
names(tab) <- c("scenario", "modal", "n")
tab <- subset(tab, modal %in% c("additive", "complementary", "redundant", "mixed"))
tot <- aggregate(n ~ scenario, tab, sum)
names(tot)[2] <- "total"
tab <- merge(tab, tot, by = "scenario")
tab$prop <- tab$n / tab$total
tab$modal <- factor(tab$modal, c("additive", "complementary", "redundant", "mixed"),
                    c("Additive", "Complementary", "Redundant", "Mixed"))
tab$label <- ifelse(tab$n == 0, "", sprintf("%d\n%s", tab$n, percent(tab$prop, accuracy = 1)))

p2 <- ggplot(tab, aes(x = modal, y = scenario, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = label), size = 3.2, lineheight = 0.95, color = "#0F172A") +
  scale_fill_gradient(low = "#EAF2FF", high = "#174EA6", labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Posterior modal regime", y = "Planted truth", fill = "Share",
       title = "Regime recovery matrix",
       subtitle = "Rows are known data-generating structures; cells show modal posterior calls") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#536173"),
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )
save_pair(p2, "fig_fullbayes_sim_regime_confusion", width = 7.4, height = 5.0)

cal$bin_mid <- seq_len(nrow(cal))
p3 <- ggplot(cal, aes(x = mean_confidence, y = accuracy, size = n)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#64748B") +
  geom_point(color = "#2364AA", alpha = 0.9) +
  geom_text(aes(label = n), vjust = -1.0, size = 3.1, color = "#263241") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0.35, 1.02)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  scale_size_continuous(range = c(3, 9), guide = "none") +
  labs(x = "Mean posterior confidence in declared regime", y = "Empirical accuracy",
       title = "Posterior regime-probability calibration",
       subtitle = "Point labels show the number of simulation rows in each confidence bin") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#536173"),
    panel.grid.minor = element_blank()
  )
save_pair(p3, "fig_fullbayes_sim_calibration", width = 6.6, height = 5.1)

cat("Wrote simulation figures to:", normalizePath(FIG, winslash = "\\", mustWork = FALSE), "\n")
