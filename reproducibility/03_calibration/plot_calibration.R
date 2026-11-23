# =============================================================================
# plot_calibration.R
#
# Figures and tables of the calibration study: the remark "Choice of alpha
# and d" (Section 5), Appendix G.2 with Figure 14, and Figures S1 to S4 of
# Section S2 of the online appendix.
#
# Run:   Rscript reproducibility/03_calibration/plot_calibration.R
#        (from the root of the repository, or from any other folder)
#
# Input (written by calibration_sweep.R, or taken from the Zenodo archive):
#     results/d_sweep_results.csv
#     results/alpha_sweep_results.csv
#
# Tables written to results/:
#     d_sweep_cell_summary.csv, alpha_sweep_cell_summary.csv
#         mean, standard deviation, median and quartiles of the SHD over the
#         100 datasets, for each setting and value
#     d_sweep_regret.csv, alpha_sweep_regret.csv
#         for each value: mean, median and largest regret over the 84
#         settings, the number of settings in which the value is best
#         (best_in, ties included), and the number in which its regret is at
#         most 10% (within_10pct)
#     d_sweep_best_d_per_scenario.csv
#         the best value of d in each setting (ties go to the smaller d)
#
# Figures written to <figure folder>/sensitivity/ (PDF and PNG):
#     sensitivity_regret          Figure 14 (Appendix G.2)
#     sensitivity_shd_vs_d        Figure S1 (online appendix)
#     sensitivity_shd_vs_d_pure   Figure S2
#     sensitivity_shd_vs_alpha    Figure S3
#     sensitivity_shd_vs_alpha_pure  Figure S4
#
# REGRET. The regret of a value v in a setting is
#     mean SHD(v) / min over v' of mean SHD(v') - 1,
# where the minimum is over the values of the same sweep. A setting whose
# smallest mean SHD is 0 gives regret 0 to every value.
#
# The script stops if a results file is incomplete or if a figure cannot be
# built.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
})

# ---- Folder of this script and shared code ----------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- file.path(.here, "..", "common")
source(file.path(COMMON_DIR, "plot_theme.R"))

# Folders of the tables and of the figures.
res_dir <- file.path(.here, "results")
fig_dir <- lewis_fig_dir("sensitivity")

# Every curve of these figures is LEWIS, so all curves are drawn with the
# line width of LEWIS in the other figures.
SWEEP_LINEWIDTH <- 1.7

# ---- Read and check one sweep -----------------------------------------------------
read_sweep <- function(name, value_col) {
  f <- file.path(res_dir, sprintf("%s_sweep_results.csv", name))
  if (!file.exists(f)) {
    stop(f, " is missing. Run calibration_sweep.R, or place the file from ",
         "the Zenodo archive (see ../ZENODO.md).")
  }
  x <- read.csv(f, stringsAsFactors = FALSE)
  key <- paste(x$graph_type, x$p, x$n, x$rho, x[[value_col]], x$rep, sep = "|")
  # 84 settings, 10 values and 100 datasets, without duplicates or missing
  # values.
  stopifnot(length(unique(x$scenario)) == 84L,
            length(unique(x[[value_col]])) == 10L,
            nrow(x) == 84L * 10L * 100L,
            !anyDuplicated(key),
            !anyNA(x$shd_rpdag))
  cat(sprintf("%s sweep: %d rows\n", name, nrow(x)))
  x
}

# Mean, standard deviation, median and quartiles of the SHD per setting and
# value.
cell_summary <- function(x, value_col) {
  x %>%
    group_by(scenario, graph_type, p, n, rho, across(all_of(value_col))) %>%
    summarize(n_reps = n(),
              shd_mean = mean(shd_rpdag),
              shd_sd = sd(shd_rpdag),
              shd_med = median(shd_rpdag),
              shd_q25 = quantile(shd_rpdag, 0.25),
              shd_q75 = quantile(shd_rpdag, 0.75),
              .groups = "drop")
}

# Regret of every value in every setting.
add_regret <- function(cell) {
  cell %>%
    group_by(scenario) %>%
    mutate(best = min(shd_mean),
           regret = ifelse(best > 0, shd_mean / best - 1, 0)) %>%
    ungroup()
}

# Summary of the regret of each value over the 84 settings.
regret_table <- function(cell, value_col) {
  add_regret(cell) %>%
    group_by(across(all_of(value_col))) %>%
    summarize(mean_regret = mean(regret),
              median_regret = median(regret),
              max_regret = max(regret),
              best_in = sum(shd_mean <= best + 1e-9),
              within_10pct = sum(regret <= 0.10),
              n_scenarios = n(),
              .groups = "drop")
}

# ---- Tables ---------------------------------------------------------------------------
dd <- read_sweep("d", "d")
aa <- read_sweep("alpha", "alpha")
ds <- cell_summary(dd, "d")
as_ <- cell_summary(aa, "alpha")
write.csv(ds, file.path(res_dir, "d_sweep_cell_summary.csv"), row.names = FALSE)
write.csv(as_, file.path(res_dir, "alpha_sweep_cell_summary.csv"), row.names = FALSE)

# The best d per setting. The rows of ds are ordered by d within a setting,
# so ties go to the smaller d.
best <- ds %>%
  group_by(scenario) %>%
  slice_min(shd_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(scenario, p, n, rho, best_d = d, best_shd_mean = shd_mean)
write.csv(best, file.path(res_dir, "d_sweep_best_d_per_scenario.csv"), row.names = FALSE)

# The regret of each value over the 84 settings.
rt_d <- regret_table(ds, "d")
rt_a <- regret_table(as_, "alpha")
write.csv(rt_d, file.path(res_dir, "d_sweep_regret.csv"), row.names = FALSE)
write.csv(rt_a, file.path(res_dir, "alpha_sweep_regret.csv"), row.names = FALSE)
cat("\nRegret by d (84 settings):\n")
print(as.data.frame(rt_d), row.names = FALSE, digits = 3)
cat("\nRegret by alpha (84 settings):\n")
print(as.data.frame(rt_a), row.names = FALSE, digits = 3)

# ---- Axes of the figures ----------------------------------------------------------------
# The d axis is logarithmic. The value d = 0 (the uniform prior) is drawn at
# a quarter of the smallest positive d and labeled "0". Only the decades are
# labeled.
dpos <- sort(unique(ds$d[ds$d > 0]))
d0x <- dpos[1] / 4
ds$d_x <- ifelse(ds$d == 0, d0x, ds$d)
d_brks <- c(d0x, dpos)
.decades <- c("0.01", "0.1", "1", "10")
d_labs <- ifelse(d_brks == d0x, "0",
                 ifelse(as.character(d_brks) %in% .decades, as.character(d_brks), ""))
D_XLAB <- "Complexity-prior scale d (log scale, leftmost point is d = 0)"
A_XLAB <- "Complexity-prior exponent α"
SHD_YLAB <- "Mean SHD (vs true RPDAG)"

# Every value of alpha has a tick. Only some are labeled, so that labels do
# not collide in narrow panels.
a_brks <- sort(unique(as_$alpha))
a_labs_sweep <- ifelse(round(a_brks, 2) %in% c(0.60, 0.75, 0.90), sprintf("%.2f", a_brks), "")
a_labs_regret <- ifelse(round(a_brks, 2) %in% c(0.55, 0.65, 0.75, 0.85, 0.99),
                        sprintf("%.2f", a_brks), "")

# Label of the sample size of a panel.
n_label <- function(n) factor(paste0("n = ", n), levels = paste0("n = ", sort(unique(n))))

# ---- Mean SHD against the value, rho = 0.5 (Figures S1 and S3) -------------------------
sweep_mix_plot <- function(cell, x_col, vline, xscale) {
  mix <- cell %>% filter(rho == 0.5)
  mix$n_lbl <- n_label(mix$n)
  ggplot(mix, aes(x = .data[[x_col]], y = shd_mean,
                  color = factor(p), group = factor(p))) +
    geom_vline(xintercept = vline, linetype = "dashed", color = "grey40",
               linewidth = 0.5) +
    geom_line(linewidth = SWEEP_LINEWIDTH) +
    geom_point(size = 2.0) +
    xscale +
    scale_color_manual(values = GROUP_COLORS) +
    expand_limits(y = 0) +
    facet_grid(graph_type ~ n_lbl) +
    labs(x = if (x_col == "d_x") D_XLAB else A_XLAB, y = SHD_YLAB,
         color = "Nodes p") +
    theme_lewis() +
    theme(legend.title = element_text(size = rel(0.95)))
}

# ---- Mean SHD against the value, rho in {0, 1} (Figures S2 and S4) ----------------------
sweep_pure_plot <- function(cell, x_col, vline, xscale) {
  pure <- cell %>% filter(rho %in% c(0, 1))
  pure$regime <- factor(ifelse(pure$rho == 0, "All Gaussian", "All non-Gaussian"),
                        levels = c("All Gaussian", "All non-Gaussian"))
  pure$n_lbl <- n_label(pure$n)
  ggplot(pure, aes(x = .data[[x_col]], y = shd_mean,
                   color = factor(p), linetype = graph_type,
                   group = interaction(p, graph_type))) +
    geom_vline(xintercept = vline, linetype = "dashed", color = "grey40",
               linewidth = 0.5) +
    geom_line(linewidth = SWEEP_LINEWIDTH) +
    geom_point(size = 2.0) +
    xscale +
    scale_color_manual(values = GROUP_COLORS) +
    scale_linetype_manual(values = c(ER = "solid", SF = "dotdash")) +
    expand_limits(y = 0) +
    facet_grid(regime ~ n_lbl, scales = "free_y") +
    labs(x = if (x_col == "d_x") D_XLAB else A_XLAB, y = SHD_YLAB,
         color = "Nodes p", linetype = "Graph") +
    theme_lewis() +
    theme(legend.title = element_text(size = rel(0.95)))
}

# Figures S1 to S4.
d_scale <- scale_x_log10(breaks = d_brks, labels = d_labs)
a_scale <- scale_x_continuous(breaks = a_brks, labels = a_labs_sweep)
save_figure(sweep_mix_plot(ds, "d_x", LEWIS_D, d_scale),
            file.path(fig_dir, "sensitivity_shd_vs_d"), width = 13, height = 6.5)
save_figure(sweep_pure_plot(ds, "d_x", LEWIS_D, d_scale),
            file.path(fig_dir, "sensitivity_shd_vs_d_pure"), width = 14, height = 8)
save_figure(sweep_mix_plot(as_, "alpha", LEWIS_ALPHA, a_scale),
            file.path(fig_dir, "sensitivity_shd_vs_alpha"), width = 13, height = 5.0)
save_figure(sweep_pure_plot(as_, "alpha", LEWIS_ALPHA, a_scale),
            file.path(fig_dir, "sensitivity_shd_vs_alpha_pure"), width = 14, height = 8)

# ---- Regret (Figure 14 of Appendix G.2) -------------------------------------------------
# One thin line per setting, colored by p, and the mean regret over the
# settings as a thick black line. Values above 200% are drawn at 200%.
regret_panel <- function(cell, x_col, xlab, vline, xscale) {
  per <- add_regret(cell)
  per$x <- per[[x_col]]
  avg <- per %>% group_by(x) %>% summarize(regret = mean(regret), .groups = "drop")
  ggplot(per, aes(x = x, y = regret)) +
    geom_vline(xintercept = vline, linetype = "dashed", color = "grey40",
               linewidth = 0.5) +
    geom_line(aes(color = factor(p), group = scenario), linewidth = 0.45, alpha = 0.45) +
    geom_line(data = avg, aes(x = x, y = regret), color = "black",
              linewidth = SWEEP_LINEWIDTH, inherit.aes = FALSE) +
    geom_point(data = avg, aes(x = x, y = regret), color = "black",
               size = 2.4, inherit.aes = FALSE) +
    xscale +
    scale_color_manual(values = GROUP_COLORS) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 2), oob = scales::squish) +
    labs(x = xlab, y = "Regret", color = "Nodes p") +
    theme_lewis() +
    theme(legend.title = element_text(size = rel(0.95)))
}
gd <- regret_panel(ds, "d_x", D_XLAB, LEWIS_D, d_scale)
ga <- regret_panel(as_, "alpha", A_XLAB, LEWIS_ALPHA,
                   scale_x_continuous(breaks = a_brks, labels = a_labs_regret))
p_reg <- (gd | ga) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_figure(p_reg, file.path(fig_dir, "sensitivity_regret"), width = 13, height = 5.6)

cat("\nAll five figures of the calibration study are written.\n")
