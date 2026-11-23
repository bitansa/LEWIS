# =============================================================================
# aggregate_and_plot.R
#
# Final step of the comparative study (Section 7.2, Appendix G.3 and
# Section S1 of the online appendix). It reads the ten files
# results/results_<method>.csv and writes
#
#     results/comparison_all.csv           all rows of all methods
#     results/comparison_cell_summary.csv  one row per setting and method:
#                                          mean, sd, median and quartiles of
#                                          the SHD, mean, sd and median of the
#                                          two SID bounds, mean run time
#
# and the figures and tables of the paper, in lewis_fig_dir("comparison"):
#
#     comparison_shd_vs_p.pdf       Figure 9 (Section 7.2)
#     comparison_shd_vs_rho.pdf     Figure 10 (Section 7.2)
#     comparison_shd_vs_n.pdf       Figure 15 (Appendix G.3)
#     comparison_shd_vs_scale.pdf   Figure 16 (Appendix G.3)
#     tables/shd_core_ER.tex, tables/shd_core_SF.tex            Tables S1, S2
#     tables/shd_mixture_ER.tex, tables/shd_mixture_SF.tex      Tables S3, S4
#     tables/shd_effect_ER.tex, tables/shd_effect_SF.tex        Tables S5, S6
#     tables/sid_lower_core_ER.tex, tables/sid_upper_core_ER.tex Tables S7, S8
#     tables/sid_lower_core_SF.tex, tables/sid_upper_core_SF.tex Tables S9, S10
#     tables/sid_mixture_ER.tex, tables/sid_mixture_SF.tex      Tables S11, S12
#     tables/sid_effect_ER.tex, tables/sid_effect_SF.tex        Tables S13, S14
#
# Each figure is written as PDF and PNG. The tables are LaTeX tabular
# environments that the online appendix inputs.
#
# CHECKS. The script stops unless the results are complete: the ten
# methods, the 54 settings, exactly 100 distinct datasets for each setting
# and method, no duplicated row, no missing value and no failed row.
#
# Each setting is computed once but may appear in several figures. The core
# settings (rho = 0.5, s = 1) give the figures against n and p, and also the
# point rho = 0.5 of the figure against rho and the point s = 1 of the
# figure against s.
#
# Run from any folder, for example from the root of the repository:
#     Rscript reproducibility/02_comparison/aggregate_and_plot.R
# The figure folder is set by the environment variable LEWIS_FIG_DIR (see
# common/paths.R).
#
# Dependencies: dplyr, tidyr, ggplot2.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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
source(file.path(.here, "comparison_cells.R"))

# Folders of the results, the figures and the tables. lewis_fig_dir()
# creates the folders of the figures and tables.
res_dir <- file.path(.here, "results")
fig_dir <- lewis_fig_dir("comparison")
tab_dir <- lewis_fig_dir("comparison", "tables")

# ---- Read the ten result files --------------------------------------------------
files <- file.path(res_dir, paste0("results_", METHOD_ORDER, ".csv"))
missing <- files[!file.exists(files)]
if (length(missing) > 0L) {
  stop("These result files are missing: ", paste(basename(missing), collapse = ", "),
       ". Run the jobs of this study or place the raw results from the archive",
       " (see ZENODO.md) in results/.")
}
all_rows <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))
cat(sprintf("Read %d files with %d rows.\n", length(files), nrow(all_rows)))

# ---- Checks of completeness -------------------------------------------------------
cells <- comparison_cells()
all_rows$cell <- cell_key(all_rows$experiment, all_rows$graph_type, all_rows$p,
                          all_rows$n, all_rows$rho, all_rows$effect_scale)
grid_keys <- cell_key(cells$experiment, cells$graph_type, cells$p, cells$n,
                      cells$rho, cells$effect_scale)
# The expected columns are present.
stopifnot(all(c("experiment", "graph_type", "p", "n", "rho", "effect_scale", "rep",
                "method", "shd_rpdag", "sid_lower", "sid_upper", "n_edges_est",
                "runtime_s", "skipped") %in% names(all_rows)))
# The ten methods and the 54 settings, and nothing else.
stopifnot(setequal(unique(all_rows$method), METHOD_ORDER))
stopifnot(setequal(unique(all_rows$cell), grid_keys))
# No duplicated row.
stopifnot(!anyDuplicated(all_rows[, c("cell", "method", "rep")]))
# Exactly the datasets 1 to 100 for each setting and method.
reps <- all_rows %>% group_by(cell, method) %>%
  summarize(n_rep = n_distinct(rep), ok = setequal(rep, seq_len(N_REPS)),
            .groups = "drop")
stopifnot(nrow(reps) == length(grid_keys) * length(METHOD_ORDER),
          all(reps$n_rep == N_REPS), all(reps$ok))
# No failed row and no missing value.
stopifnot(!any(as.logical(all_rows$skipped)))
stopifnot(!anyNA(all_rows[, c("shd_rpdag", "sid_lower", "sid_upper", "n_edges_est")]))
cat("Checks passed: 10 methods, 54 settings, 100 datasets each, no missing value.\n")

# Write all rows of all methods in one file.
all_rows$cell <- NULL
write.csv(all_rows, file.path(res_dir, "comparison_all.csv"), row.names = FALSE)

# ---- Summary of each setting and method -----------------------------------------
# The checks above ensure 100 datasets and no missing value, so no na.rm.
summ <- all_rows %>%
  group_by(experiment, graph_type, p, n, rho, effect_scale, method) %>%
  summarize(n_reps = n(),
            n_ok = sum(!skipped),
            shd_mean = mean(shd_rpdag),
            shd_sd = sd(shd_rpdag),
            shd_med = median(shd_rpdag),
            shd_q25 = quantile(shd_rpdag, 0.25),
            shd_q75 = quantile(shd_rpdag, 0.75),
            sid_lo_mean = mean(sid_lower),
            sid_lo_sd = sd(sid_lower),
            sid_up_mean = mean(sid_upper),
            sid_up_sd = sd(sid_upper),
            sid_lo_med = median(sid_lower),
            sid_up_med = median(sid_upper),
            runtime_mean_s = mean(runtime_s),
            .groups = "drop")
write.csv(summ, file.path(res_dir, "comparison_cell_summary.csv"), row.names = FALSE)
cat(sprintf("Wrote comparison_cell_summary.csv (%d rows).\n", nrow(summ)))

# Labels of the methods and of the panels, shared by the four figures.
summ$Method <- method_factor(summ$method)
lab_p <- labeller(p = function(v) paste0("p = ", v))
lab_n <- labeller(n = function(v) paste0("n = ", v))
Y_LAB <- "Mean SHD (vs true RPDAG)"

# The aesthetics shared by every figure.
method_aes <- aes(color = Method, shape = Method, linetype = Method,
                  linewidth = Method, size = Method, group = Method)

# ---- Figure 15: core settings, mean SHD against n --------------------------
# The rows are the values of p and the columns are the graph types. The two
# panels of one p share one y-axis from zero.
core <- summ %>% filter(experiment == "core")
p_n <- ggplot(core, aes(x = n, y = shd_mean)) + method_aes +
  geom_line() + geom_point() +
  scale_x_log10(breaks = sort(unique(core$n))) +
  expand_limits(y = 0) +
  scale_method_all() +
  facet_grid(p ~ graph_type, scales = "free_y", labeller = lab_p) +
  labs(x = "Sample size n (log scale)", y = Y_LAB) +
  theme_lewis()
save_figure(p_n, file.path(fig_dir, "comparison_shd_vs_n"), width = 10, height = 13)

# ---- Figure 9: core settings, mean SHD against p ----------------------------
# The same settings at n = 500 and n = 1000, the sample sizes present for
# every p. All panels share one y-axis.
core_p <- core %>% filter(n %in% c(500L, 1000L))
p_p <- ggplot(core_p, aes(x = p, y = shd_mean)) + method_aes +
  geom_line() + geom_point() +
  scale_x_log10(breaks = sort(unique(core_p$p))) +
  expand_limits(y = 0) +
  scale_method_all() +
  facet_grid(graph_type ~ n, labeller = lab_n) +
  labs(x = "Number of nodes p (log scale)", y = Y_LAB) +
  theme_lewis()
save_figure(p_p, file.path(fig_dir, "comparison_shd_vs_p"), width = 13, height = 6.5)

# ---- Figure 10: mean SHD against rho -----------------------------------------
mix <- summ %>%
  filter(experiment %in% c("mixture", "core"), n == 500L,
         p %in% c(20L, 50L), effect_scale == 1.0)
p_rho <- ggplot(mix, aes(x = rho, y = shd_mean)) + method_aes +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  expand_limits(y = 0) +
  scale_method_all() +
  facet_grid(p ~ graph_type, scales = "free_y", labeller = lab_p) +
  labs(x = "Fraction of non-Gaussian nodes \u03c1", y = Y_LAB) +
  theme_lewis()
save_figure(p_rho, file.path(fig_dir, "comparison_shd_vs_rho"), width = 10, height = 8)

# ---- Figure 16: mean SHD against s ---------------------------------------------
esc <- summ %>%
  filter(experiment %in% c("effect_scale", "core"), n == 500L,
         p %in% c(20L, 50L), rho == 0.5)
p_s <- ggplot(esc, aes(x = effect_scale, y = shd_mean)) + method_aes +
  geom_line() + geom_point() +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2)) +
  expand_limits(y = 0) +
  scale_method_all() +
  facet_grid(p ~ graph_type, scales = "free_y", labeller = lab_p) +
  labs(x = "Effect-size scale s (log scale)", y = Y_LAB) +
  theme_lewis()
save_figure(p_s, file.path(fig_dir, "comparison_shd_vs_scale"), width = 10, height = 8)

# =============================================================================
# LaTeX tables: mean (sd) of the SHD and of the two SID bounds.
# =============================================================================
# Each table is a tabular environment made of stacked blocks, one block for
# each p (or each bound and p). A block has a title row, its own column
# header and one row per method. Every entry is "mean (sd)" over the 100
# datasets, with one decimal. In each column of a block, the smallest mean
# at the printed precision is in bold. The methods appear in METHOD_ORDER.

# Remove old tables so that the folder holds only the current ones.
unlink(list.files(tab_dir, pattern = "\\.tex$", full.names = TRUE))

# Format "mean (sd)" with one decimal.
fmt_cell <- function(m, s) sprintf("%.1f (%.1f)", m, s)
# Escape underscores for LaTeX.
tex_escape <- function(x) gsub("_", "\\\\_", x)
# Round to the printed precision.
round_print <- function(v) as.numeric(sprintf("%.1f", v))

# The rows of one block: the method names and a character matrix of
# formatted entries (methods x columns), with the best entry in bold.
block_rows <- function(dat, col_levels) {
  dat <- dat %>% mutate(col_label = factor(col_label, levels = col_levels))
  stopifnot(!anyDuplicated(dat[, c("Method", "col_label")]), !anyNA(dat$mean))
  wide_m <- dat %>% select(Method, col_label, mean) %>% arrange(Method) %>%
    pivot_wider(names_from = col_label, values_from = mean)
  wide_s <- dat %>% select(Method, col_label, sd) %>% arrange(Method) %>%
    pivot_wider(names_from = col_label, values_from = sd)
  stopifnot(identical(wide_m$Method, wide_s$Method))
  cols <- col_levels[col_levels %in% names(wide_m)]
  stopifnot(length(cols) == length(col_levels))
  # The smallest printed mean of each column.
  best <- vapply(cols, function(cc) min(round_print(wide_m[[cc]])), numeric(1))
  cells <- matrix("", nrow(wide_m), length(cols))
  for (i in seq_len(nrow(wide_m))) for (j in seq_along(cols)) {
    m_ij <- wide_m[[cols[j]]][i]
    txt <- fmt_cell(m_ij, wide_s[[cols[j]]][i])
    cells[i, j] <- if (round_print(m_ij) == best[j]) paste0("\\textbf{", txt, "}") else txt
  }
  list(methods = tex_escape(as.character(wide_m$Method)), cols = cols, cells = cells)
}

# Write one table of stacked blocks. The table is as wide as the widest
# block. Narrower blocks leave their last columns empty.
write_block_table <- function(blocks, file, metric_name) {
  rows <- lapply(blocks, function(bl) block_rows(bl$dat, bl$col_levels))
  K <- max(vapply(rows, function(b) length(b$cols), integer(1)))
  pad <- function(x) c(x, rep("", K - length(x)))
  lines <- c(
    sprintf("%% %s: mean (sd) over 100 datasets. The best mean in each column of a block is in bold.",
            metric_name),
    sprintf("\\begin{tabular}{l%s}", strrep("r", K)),
    "\\toprule")
  for (b in seq_along(rows)) {
    r <- rows[[b]]
    k <- length(r$cols)
    if (b > 1L) lines <- c(lines, "\\addlinespace[6pt]")
    lines <- c(lines,
      sprintf("& \\multicolumn{%d}{c}{%s}%s \\\\", k, blocks[[b]]$title,
              strrep(" &", K - k)),
      sprintf("\\cmidrule(lr){2-%d}", k + 1L),
      paste0("Method & ", paste(pad(r$cols), collapse = " & "), " \\\\"),
      "\\midrule",
      paste0(r$methods, " & ",
             apply(r$cells, 1, function(z) paste(pad(z), collapse = " & ")),
             " \\\\"))
    if (b < length(rows)) lines <- c(lines, "\\midrule")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file)
  cat(sprintf("Wrote %s (%d blocks, up to %d columns).\n", basename(file),
              length(rows), K))
}

# The mean and sd of one metric for every setting and method.
metric_long <- function(d, mean_col, sd_col) {
  d %>% transmute(Method, experiment, graph_type, p, n, rho, effect_scale,
                  mean = .data[[mean_col]], sd = .data[[sd_col]])
}
ML <- list(shd = metric_long(summ, "shd_mean", "shd_sd"),
           sid_lower = metric_long(summ, "sid_lo_mean", "sid_lo_sd"),
           sid_upper = metric_long(summ, "sid_up_mean", "sid_up_sd"))
METRIC_NAMES <- c(shd = "SHD from the true RPDAG",
                  sid_lower = "Lower SID bound against the true DAG",
                  sid_upper = "Upper SID bound against the true DAG")
BOUND_TITLE <- c(sid_lower = "Lower bound", sid_upper = "Upper bound")

# Block of the core settings of one graph type and one p. The columns are n.
blk_core <- function(ml, gt, pp, title) {
  d <- ml %>% filter(experiment == "core", graph_type == gt, p == pp) %>%
    mutate(col_label = sprintf("$n = %d$", n))
  list(title = title, dat = d, col_levels = sprintf("$n = %d$", sort(unique(d$n))))
}
# Block of the settings that vary rho. The columns are rho. The column
# rho = 0.5 is a core setting.
blk_mix <- function(ml, gt, pp, title) {
  d <- ml %>% filter(experiment %in% c("mixture", "core"), n == 500L,
                     graph_type == gt, p == pp, effect_scale == 1.0) %>%
    mutate(col_label = sprintf("$\\rho = %g$", rho))
  list(title = title, dat = d,
       col_levels = sprintf("$\\rho = %g$", sort(unique(d$rho))))
}
# Block of the settings that vary s. The columns are s. The column s = 1 is
# a core setting.
blk_eff <- function(ml, gt, pp, title) {
  d <- ml %>% filter(experiment %in% c("effect_scale", "core"), n == 500L,
                     graph_type == gt, p == pp, rho == 0.5) %>%
    mutate(col_label = sprintf("$s = %g$", effect_scale))
  list(title = title, dat = d,
       col_levels = sprintf("$s = %g$", sort(unique(d$effect_scale))))
}
P_CORE <- c(10L, 20L, 50L, 100L)
P_FAM <- c(20L, 50L)

# ---- SHD tables (Tables S1 to S6) -------------------------------------------------
for (gt in c("ER", "SF")) {
  write_block_table(lapply(P_CORE, function(pp) blk_core(ML$shd, gt, pp, sprintf("$p = %d$", pp))),
                    file.path(tab_dir, sprintf("shd_core_%s.tex", gt)),
                    sprintf("%s, rho = 0.5, s = 1, %s graphs", METRIC_NAMES[["shd"]], gt))
  write_block_table(lapply(P_FAM, function(pp) blk_mix(ML$shd, gt, pp, sprintf("$p = %d$", pp))),
                    file.path(tab_dir, sprintf("shd_mixture_%s.tex", gt)),
                    sprintf("%s, n = 500, s = 1, each rho, %s graphs", METRIC_NAMES[["shd"]], gt))
  write_block_table(lapply(P_FAM, function(pp) blk_eff(ML$shd, gt, pp, sprintf("$p = %d$", pp))),
                    file.path(tab_dir, sprintf("shd_effect_%s.tex", gt)),
                    sprintf("%s, n = 500, rho = 0.5, each s, %s graphs", METRIC_NAMES[["shd"]], gt))
}

# ---- SID tables (Tables S7 to S14) ---------------------------------------------------
for (gt in c("ER", "SF")) {
  # Core settings: one table for each bound, with one block for each p.
  for (mk in c("sid_lower", "sid_upper")) {
    write_block_table(lapply(P_CORE, function(pp) blk_core(ML[[mk]], gt, pp, sprintf("$p = %d$", pp))),
                      file.path(tab_dir, sprintf("%s_core_%s.tex", mk, gt)),
                      sprintf("%s, rho = 0.5, s = 1, %s graphs", METRIC_NAMES[[mk]], gt))
  }
  # Settings that vary rho or s: both bounds in one table, four blocks.
  blocks_for <- function(builder) {
    unlist(lapply(c("sid_lower", "sid_upper"), function(mk) lapply(P_FAM, function(pp)
      builder(ML[[mk]], gt, pp, sprintf("%s, $p = %d$", BOUND_TITLE[[mk]], pp)))),
      recursive = FALSE)
  }
  write_block_table(blocks_for(blk_mix), file.path(tab_dir, sprintf("sid_mixture_%s.tex", gt)),
                    sprintf("SID bounds against the true DAG, n = 500, s = 1, each rho, %s graphs", gt))
  write_block_table(blocks_for(blk_eff), file.path(tab_dir, sprintf("sid_effect_%s.tex", gt)),
                    sprintf("SID bounds against the true DAG, n = 500, rho = 0.5, each s, %s graphs", gt))
}

# ---- Check that all 14 table files exist ----------------------------------------------
expected_tex <- c(sprintf("shd_%s_%s.tex", rep(c("core", "mixture", "effect"), each = 2), c("ER", "SF")),
                  sprintf("%s_core_%s.tex", rep(c("sid_lower", "sid_upper"), 2), rep(c("ER", "SF"), each = 2)),
                  sprintf("sid_%s_%s.tex", rep(c("mixture", "effect"), each = 2), c("ER", "SF")))
stopifnot(all(file.exists(file.path(tab_dir, expected_tex))))
cat(sprintf("\nFigures in %s\nTables in %s\n", fig_dir, tab_dir))
