# =============================================================================
# plot_motivation.R
#
# Figure 1 of the paper, and the two summary tables of the study.
#
# Reads results/motivation_results.csv (the per-dataset results, from the
# Zenodo archive or from motivation_experiment.R and
# pairwise_lingam/import_motivation_pairwise.R) and writes
#
#     <figure folder>/motivation/motivation_intro_2x3.pdf and .png
#         Figure 1. Mean SHD over 100 datasets from the true RPDAG against
#         the number k of nodes with t5 errors (k decreases from left to
#         right). Columns: n = 100, 250, 500. Top row: the other nodes have
#         Gaussian errors. Bottom row: the other nodes have t15 errors.
#         Of the two versions of LEWIS, the figure shows only LEWIS
#         (complexity), with the legend entry "LEWIS". It keeps the color,
#         line type and marker of LEWIS (complexity) in the other figures.
#     results/motivation_cell_summary.csv
#         mean, standard deviation, median and quartiles of the SHD for
#         each regime, n and method
#     results/motivation_gaussianity_table.csv
#         median Kolmogorov-Smirnov p-value and mean excess kurtosis of the
#         errors for each regime and n
#
# All six panels share one vertical axis that starts at zero. The upper
# limit is computed from all methods of both rows, including LEWIS
# (uniform), as in the figure of the paper.
#
# The script stops if a check fails or if the figure cannot be built.
#
# Run from any folder, for example from the root of the repository:
#     Rscript reproducibility/01_motivation/plot_motivation.R
# The figure folder is set by the environment variable LEWIS_FIG_DIR (see
# common/paths.R).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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
suppressPackageStartupMessages({
  source(file.path(COMMON_DIR, "scenario_instances.R"))
  source(file.path(COMMON_DIR, "paper12_instance.R"))
  source(file.path(.here, "motivation_config.R"))
})

# Folders of the results and of the figure. lewis_fig_dir() creates the
# folder of the figure.
res_dir <- file.path(.here, "results")
fig_dir <- lewis_fig_dir("motivation")

# ---- Read the per-dataset results --------------------------------------------
f_raw <- file.path(res_dir, "motivation_results.csv")
if (!file.exists(f_raw)) {
  stop("The file ", f_raw, " is missing. Get it from the Zenodo archive ",
       "(see reproducibility/ZENODO.md) or run motivation_experiment.R and ",
       "pairwise_lingam/import_motivation_pairwise.R.")
}
d <- read.csv(f_raw, stringsAsFactors = FALSE)
# Keep the 25 regimes of Figure 1.
d <- d[d$regime %in% REGIME_NAMES, ]
cat(sprintf("%d rows, methods: %s\n", nrow(d),
            paste(sort(unique(d$method)), collapse = ", ")))

# ---- Checks: 100 distinct datasets per setting and method, nothing missing ---
key <- paste(d$regime, d$n, d$method, d$rep, sep = "|")
stopifnot("duplicated rows" = !anyDuplicated(key))
stopifnot("failed rows" = !any(as.logical(d$skipped)))
stopifnot("missing SHD values" = !anyNA(d$shd_rpdag))
stopifnot("not all ten methods are present" = setequal(unique(d$method), METHOD_ORDER))
reps <- table(paste(d$regime, d$n, d$method, sep = "|"))
stopifnot("a setting and method without 100 datasets" =
            length(reps) == length(REGIME_NAMES) * length(N_VALUES) * length(METHOD_ORDER),
          all(reps == N_REPS))

# Labels of the methods and of the sample sizes.
d$Method <- method_factor(d$method)
d$n_lbl <- factor(paste0("n = ", d$n), levels = paste0("n = ", sort(unique(d$n))))

# ---- Table of the properties of the errors -----------------------------------
gauss_tab <- d %>%
  group_by(regime, family, k_ng, df, n) %>%
  summarize(ks_p_median = median(ks_pvalue_median),
            kurtosis_mean = mean(excess_kurtosis),
            .groups = "drop") %>%
  arrange(n, k_ng, df)
write.csv(gauss_tab, file.path(res_dir, "motivation_gaussianity_table.csv"),
          row.names = FALSE)
cat("Wrote results/motivation_gaussianity_table.csv\n")

# ---- Summary of the SHD for each regime, n and method ------------------------
summ <- d %>%
  group_by(regime, family, k_ng, df, n, n_lbl, method, Method) %>%
  summarize(n_reps = n(),
            res_mean = mean(shd_rpdag),
            res_sd = sd(shd_rpdag),
            res_med = median(shd_rpdag),
            res_q25 = quantile(shd_rpdag, 0.25),
            res_q75 = quantile(shd_rpdag, 0.75),
            cp_med = median(shd_vs_cpdag),
            ks_p = median(ks_pvalue_median),
            .groups = "drop")
write.csv(summ %>% select(-Method, -n_lbl),
          file.path(res_dir, "motivation_cell_summary.csv"), row.names = FALSE)
cat(sprintf("Wrote results/motivation_cell_summary.csv (%d rows)\n", nrow(summ)))

# ---- Data of the two rows ---------------------------------------------------
# The k of a row comes from INTRO_ROWS. The regime all_t5 is k = 12 in both
# rows.
intro_row_data <- function(key) {
  summ %>%
    inner_join(INTRO_ROWS[[key]], by = "regime") %>%
    mutate(k_ng = k_ng.y) %>%
    select(-k_ng.x, -k_ng.y)
}
r1 <- intro_row_data("t5_gauss")
r2 <- intro_row_data("t5_t15")
stopifnot(nrow(distinct(r1, regime, n)) == 39L, nrow(distinct(r2, regime, n)) == 39L)

# Common vertical axis: from zero to the largest mean SHD of both rows plus
# 6 percent.
ylim_for <- function(...) {
  m <- max(unlist(lapply(list(...), function(z) z$res_mean)))
  c(0, ceiling(m * 1.06))
}

# ---- Legend of Figure 1: LEWIS (complexity) is shown as "LEWIS" -------------
.intro_rename <- function(x) {
  x[x == "LEWIS (complexity)"] <- "LEWIS"
  x
}
# Drop LEWIS (uniform) and rename LEWIS (complexity).
intro_only_lewis <- function(dat) {
  dat <- dat[as.character(dat$Method) != "LEWIS (uniform)", , drop = FALSE]
  lv <- .intro_rename(setdiff(levels(dat$Method), "LEWIS (uniform)"))
  dat$Method <- factor(.intro_rename(as.character(dat$Method)), levels = lv)
  dat
}
# The method scales of common/plot_theme.R without LEWIS (uniform).
.intro_keys <- setdiff(names(METHOD_COLORS), "LEWIS (uniform)")
.intro_vals <- function(v) {
  v <- v[.intro_keys]
  names(v) <- .intro_rename(names(v))
  v
}
scale_method_intro <- function() {
  g <- guide_legend(nrow = 2)
  list(scale_color_manual(values = .intro_vals(METHOD_COLORS)),
       scale_shape_manual(values = .intro_vals(METHOD_SHAPES)),
       scale_linetype_manual(values = .intro_vals(METHOD_LINETYPES)),
       scale_linewidth_manual(values = .intro_vals(METHOD_LINEWIDTHS)),
       scale_size_manual(values = .intro_vals(METHOD_POINTSIZES)),
       guides(color = g, shape = g, linetype = g, linewidth = g, size = g))
}

# Row labels with the degrees of freedom as subscripts.
INTRO_ROW_EXPR <- list(
  t5_gauss = expression(paste("Non-Gaussian nodes ", t[5], ", remaining nodes exactly Gaussian")),
  t5_t15 = expression(paste("Non-Gaussian nodes ", t[5], ", remaining nodes ", t[15])))

# One row of three panels. The horizontal axis is reversed, so that k
# decreases from left to right. The figure is drawn at 16 x 11 inches and
# scaled down in the paper, so the base font size is 18.
panel_intro <- function(dat, ylim, base_size = 18) {
  dat <- intro_only_lewis(dat)
  ggplot(dat, aes(x = k_ng, y = res_mean, color = Method, shape = Method,
                  linetype = Method, linewidth = Method, size = Method,
                  group = Method)) +
    geom_line() +
    geom_point() +
    scale_x_reverse(breaks = seq(12, 0, -2)) +
    coord_cartesian(ylim = ylim) +
    scale_method_intro() +
    facet_wrap(~ n_lbl, nrow = 1) +
    labs(x = expression(paste("Number of Student-", t[5], " nodes ", italic(k))),
         y = "Mean SHD (vs true RPDAG)") +
    theme_lewis(base_size = base_size)
}

# ---- Figure 1 -----------------------------------------------------------------
# The two rows share one vertical axis and one legend at the bottom.
yl <- ylim_for(r1, r2)
top <- panel_intro(r1, yl) + labs(subtitle = INTRO_ROW_EXPR[["t5_gauss"]])
bot <- panel_intro(r2, yl) + labs(subtitle = INTRO_ROW_EXPR[["t5_t15"]])
fig1 <- (top / bot) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_figure(fig1, file.path(fig_dir, "motivation_intro_2x3"), width = 16, height = 11)
stopifnot(file.exists(file.path(fig_dir, "motivation_intro_2x3.pdf")))

# ---- Tables of the mean SHD shown in Figure 1 ---------------------------------
# One printed table per row and n, with one column per method.
for (key in names(INTRO_ROWS)) {
  d23 <- intro_row_data(key)
  for (nn in sort(unique(d23$n))) {
    cat(sprintf("\n%s, n = %d: mean SHD from the true RPDAG\n",
                INTRO_ROW_LABELS[[key]], nn))
    w <- d23 %>% filter(n == nn) %>% select(k_ng, Method, res_mean) %>%
      arrange(desc(k_ng)) %>%
      pivot_wider(names_from = Method, values_from = res_mean)
    print(as.data.frame(w), row.names = FALSE, digits = 3)
  }
}
