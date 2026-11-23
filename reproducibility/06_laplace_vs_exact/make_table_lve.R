# =============================================================================
# make_table_lve.R
#
# Builds the numbers of Table 1 of the paper (Appendix F) from
# results/exact_vs_laplace.csv, which laplace_vs_exact.R writes.
#
# Output:
#   results/exact_vs_laplace_summary.csv
#       one row per (p, n, prior, chain) with the mean and standard
#       deviation of shd_cpdag over the 100 datasets (shd_cpdag_mean,
#       shd_cpdag_sd), the mean SHD between the CPDAGs of the two chains
#       (agree_mean), and the mean and standard deviation of the computing
#       time in seconds (rt_mean_s, rt_sd_s).
#   results/exact_vs_laplace_runtime_ratio.csv
#       one row per (p, n, prior) with the mean computing times of the two
#       chains (columns Exact and Laplace) and their ratio
#       runtime_ratio_exact_over_laplace = mean time of the exact chain
#       divided by mean time of the Laplace chain, rounded to one decimal.
#       The time ratio of Table 1 is this ratio of MEAN computing times. It
#       is not the mean of the ratios of the datasets.
#   <figure folder>/laplace_vs_exact/table_lve_body.tex
#       the body of Table 1 (the rows between \midrule and \bottomrule).
#
# The script stops if the per-dataset file is incomplete or has duplicated
# or missing values.
#
# Run:   Rscript reproducibility/06_laplace_vs_exact/make_table_lve.R
#        (from the root of the repository, or from any other folder)
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
source(file.path(COMMON_DIR, "paths.R"))

# ---- Read and check the per-dataset results -----------------------------------
res_dir <- file.path(.here, "results")
f <- file.path(res_dir, "exact_vs_laplace.csv")
if (!file.exists(f)) {
  stop("results/exact_vs_laplace.csv is missing. Run laplace_vs_exact.R, ",
       "or place the file from the Zenodo archive (see ../ZENODO.md).")
}
res <- read.csv(f, stringsAsFactors = FALSE)

# 13 pairs (p, n), 100 datasets, 2 priors and 2 chains, without duplicates
# or missing values.
key <- paste(res$p, res$n, res$rep, res$prior, res$method, sep = "|")
stopifnot(nrow(res) == 13L * 100L * 2L * 2L,
          !anyDuplicated(key),
          !anyNA(res[, c("shd_cpdag", "shd_between", "runtime_s")]),
          all(table(paste(res$p, res$n, res$prior, res$method)) == 100L))

# ---- Summary per (p, n, prior, chain) ----------------------------------------------
summary_tbl <- res %>%
  group_by(p, n, n_over_p, prior, method) %>%
  summarize(shd_cpdag_mean = mean(shd_cpdag),
            shd_cpdag_sd = sd(shd_cpdag),
            agree_mean = mean(shd_between),
            rt_mean_s = mean(runtime_s),
            rt_sd_s = sd(runtime_s),
            .groups = "drop") %>%
  arrange(n_over_p, prior, method)
write.csv(summary_tbl, file.path(res_dir, "exact_vs_laplace_summary.csv"),
          row.names = FALSE)

# ---- Ratio of the mean computing times ---------------------------------------------
ratio_tbl <- summary_tbl %>%
  select(p, n, n_over_p, prior, method, rt_mean_s) %>%
  pivot_wider(names_from = method, values_from = rt_mean_s) %>%
  mutate(runtime_ratio_exact_over_laplace = round(Exact / Laplace, 1))
write.csv(ratio_tbl, file.path(res_dir, "exact_vs_laplace_runtime_ratio.csv"),
          row.names = FALSE)

# ---- Table 1 -------------------------------------------------------------------------
# One row per (p, n): mean SHD of L and E and the time ratio, first for the
# complexity prior and then for the uniform prior.
shd_of <- function(pp, nn, pr, m) {
  summary_tbl$shd_cpdag_mean[summary_tbl$p == pp & summary_tbl$n == nn &
                             summary_tbl$prior == pr & summary_tbl$method == m]
}
ratio_of <- function(pp, nn, pr) {
  ratio_tbl$runtime_ratio_exact_over_laplace[ratio_tbl$p == pp & ratio_tbl$n == nn &
                                             ratio_tbl$prior == pr]
}
pn <- unique(summary_tbl[, c("p", "n")])
pn <- pn[order(pn$p, pn$n), ]
lines <- character(0)
for (k in seq_len(nrow(pn))) {
  pp <- pn$p[k]
  nn <- pn$n[k]
  lines <- c(lines, sprintf("%d & %d & %.2f & %.2f & %.1f & %.2f & %.2f & %.1f \\\\",
                            pp, nn,
                            shd_of(pp, nn, "complexity", "Laplace"),
                            shd_of(pp, nn, "complexity", "Exact"),
                            ratio_of(pp, nn, "complexity"),
                            shd_of(pp, nn, "uniform", "Laplace"),
                            shd_of(pp, nn, "uniform", "Exact"),
                            ratio_of(pp, nn, "uniform")))
}
out_tex <- file.path(lewis_fig_dir("laplace_vs_exact"), "table_lve_body.tex")
writeLines(lines, out_tex)
cat("Table 1 (p, n, SHD (L), SHD (E), time ratio, complexity then uniform prior):\n")
cat(lines, sep = "\n")
cat(sprintf("\nWrote %s and the two summary files in results/.\n", out_tex))
