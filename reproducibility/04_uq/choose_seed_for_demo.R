# =============================================================================
# choose_seed_for_demo.R
#
# Selects the dataset of the uncertainty quantification study by the rule of
# Appendix H, fixed in advance:
#
#   Generate 10 datasets with the seeds 1 to 10. Run LEWIS (complexity) on
#   each. Take the dataset whose point estimate has the SHD from the true
#   RPDAG closest to the median SHD of LEWIS (complexity) over the 100
#   datasets of the same setting in the comparative study. Ties go to the
#   smaller seed.
#
# The median is read from 02_comparison/results/comparison_cell_summary.csv
# (method lewis_complexity, ER, p = 10, n = 500, rho = 0.5, s = 1, column
# shd_med). It is 2.
#
# Output:
#   results/seed_sweep/seed_sweep_summary.csv    one row per seed (the
#       spread quoted in Appendix H: SHD from the truth 0 to 9, radius 1 to
#       9, truth in the ball in 6 of the 10 datasets)
#   results/seed_sweep/seed_XX/tables/*.csv      the tables of each seed
#
# The run of each seed is exactly the run of uq_decision_rule_demo.R for
# that seed (both call run_uq() of uq_common.R).
#
# Run:   Rscript reproducibility/04_uq/choose_seed_for_demo.R
#        (from the root of the repository, or from any other folder)
# Time:  about two minutes on one core.
# =============================================================================

# ---- Folder of this script and shared code ----------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- file.path(.here, "..", "common")
source(file.path(.here, "uq_common.R"))

# The seeds of the ten datasets of the rule.
SEEDS <- 1:10

# ---- The benchmark median of the setting ------------------------------------
# Read the median SHD of LEWIS (complexity) in the same setting of the
# comparative study.
f_bench <- file.path(.here, "..", "02_comparison", "results", "comparison_cell_summary.csv")
if (!file.exists(f_bench)) stop("File not found: ", f_bench)
bench <- read.csv(f_bench, stringsAsFactors = FALSE)
if (!"shd_med" %in% names(bench)) stop("The column shd_med is missing in ", f_bench)
row <- bench[bench$method == "lewis_complexity" & bench$experiment == "core" &
               bench$graph_type == "ER" & bench$p == UQ_P & bench$n == UQ_N &
               bench$rho == UQ_RHO & bench$effect_scale == UQ_S, ]
if (nrow(row) != 1L) stop("Expected one row of the setting in ", f_bench, ", found ", nrow(row))
BENCH_MEDIAN <- row$shd_med
cat(sprintf("Median SHD of LEWIS (complexity) over the %d datasets of the setting: %g\n",
            row$n_reps, BENCH_MEDIAN))

# ---- Run every seed ------------------------------------------------------------
rows <- list()
for (s in SEEDS) {
  cat(sprintf("\nSeed %d of %d\n", s, length(SEEDS)))
  res <- run_uq(s)
  # Write the tables of this seed and keep its row of the summary.
  write_uq_tables(res, file.path(.here, "results", "seed_sweep", sprintf("seed_%02d", s), "tables"))
  rows[[length(rows) + 1L]] <- res$sweep_row
}
summary_tbl <- do.call(rbind, rows)

# ---- Apply the rule -----------------------------------------------------------
summary_tbl$dist_to_bench_median <- abs(summary_tbl$shd_to_truth - BENCH_MEDIAN)
chosen <- summary_tbl$seed[order(summary_tbl$dist_to_bench_median, summary_tbl$seed)][1]

# ---- Report -----------------------------------------------------------------
cat("\nSeed sweep of the uncertainty quantification study (ER, p = 10, rho = 0.5, n = 500)\n")
print(summary_tbl, row.names = FALSE)
cat(sprintf("\nSHD from the truth: %d to %d. Radius of the ball: %d to %d.\n",
            min(summary_tbl$shd_to_truth), max(summary_tbl$shd_to_truth),
            min(summary_tbl$eps_alpha), max(summary_tbl$eps_alpha)))
cat(sprintf("The ball contains the true RPDAG in %d of the %d datasets.\n",
            sum(summary_tbl$truth_in_ball), nrow(summary_tbl)))
cat(sprintf("Selected seed (SHD closest to %g, ties to the smaller seed): %d\n",
            BENCH_MEDIAN, chosen))

# ---- Write the summary ------------------------------------------------------
out_csv <- file.path(.here, "results", "seed_sweep", "seed_sweep_summary.csv")
write.csv(summary_tbl, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "\n")
