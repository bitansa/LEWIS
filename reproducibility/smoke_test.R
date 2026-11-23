# =============================================================================
# smoke_test.R
#
# A quick check, in a few minutes on a laptop, that the package LEWIS and
# the shared code of this folder are installed and give the results of the
# paper. Run it first on any new machine:
#
#     Rscript reproducibility/smoke_test.R
#
# It checks six things.
#   1. The seeds and the fixed instances are those of the paper.
#   2. The datasets are reproducible.
#   3. lewis() with center = FALSE reproduces stored results of the
#      comparative study (Section 7.2) on three datasets of the setting
#      ER, p = 10, n = 50, rho = 0.5, s = 1, and of the calibration study
#      (Appendix G.2) on one dataset.
#   4. The RPDAG, the SHD and the CPDAG of pcalg agree.
#   5. The evaluation of the comparative study (SHD and SID) reproduces two
#      stored rows, if the packages pcalg and SID are installed.
#   6. The run of the uncertainty quantification study (Section 7.3) on the
#      selected dataset gives the numbers of the paper.
#
# The stored values are in the per-dataset files of the Zenodo archive (see
# ZENODO.md). They are written in this script, so the archive is not
# needed to run it.
#
# The script ends with "All checks passed." and exit status 0, or lists the
# failed checks and exits with status 1.
# =============================================================================

# ---- Folder of this script and shared code ----------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- file.path(.here, "common")
t_start <- Sys.time()

# Load the package and the shared code.
suppressPackageStartupMessages(library(LEWIS))
source(file.path(COMMON_DIR, "settings.R"))
source(file.path(COMMON_DIR, "scenario_instances.R"))
source(file.path(COMMON_DIR, "graph_tools.R"))
cat(sprintf("LEWIS %s, R %s\n\n", utils::packageVersion("LEWIS"), getRversion()))

# ---- Helper: record one check ------------------------------------------------
failed <- character(0)
check <- function(label, ok) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAILED", label))
  if (!ok) failed <<- c(failed, label)
  invisible(ok)
}

# ---- 1. Seeds and fixed instances --------------------------------------------
cat("1. Seeds and fixed instances\n")
check("stable_seed(\"instance|ER|p10\") is 803763713",
      stable_seed("instance|ER|p10") == 803763713L)
inst <- build_instance("ER", 10L)
view <- instance_view(inst, rho = 0.5, effect_scale = 1)
check("the ER graph with p = 10 has 9 edges (Section 7.3)", sum(inst$dag) == 9L)
check("its non-Gaussian nodes at rho = 0.5 are 1, 2, 3, 7, 10 (Section 7.3)",
      identical(as.integer(view$nG_star), c(1L, 2L, 3L, 7L, 10L)))
check("build_instance() leaves the random numbers of the caller unchanged", {
  set.seed(1)
  a <- runif(1)
  set.seed(1)
  .instance_memo[["instance|SF|p20"]] <- NULL
  build_instance("SF", 20L)
  identical(a, runif(1))
})

# ---- 2. Reproducible datasets -------------------------------------------------
cat("2. Datasets\n")
X1 <- draw_dataset(view, n = 50L, rep = 1L, regime_tag = "mixture")
X1b <- draw_dataset(view, n = 50L, rep = 1L, regime_tag = "mixture")
check("a dataset is the same when it is drawn twice", identical(X1, X1b))
check("the dataset has 50 rows and 10 columns", identical(dim(X1), c(50L, 10L)))

# ---- 3. lewis() reproduces stored results ------------------------------------
# The stored values are the SHD from the true RPDAG in
# 02_comparison/results/results_lewis_complexity.csv and
# results_lewis_uniform.csv (setting core, ER, p = 10, n = 50, datasets 1 to
# 3), and in 03_calibration/results/d_sweep_results.csv (ER, p = 10, n = 50,
# rho = 0, d = 0.03, dataset 1). As in the studies, lewis() is called right
# after draw_dataset() and without a seed.
cat("3. lewis() on stored datasets\n")
truth <- unname(rpdag(view$dag, view$nG_star))
expected <- list(lewis_complexity = c(7L, 9L, 9L), lewis_uniform = c(7L, 9L, 10L))
for (m in names(expected)) {
  got <- integer(0)
  for (r in 1:3) {
    X <- draw_dataset(view, n = 50L, rep = r, regime_tag = "mixture")
    fit <- if (m == "lewis_complexity") {
      lewis(X, d = LEWIS_D, alpha = LEWIS_ALPHA, center = FALSE, ess_min = 0)
    } else {
      lewis(X, d = 0, center = FALSE, ess_min = 0)
    }
    got <- c(got, LEWIS::shd(unname(fit$D_hat), truth))
  }
  check(sprintf("%s: SHD %s on datasets 1 to 3 (stored %s)", METHOD_LABELS[[m]],
                paste(got, collapse = ", "), paste(expected[[m]], collapse = ", ")),
        identical(got, expected[[m]]))
}
view0 <- instance_view(inst, rho = 0, effect_scale = 1)
X <- draw_dataset(view0, n = 50L, rep = 1L, regime_tag = "mixture")
fit <- lewis(X, d = LEWIS_D, alpha = LEWIS_ALPHA, center = FALSE, ess_min = 0)
got <- LEWIS::shd(unname(fit$D_hat), unname(rpdag(view0$dag, view0$nG_star)))
check(sprintf("calibration, rho = 0, d = 0.03: SHD %d on dataset 1 (stored 6)", got),
      got == 6L)

# ---- 4. RPDAG, SHD and CPDAG --------------------------------------------------
cat("4. RPDAG, SHD and CPDAG\n")
check("the RPDAG with all nodes non-Gaussian is the DAG",
      identical(unname(rpdag(inst$dag, seq_len(10L))), inst$dag))
check("the SHD of a graph from itself is 0", LEWIS::shd(truth, truth) == 0L)
if (requireNamespace("pcalg", quietly = TRUE)) {
  check("the CPDAG of pcalg equals the RPDAG with no non-Gaussian node",
        identical(to_cpdag(inst$dag), unname(rpdag(inst$dag))))
} else {
  cat("  [skipped] the package pcalg is not installed\n")
}

# ---- 5. Evaluation of the comparative study ----------------------------------
# Stored row of results_lewis_complexity.csv for dataset 1: SHD 7, SID bounds
# 10 and 14, 5 edges.
cat("5. Evaluation of the comparative study\n")
if (requireNamespace("pcalg", quietly = TRUE) && requireNamespace("SID", quietly = TRUE)) {
  source(file.path(COMMON_DIR, "method_wrappers.R"))
  source(file.path(COMMON_DIR, "scoring.R"))
  X <- draw_dataset(view, n = 50L, rep = 1L, regime_tag = "mixture")
  out <- run_method("lewis_complexity", X)
  row <- score_estimate("lewis_complexity", out, view)
  check(sprintf("LEWIS (complexity), dataset 1: SHD %g, SID %g to %g, %d edges (stored 7, 10 to 14, 5)",
                row$shd_rpdag, row$sid_lower, row$sid_upper, row$n_edges_est),
        row$shd_rpdag == 7 && row$sid_lower == 10 && row$sid_upper == 14 &&
          row$n_edges_est == 5L)
  X <- draw_dataset(view, n = 50L, rep = 1L, regime_tag = "mixture")
  out <- run_method("pc", X)
  row <- score_estimate("pc", out, view)
  check(sprintf("PC, dataset 1: SHD %g, SID %g to %g (stored 9, 19 to 19)",
                row$shd_rpdag, row$sid_lower, row$sid_upper),
        row$shd_rpdag == 9 && row$sid_lower == 19 && row$sid_upper == 19)
} else {
  cat("  [skipped] the packages pcalg and SID are needed\n")
}

# ---- 6. Uncertainty quantification --------------------------------------------
# The selected dataset of Section 7.3 is seed 2.
cat("6. Uncertainty quantification on the selected dataset\n")
source(file.path(.here, "04_uq", "uq_common.R"))
res <- run_uq(2L)
s <- res$summary
check(sprintf("acceptance rate %.4f (paper 6.4%%)", s$acceptance_rate),
      abs(s$acceptance_rate - 0.0642) < 1e-4)
check(sprintf("%d distinct DAGs and %d distinct RPDAGs (paper 32 and 24)",
              s$distinct_dags, s$distinct_rpdags),
      s$distinct_dags == 32L && s$distinct_rpdags == 24L)
check(sprintf("D_hat has %d edges, SHD %d from the truth, mass %.6f (paper 8, 2, 42.8%%)",
              s$D_hat_edges, s$shd_Dhat_truth, res$uq$mass[res$uq$is_D_hat][1]),
      s$D_hat_edges == 8L && s$shd_Dhat_truth == 2L &&
        abs(res$uq$mass[res$uq$is_D_hat][1] - 0.427625) < 1e-9)
check(sprintf("nG_hat is {%s} (paper {1, 2, 3, 7, 8, 10})", s$nG_hat),
      identical(s$nG_hat, "1,2,3,7,8,10"))
check(sprintf("credible ball radius %d (paper 3)", s$eps_alpha), s$eps_alpha == 3L)

# ---- Result ---------------------------------------------------------------------
cat(sprintf("\nRun time: %.0f seconds.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
if (length(failed) > 0L) {
  cat(sprintf("%d checks failed:\n", length(failed)))
  for (f in failed) cat("  -", f, "\n")
  quit(status = 1)
}
cat("All checks passed.\n")
