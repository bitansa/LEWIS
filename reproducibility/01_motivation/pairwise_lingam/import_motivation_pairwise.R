# =============================================================================
# import_motivation_pairwise.R   (Figure 1, Pairwise-LiNGAM stage 3 of 3)
#
# Reads manifest.csv, the datasets in inputs/ and the likelihood-ratio
# matrices that MATLAB wrote to outputs/ (stage 2). For every regime, n and
# dataset it builds the estimate of Pairwise-LiNGAM, evaluates it against
# the true RPDAG of the regime, and adds a row with method =
# "pairwise_lingam" to ../results/motivation_results.csv.
#
# THE ESTIMATE (Appendix G.1), as in 02_comparison/pairwise_lingam:
#   1. Skeleton: partial correlations from the inverse covariance matrix,
#      with a Fisher z-test at level 0.05 for each pair.
#   2. Orientation: each edge {i, j} of the skeleton is oriented by the sign
#      of the likelihood ratio LR[i, j] of pwling: i -> j if it is positive,
#      j -> i if it is negative. An exact zero leaves the edge undirected.
# The method assumes that all errors are non-Gaussian, so the graph is used
# as returned. It may contain a directed cycle.
#
# WHY inputs/ IS NEEDED. The skeleton is computed from the dataset. The
# script reads the dataset from inputs/, the same file that MATLAB read.
# If inputs/ is missing, run export_motivation_for_matlab.R again. It
# regenerates the same files from their seeds in a few minutes.
#
# The properties of the errors (columns ks_pvalue_median and
# excess_kurtosis) belong to the dataset. They are copied from the rows that
# motivation_experiment.R wrote for the same dataset, so that script must
# run first. The run time comes from outputs/runtimes.csv, matched on
# LR_file.
#
# The script can be restarted. Rows already present without failure are
# skipped. A failed row is computed again and replaces the old row. The
# script exits with a non-zero status if any dataset fails.
#
# Run from any folder, for example from the root of the repository:
#     Rscript reproducibility/01_motivation/pairwise_lingam/import_motivation_pairwise.R
# On a cluster, see slurm/motivation_pairwise_import.slurm.
# =============================================================================

ALPHA_SKEL <- 0.05
METHOD_NAME <- "pairwise_lingam"

# ---- Folder of this script and shared code ----------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- file.path(.here, "..", "..", "common")
suppressPackageStartupMessages({
  library(LEWIS)
  source(file.path(COMMON_DIR, "scenario_instances.R"))
  source(file.path(COMMON_DIR, "paper12_instance.R"))
  source(file.path(COMMON_DIR, "graph_tools.R"))
  source(file.path(.here, "..", "motivation_config.R"))
})

# ---- Inputs and outputs ------------------------------------------------------
# Stop at once if an input of this stage is missing.
if (!file.exists(file.path(.here, "manifest.csv"))) {
  stop("The file manifest.csv is missing. Run export_motivation_for_matlab.R first.")
}
manifest <- read.csv(file.path(.here, "manifest.csv"), stringsAsFactors = FALSE)
in_dir <- file.path(.here, "inputs")
out_dir <- file.path(.here, "outputs")
csv_out <- file.path(.here, "..", "results", "motivation_results.csv")
if (!file.exists(csv_out)) {
  stop("Run motivation_experiment.R first. The file ", csv_out, " is missing.")
}
if (!dir.exists(in_dir) || length(list.files(in_dir, pattern = "_X\\.csv$")) == 0L) {
  stop("The folder inputs/ is missing or empty. Run ",
       "export_motivation_for_matlab.R to regenerate it.")
}
if (!dir.exists(out_dir)) {
  stop("The folder outputs/ is missing. Run MATLAB (stage 2) or get the ",
       "outputs from the Zenodo archive.")
}

# ---- Run times of MATLAB --------------------------------------------------------
# runtimes.csv is appended by every MATLAB run, so a file can have several
# rows. The last row with status OK belongs to the file on disk.
rt_path <- file.path(out_dir, "runtimes.csv")
runtimes <- if (file.exists(rt_path)) read.csv(rt_path, stringsAsFactors = FALSE) else
  data.frame(LR_file = character(), runtime_s = numeric(), status = character())
runtime_of <- function(lr_file) {
  sel <- which(runtimes$LR_file == lr_file & runtimes$status == "OK")
  if (length(sel) > 0L) as.numeric(runtimes$runtime_s[sel[length(sel)]]) else NA_real_
}

# ---- References of the evaluation ---------------------------------------------
REG_BY_NAME <- setNames(REGIMES, REGIME_NAMES)
stopifnot(all(manifest$experiment %in% REGIME_NAMES))
TRUE_CPDAG <- to_cpdag(VIEW12$dag)
REFS <- lapply(REG_BY_NAME, function(rg) unname(rpdag(VIEW12$dag, rg$nG)))

# ---- Skeleton from partial correlations ----------------------------------------
# The generalized inverse is used only if the covariance matrix is singular.
build_skeleton_pcorr <- function(X, alpha_skel = 0.05) {
  n_obs <- nrow(X)
  p <- ncol(X)
  Sigma <- cov(X)
  Omega <- tryCatch(solve(Sigma), error = function(e) MASS::ginv(Sigma))
  # Partial correlations.
  pcorr <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) if (i != j) {
    pcorr[i, j] <- -Omega[i, j] / sqrt(Omega[i, i] * Omega[j, j])
  }
  # Fisher z-test of each pair.
  z_stat <- 0.5 * log((1 + pcorr) / (1 - pcorr))
  se_z <- 1 / sqrt(max(n_obs - p - 1, 1))
  pvals <- 2 * pnorm(-abs(z_stat / se_z))
  diag(pvals) <- 1
  skel <- (pvals < alpha_skel)
  diag(skel) <- FALSE
  skel
}

# ---- Orientation by the sign of the likelihood ratio --------------------------
orient_by_LR <- function(LR, skeleton) {
  p <- nrow(LR)
  A <- matrix(0L, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) if (i < j) {
    if (skeleton[i, j] || skeleton[j, i]) {
      if (LR[i, j] > 0) {
        A[i, j] <- 1L
      } else if (LR[i, j] < 0) {
        A[j, i] <- 1L
      } else {
        A[i, j] <- 1L
        A[j, i] <- 1L
      }
    }
  }
  A
}

# ---- Rows already imported -----------------------------------------------------
# A row counts as done only if it did not fail.
prev <- read.csv(csv_out, stringsAsFactors = FALSE)
have <- prev$method == METHOD_NAME
prev_failed <- as.logical(prev$skipped)
prev_key <- sprintf("%s|n%d|r%d", prev$regime, prev$n, prev$rep)
done <- prev_key[have & !prev_failed]

# ---- Evaluate every dataset -------------------------------------------------------
rows <- list()
for (k in seq_len(nrow(manifest))) {
  mr <- manifest[k, ]
  if (sprintf("%s|n%d|r%d", mr$experiment, mr$n, mr$rep_id) %in% done) next
  rg <- REG_BY_NAME[[mr$experiment]]
  ref <- REFS[[mr$experiment]]

  # Build the estimate. A failure is printed with its message and flagged.
  est <- tryCatch({
    lr_path <- file.path(out_dir, mr$LR_file)
    if (!file.exists(lr_path)) stop("the file ", mr$LR_file, " is missing")
    X <- as.matrix(read.csv(file.path(in_dir, mr$X_file), header = FALSE))
    LR <- as.matrix(read.csv(lr_path, header = FALSE))
    storage.mode(LR) <- "double"
    stopifnot(nrow(X) == mr$n, ncol(X) == P12, dim(LR) == c(P12, P12))
    orient_by_LR(LR, build_skeleton_pcorr(X, ALPHA_SKEL))
  }, error = function(e) {
    message(sprintf("Dataset %s failed: %s", mr$LR_file, conditionMessage(e)))
    NULL
  })
  bad <- is.null(est)

  # Properties of the errors, copied from the rows of the same dataset.
  src <- prev[prev$regime == rg$regime & prev$n == mr$n & prev$rep == mr$rep_id, ]
  if (nrow(src) == 0L) {
    stop("No row of motivation_experiment.R for ", mr$LR_file,
         ". Run motivation_experiment.R first.")
  }

  rows[[length(rows) + 1L]] <- data.frame(
    regime = rg$regime, family = rg$family, k_ng = rg$k_ng, df = rg$df,
    p = P12, n = mr$n, rep = mr$rep_id, method = METHOD_NAME,
    shd_rpdag = if (bad) NA_real_ else LEWIS::shd(est, ref),
    shd_vs_cpdag = if (bad) NA_real_ else LEWIS::shd(est, TRUE_CPDAG),
    n_edges_est = if (bad) NA_integer_ else n_edges_partial(est),
    runtime_s = runtime_of(mr$LR_file), skipped = bad,
    ks_pvalue_median = src$ks_pvalue_median[1],
    excess_kurtosis = src$excess_kurtosis[1],
    stringsAsFactors = FALSE)
  if (k %% 500L == 0L) cat(sprintf("Evaluated %d of %d datasets\n", k, nrow(manifest)))
}

# ---- Write the new rows ------------------------------------------------------------
# Append the new rows. Failed rows of an earlier run are replaced.
n_failed <- 0L
if (length(rows) > 0L) {
  block <- do.call(rbind, rows)
  n_failed <- sum(block$skipped)
  # Remove the failed rows of earlier runs that the new rows replace.
  block_key <- sprintf("%s|n%d|r%d", block$regime, block$n, block$rep)
  drop_old <- have & prev_failed & prev_key %in% block_key
  if (any(drop_old)) {
    write.csv(prev[!drop_old, , drop = FALSE], csv_out, row.names = FALSE)
    cat(sprintf("Replaced %d failed rows of an earlier run.\n", sum(drop_old)))
  }
  write.table(block, csv_out, sep = ",", row.names = FALSE,
              col.names = FALSE, append = TRUE)
  cat(sprintf("Added %d Pairwise-LiNGAM rows, %d of them failed.\n",
              nrow(block), n_failed))
} else {
  cat("Nothing to import. All rows are present.\n")
}

# ---- Final check ---------------------------------------------------------------------
# Every dataset of the manifest has one row, and no row failed.
res <- read.csv(csv_out, stringsAsFactors = FALSE)
res <- res[res$method == METHOD_NAME & res$regime %in% REGIME_NAMES, ]
rkey <- sprintf("%s|n%d|r%d", res$regime, res$n, res$rep)
mkey <- sprintf("%s|n%d|r%d", manifest$experiment, manifest$n, manifest$rep_id)
n_missing <- sum(!(mkey %in% rkey))
cat(sprintf("Pairwise-LiNGAM rows: %d, missing %d, failed %d, duplicated %d.\n",
            nrow(res), n_missing, sum(as.logical(res$skipped)), sum(duplicated(rkey))))
if (n_failed > 0L || n_missing > 0L || any(as.logical(res$skipped)) || anyDuplicated(rkey)) {
  message("The import is not complete. Fix the errors above and run the script again.")
  quit(status = 1)
}
cat("All Pairwise-LiNGAM rows are present and none failed.\n")
