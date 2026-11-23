# =============================================================================
# import_pairwise_results.R   (Pairwise-LiNGAM, stage 3 of 3)
#
# Reads manifest.csv, the datasets in inputs/ (stage 1) and the
# likelihood-ratio matrices in outputs/ (stage 2, MATLAB). It builds the
# Pairwise-LiNGAM estimate of every dataset, evaluates it as every other
# method, and writes
#
#     ../results/results_pairwise_lingam.csv
#
# with the columns of the other result files.
#
# ESTIMATE (Appendix G.1, following Hyvarinen and Smith 2013):
#   1. Skeleton. The partial correlations are computed from the inverse of
#      the sample covariance. A pair is joined when the Fisher z test of its
#      partial correlation rejects at level 0.05.
#   2. Orientation. Each edge {i, j} of the skeleton is oriented by the sign
#      of the likelihood ratio LR[i, j] of pwling: i -> j if it is positive,
#      j -> i if it is negative. An exact zero leaves the edge undirected.
# Each edge is oriented on its own, so the graph may contain a directed
# cycle. The method assumes that all errors are non-Gaussian, so its graph
# is compared with the true RPDAG as it is.
#
# The run times come from outputs/runtimes.csv, matched on LR_file. When
# MATLAB ran a dataset more than once, the last successful run is used.
#
# FAILURES. A dataset whose LR file is missing or cannot be read gives a row
# with skipped = TRUE, and its message is printed. An error while
# evaluating a dataset is printed and gives no row. The script then writes
# nothing and exits with status 1.
#
# Run from any folder, for example from the root of the repository:
#     MC_CORES=100 Rscript reproducibility/02_comparison/pairwise_lingam/import_pairwise_results.R
# On a cluster, see slurm/import_pairwise_comparison.slurm. The SID at
# p = 50 and p = 100 takes hours on one core, so use many cores.
# =============================================================================

ALPHA_SKEL <- 0.05
METHOD_NAME <- "pairwise_lingam"

# The datasets are evaluated in parallel.
suppressPackageStartupMessages(library(parallel))

# ---- Folder of this script and shared code ----------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- file.path(.here, "..", "..", "common")
source(file.path(COMMON_DIR, "scenario_instances.R"))
source(file.path(COMMON_DIR, "method_wrappers.R"))
source(file.path(COMMON_DIR, "scoring.R"))
source(file.path(.here, "..", "comparison_cells.R"))

# ---- Inputs -----------------------------------------------------------------------
# Stop at once if an input of this stage is missing.
if (!file.exists(file.path(.here, "manifest.csv"))) {
  stop("The file manifest.csv is missing. Run export_datasets_for_matlab.R first.")
}
manifest <- read.csv(file.path(.here, "manifest.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == nrow(comparison_cells()) * N_REPS)
out_dir <- file.path(.here, "outputs")
in_dir <- file.path(.here, "inputs")
res_dir <- file.path(.here, "..", "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
csv_out <- file.path(res_dir, "results_pairwise_lingam.csv")

# The skeleton is built from the exported datasets, so inputs/ is needed.
if (!dir.exists(in_dir) || length(list.files(in_dir, pattern = "_X\\.csv$")) == 0L) {
  stop("The folder inputs/ is missing or empty. Run export_datasets_for_matlab.R first.")
}
rt_path <- file.path(out_dir, "runtimes.csv")
if (!file.exists(rt_path)) stop("outputs/runtimes.csv is missing. Run stage 2 first.")
runtimes <- read.csv(rt_path, stringsAsFactors = FALSE)

# Run time of the last successful MATLAB run of one dataset.
runtime_of <- function(lr_file) {
  sel <- which(runtimes$LR_file == lr_file & runtimes$status == "OK")
  if (length(sel) > 0L) as.numeric(runtimes$runtime_s[sel[length(sel)]]) else NA_real_
}

# ---- Skeleton by Fisher z tests of the partial correlations ----------------------
build_skeleton_pcorr <- function(X, alpha_skel = 0.05) {
  n_obs <- nrow(X)
  p <- ncol(X)
  Sigma <- cov(X)
  # Inverse covariance. The generalized inverse is used if Sigma is singular.
  Omega <- tryCatch(solve(Sigma), error = function(e) {
    message("Singular covariance matrix, the generalized inverse is used: ",
            conditionMessage(e))
    MASS::ginv(Sigma)
  })
  # Partial correlations.
  pcorr <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) if (i != j) {
    pcorr[i, j] <- -Omega[i, j] / sqrt(Omega[i, i] * Omega[j, j])
  }
  # Fisher z test.
  z_stat <- 0.5 * log((1 + pcorr) / (1 - pcorr))
  se_z <- 1 / sqrt(max(n_obs - p - 1, 1))
  pvals <- 2 * pnorm(-abs(z_stat / se_z))
  diag(pvals) <- 1
  skel <- (pvals < alpha_skel)
  diag(skel) <- FALSE
  skel
}

# ---- Orientation by the sign of the likelihood ratio -------------------------------
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

# ---- Workers --------------------------------------------------------------------
# The environment variable MC_CORES, or all cores but one.
MC_CORES <- {
  env <- suppressWarnings(as.integer(Sys.getenv("MC_CORES", "")))
  if (is.finite(env) && env >= 1L) env else max(1L, detectCores() - 1L)
}
cat(sprintf("Evaluating %d datasets on %d workers.\n", nrow(manifest), MC_CORES))

# ---- One dataset ------------------------------------------------------------------
score_row <- function(k) {
  mr <- manifest[k, ]
  inst <- build_instance(mr$graph_type, mr$p)
  view <- instance_view(inst, rho = mr$rho, effect_scale = mr$effect_scale)
  rt <- runtime_of(mr$LR_file)

  # Build the estimate. A missing or unreadable file gives a flagged row.
  out <- tryCatch({
    lr_path <- file.path(out_dir, mr$LR_file)
    if (!file.exists(lr_path)) stop("the file ", mr$LR_file, " is missing")
    LR <- as.matrix(read.csv(lr_path, header = FALSE))
    storage.mode(LR) <- "double"
    X <- as.matrix(read.csv(file.path(in_dir, mr$X_file), header = FALSE))
    skel <- build_skeleton_pcorr(X, alpha_skel = ALPHA_SKEL)
    list(dag = orient_by_LR(LR, skel), runtime = rt, skipped = FALSE, error = "")
  }, error = function(e) {
    message(sprintf("Dataset %s: %s", mr$LR_file, conditionMessage(e)))
    list(dag = matrix(0L, mr$p, mr$p), runtime = rt, skipped = TRUE,
         error = conditionMessage(e))
  })

  # Evaluate the estimate.
  met <- score_estimate(METHOD_NAME, out, view)
  cbind(data.frame(experiment = mr$experiment, graph_type = mr$graph_type,
                   p = mr$p, n = mr$n, rho = mr$rho, effect_scale = mr$effect_scale,
                   rep = mr$rep_id, method = METHOD_NAME, stringsAsFactors = FALSE),
        met)
}

# ---- All datasets ----------------------------------------------------------------
# The slow datasets (p = 100) are spread over the workers by interleaving
# the manifest.
order_k <- order(seq_len(nrow(manifest)) %% MC_CORES, seq_len(nrow(manifest)))
rows <- mclapply(order_k, function(k) {
  tryCatch(score_row(k),
           error = function(e) sprintf("row %d of the manifest: %s", k, conditionMessage(e)))
}, mc.cores = MC_CORES, mc.preschedule = FALSE)

# ---- Report the failures and write the file ------------------------------------------
# The file is written only if every dataset gave a row. The rows are put
# back in the order of the manifest.
failed <- which(!vapply(rows, is.data.frame, logical(1)))
for (i in failed) {
  message(if (is.character(rows[[i]])) rows[[i]] else
    sprintf("row %d of the manifest: no result (the worker process was killed)", order_k[i]))
}
if (length(failed) > 0L) {
  message(sprintf("%d of %d datasets gave no result. No file was written.",
                  length(failed), nrow(manifest)))
  quit(status = 1)
}
res <- do.call(rbind, rows[order(order_k)])
write.csv(res, csv_out, row.names = FALSE)
n_skip <- sum(res$skipped)
cat(sprintf("Wrote %s (%d rows, %d with skipped = TRUE).\n", csv_out, nrow(res), n_skip))
if (n_skip > 0L) {
  message("Some datasets have no valid LR file. Run stage 2 again for them, then this script.")
  quit(status = 1)
}
