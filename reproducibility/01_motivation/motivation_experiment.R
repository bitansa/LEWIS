# =============================================================================
# motivation_experiment.R
#
# The computation of Figure 1 (design in Appendix G.1, "Design of
# Figure 1").
#
# For each of the 25 regimes of motivation_config.R, each sample size
# n in {100, 250, 500} and each of 100 datasets, the script generates the
# dataset on the 12-node DAG of the paper and runs nine methods on it, one
# after the other:
#
#     LEWIS (complexity), LEWIS (uniform), ICA-LiNGAM, DirectLiNGAM, RCD,
#     TL, NOTEARS, PC, GES
#
# Pairwise-LiNGAM runs on the same datasets in MATLAB (see the folder
# pairwise_lingam/), and its rows are added to the same results file.
#
# RANDOM NUMBERS. draw_dataset() sets the seed of the dataset. The nine
# methods then run in the order above and share one random stream: the
# chain of LEWIS (complexity) starts from the state left by the data
# generation, the chain of LEWIS (uniform) from the state left by LEWIS
# (complexity), and the seed of ICA-LiNGAM is drawn from the state left by
# the two chains. The methods are therefore not seeded one by one. A
# separate seed for each method would change the results of every method
# except the first. The datasets and the results of the paper depend only
# on the regime name, n, the dataset number and this order of the methods.
#
# EVALUATION. For each method, the SHD between its graph (see
# est_partial_dag() in common/method_wrappers.R) and the true RPDAG of the
# regime (column shd_rpdag), and the SHD from the true CPDAG (column
# shd_vs_cpdag, a diagnostic that the paper does not use). The SID is not
# computed here. Two properties of the errors of each dataset are also
# stored: the median over the nodes of the Kolmogorov-Smirnov p-value
# against N(0, 1), and the pooled excess kurtosis.
#
# OUTPUT: results/motivation_results.csv, one row per regime, n, method and
# dataset. The script can be restarted. It skips every (regime, n, dataset)
# whose nine rows are present without failure. A dataset with a failed
# method is run again from the start, so that the random stream is the same.
#
# Run from any folder, for example from the root of the repository:
#     MC_CORES=100 Rscript reproducibility/01_motivation/motivation_experiment.R
# MC_CORES is the number of datasets run at the same time (default: all
# cores but one). On a cluster, see slurm/motivation.slurm.
# =============================================================================

suppressPackageStartupMessages(library(parallel))

# ---- Folder of this script and shared code ----------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- file.path(.here, "..", "common")
source(file.path(COMMON_DIR, "scenario_instances.R"))
source(file.path(COMMON_DIR, "paper12_instance.R"))
source(file.path(COMMON_DIR, "method_wrappers.R"))
source(file.path(.here, "motivation_config.R"))

# ---- The methods, in the order in which they run on each dataset ------------
MOT_METHODS <- c("lewis_complexity", "lewis_uniform",
                 "lingam", "direct_lingam", "rcd",
                 "tl_zhao", "notears", "pc", "ges")

# Stop at once if the software of a method is missing.
for (m in MOT_METHODS) require_method(m)

# ---- Number of cores ----------------------------------------------------------
# The environment variable MC_CORES, or all cores but one.
MC_CORES <- {
  env <- suppressWarnings(as.integer(Sys.getenv("MC_CORES", "")))
  if (is.finite(env) && env >= 1L) env else max(1L, detectCores() - 1L)
}

# ---- Output file ------------------------------------------------------------
# Create results/ if needed.
res_dir <- file.path(.here, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
csv_out <- file.path(res_dir, "motivation_results.csv")

# Print the size of the study.
cat(sprintf("12-node DAG, n in {%s}, %d regimes, %d methods, %d datasets, %d cores\n",
            paste(N_VALUES, collapse = ", "), length(REGIMES),
            length(MOT_METHODS), N_REPS, MC_CORES))

# ---- References of the evaluation ---------------------------------------------
# The true CPDAG, and the true RPDAG of each regime.
TRUE_CPDAG <- to_cpdag(VIEW12$dag)
REFS <- lapply(REGIMES, function(rg) unname(rpdag(VIEW12$dag, rg$nG)))

# ---- Properties of the errors -------------------------------------------------
# Median over the nodes of the Kolmogorov-Smirnov p-value against N(0, 1).
ks_median_pvalue <- function(Z) {
  pv <- apply(Z, 2, function(z) suppressWarnings(ks.test(z, "pnorm"))$p.value)
  median(pv)
}
# Pooled excess kurtosis of all the unit-variance errors.
excess_kurtosis <- function(Z) {
  z <- as.numeric(Z)
  z <- z - mean(z)
  mean(z^4) / mean(z^2)^2 - 3
}

# ---- Keys of the rows already present ---------------------------------------
# Only rows without failure count as done.
.key <- function(regime, n, method, rep) sprintf("%s|n%d|%s|r%d", regime, n, method, rep)
done_keys <- character(0)
if (file.exists(csv_out)) {
  prev <- read.csv(csv_out, stringsAsFactors = FALSE)
  needed <- c("regime", "n", "method", "rep", "shd_rpdag", "skipped")
  if (!all(needed %in% names(prev))) {
    stop("The file ", csv_out, " lacks the columns ",
         paste(setdiff(needed, names(prev)), collapse = ", "), ".")
  }
  ok <- !as.logical(prev$skipped)
  done_keys <- .key(prev$regime[ok], prev$n[ok], prev$method[ok], prev$rep[ok])
  cat(sprintf("Resume: %d rows present, %d without failure.\n", nrow(prev), sum(ok)))
}

# ---- One dataset: the nine methods on the same data ---------------------------
run_one <- function(ri, n, r) {
  rg <- REGIMES[[ri]]
  ref <- REFS[[ri]]
  # Generate the dataset and keep the unit-variance errors.
  X <- draw_motivation_dataset(rg, n = n, rep = r, return_errors = TRUE)
  Z <- attr(X, "Z")
  attr(X, "Z") <- NULL
  ks_med <- ks_median_pvalue(Z)
  exkurt <- excess_kurtosis(Z)

  # Run the methods in the fixed order and evaluate each estimate.
  rows <- lapply(MOT_METHODS, function(m) {
    out <- run_method(m, X)
    bad <- isTRUE(out$skipped)
    est <- if (bad) NULL else est_partial_dag(m, out$dag)
    data.frame(regime = rg$regime, family = rg$family,
               k_ng = rg$k_ng, df = rg$df,
               p = P12, n = n, rep = r, method = m,
               shd_rpdag = if (bad) NA_real_ else LEWIS::shd(est, ref),
               shd_vs_cpdag = if (bad) NA_real_ else LEWIS::shd(est, TRUE_CPDAG),
               n_edges_est = if (bad) NA_integer_ else n_edges_partial(est),
               runtime_s = out$runtime,
               skipped = bad,
               ks_pvalue_median = ks_med,
               excess_kurtosis = exkurt,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---- Main loop: settings outside, datasets in parallel inside ----------------
t_start <- Sys.time()
n_errors <- 0L
for (n in N_VALUES) {
  for (ri in seq_along(REGIMES)) {
    rg <- REGIMES[[ri]]
    # The datasets with at least one method not done.
    to_run <- Filter(function(r) !all(.key(rg$regime, n, MOT_METHODS, r) %in% done_keys),
                   seq_len(N_REPS))
    if (length(to_run) == 0L) next
    cat(sprintf("n = %d, regime %d of %d (%s), %d datasets\n",
                n, ri, length(REGIMES), rg$regime, length(to_run)))

    # Run the datasets. An error is printed with its message and the
    # dataset is left for the next run.
    blocks <- mclapply(to_run, function(r) {
      tryCatch(run_one(ri, n, r), error = function(e) {
        message(sprintf("Regime %s, n = %d, dataset %d failed: %s",
                        rg$regime, n, r, conditionMessage(e)))
        NULL
      })
    }, mc.cores = MC_CORES, mc.preschedule = FALSE)
    n_lost <- sum(vapply(blocks, is.null, logical(1)))
    if (n_lost > 0L) {
      message(sprintf("Regime %s, n = %d: %d of %d datasets gave no result.",
                      rg$regime, n, n_lost, length(blocks)))
      n_errors <- n_errors + n_lost
    }
    blocks <- Filter(Negate(is.null), blocks)
    if (length(blocks) == 0L) next

    # Keep the new rows whose key is not done yet.
    block <- do.call(rbind, blocks)
    bkey <- .key(block$regime, block$n, block$method, block$rep)
    block <- block[!(bkey %in% done_keys), , drop = FALSE]
    if (nrow(block) == 0L) next

    # Remove failed rows of earlier runs that the new rows replace.
    if (file.exists(csv_out)) {
      prev <- read.csv(csv_out, stringsAsFactors = FALSE)
      pkey <- .key(prev$regime, prev$n, prev$method, prev$rep)
      drop <- as.logical(prev$skipped) &
        pkey %in% .key(block$regime, block$n, block$method, block$rep)
      if (any(drop)) {
        write.csv(prev[!drop, , drop = FALSE], csv_out, row.names = FALSE)
      }
    }

    # Append the new rows.
    new_file <- !file.exists(csv_out)
    write.table(block, csv_out, sep = ",", row.names = FALSE,
                col.names = new_file, append = !new_file)
    ok <- !block$skipped
    done_keys <- c(done_keys, .key(block$regime[ok], block$n[ok],
                                   block$method[ok], block$rep[ok]))
    cat(sprintf("  written: %d rows, %d failed methods (%.1f h elapsed)\n",
                nrow(block), sum(!ok),
                as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
  }
}

# ---- Final check: every row present, none failed -------------------------------
if (!file.exists(csv_out)) {
  message("No result was written. Fix the errors above and run the script again.")
  quit(status = 1)
}
res <- read.csv(csv_out, stringsAsFactors = FALSE)
res <- res[res$method %in% MOT_METHODS, ]
expected <- expand.grid(regime = REGIME_NAMES, n = N_VALUES, method = MOT_METHODS,
                        rep = seq_len(N_REPS), stringsAsFactors = FALSE)
ekey <- .key(expected$regime, expected$n, expected$method, expected$rep)
rkey <- .key(res$regime, res$n, res$method, res$rep)
n_missing <- sum(!(ekey %in% rkey))
n_failed <- sum(as.logical(res$skipped))
n_dup <- sum(duplicated(rkey))
cat(sprintf("Rows expected %d, missing %d, failed %d, duplicated %d.\n",
            length(ekey), n_missing, n_failed, n_dup))
if (n_errors > 0L || n_missing > 0L || n_failed > 0L || n_dup > 0L) {
  message("The run is not complete. Fix the errors above and run the script again.")
  quit(status = 1)
}
cat("All rows are present and no method failed.\n")
