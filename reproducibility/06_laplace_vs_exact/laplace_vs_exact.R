# =============================================================================
# laplace_vs_exact.R
#
# Appendix F of the paper, Table 1: the chain of LEWIS with the Laplace
# approximation (L) against the exact chain (E).
#
# DESIGN (Appendix F).
#   * Data: the ER graphs of the comparative study with rho = 0.5 and
#     s = 1, and the 13 pairs (p, n) of Table 1. The datasets are generated
#     by draw_dataset() with the regime tag "mixture", so they are the
#     datasets of the comparative study where the pairs (p, n) agree.
#   * 100 datasets for each pair (p, n).
#   * For each dataset, both chains start at the estimate of the TL method,
#     computed once per dataset, and run for the same number of iterations
#     (lewis_schedule(p)): 10,000 for p <= 20, 62,500 for p = 50 and
#     250,000 for p = 100.
#   * Both priors: the complexity prior with alpha = 0.75 and d = 0.03, and
#     the uniform prior (d = 0).
#   * The exact chain uses 300 importance draws per node for p <= 20 and
#     150 for p >= 50.
#   * For each chain, the highest-scoring DAG is compared with the true DAG
#     through their CPDAGs (no post-processing).
#
# ORDER OF THE RANDOM NUMBERS within one dataset, as in the runs of the
# paper: the data (draw_dataset() sets the seed), the start graph, then the
# Laplace chain and the exact chain under the uniform prior, then the
# Laplace chain and the exact chain under the complexity prior.
#
# OUTPUT: results/exact_vs_laplace.csv, one row per dataset, prior and
# chain, with the columns
#     rep, method ("Laplace" or "Exact"), prior, shd_cpdag (SHD between the
#     CPDAG of the highest-scoring DAG and the CPDAG of the true DAG),
#     shd_between (SHD between the CPDAGs of the two chains), sel_edges (the
#     highest-scoring DAG), runtime_s (computing time of the chain in
#     seconds), p, n, n_over_p, config, m_is, n_iter.
# The file is written after each pair (p, n). A new run skips the datasets
# that are complete in the file. The script exits with status 1 if any
# dataset failed.
#
# RUN. From the root of the repository (or from any other folder):
#     MC_CORES=100 Rscript reproducibility/06_laplace_vs_exact/laplace_vs_exact.R
# The datasets of one pair (p, n) run in parallel on MC_CORES cores. Without
# MC_CORES, all cores but one are used. slurm/laplace_vs_exact.slurm runs the
# same command. Then run make_table_lve.R for Table 1.
# =============================================================================

suppressPackageStartupMessages({
  library(LEWIS)
  library(parallel)
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
source(file.path(COMMON_DIR, "settings.R"))
source(file.path(COMMON_DIR, "scenario_instances.R"))
source(file.path(COMMON_DIR, "graph_tools.R"))
source(file.path(COMMON_DIR, "laplace_exact.R"))

# Compile the exact chain once, before the workers are forked.
load_exact_engine()

# ---- Design -------------------------------------------------------------------
CONFIGS <- data.frame(
  p = c(10L, 10L, 10L, 10L, 20L, 20L, 20L, 20L, 50L, 50L, 50L, 100L, 100L),
  n = c(100L, 250L, 500L, 1000L, 100L, 250L, 500L, 1000L, 250L, 500L, 1000L,
        500L, 1000L))
N_REPS <- 100L
PRIORS <- c("uniform", "complexity")
METHODS <- c("Laplace", "Exact")

# Importance draws per node, and the number of iterations of both chains.
m_is_for_p <- function(p) if (p <= 20L) 300L else 150L
n_iter_for_p <- function(p) lewis_schedule(p)$n_iter

# Number of parallel workers: MC_CORES, or all cores but one.
MC_CORES <- {
  env <- suppressWarnings(as.integer(Sys.getenv("MC_CORES", "")))
  if (is.finite(env) && env >= 1L) env else max(1L, detectCores() - 1L)
}

# ---- Output file and resume ----------------------------------------------------
res_dir <- file.path(.here, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
csv_out <- file.path(res_dir, "exact_vs_laplace.csv")

# Key of one row.
.key <- function(p, n, rep, prior, method) {
  sprintf("p%d|n%d|r%d|%s|%s", as.integer(p), as.integer(n), as.integer(rep),
          prior, method)
}

# Rows already in the file. A dataset is complete when its four rows exist.
prev <- if (file.exists(csv_out)) read.csv(csv_out, stringsAsFactors = FALSE) else NULL
if (!is.null(prev)) {
  stopifnot(!anyDuplicated(.key(prev$p, prev$n, prev$rep, prev$prior, prev$method)))
  cat(sprintf("Resume: %d rows are in %s\n", nrow(prev), basename(csv_out)))
}
done_keys <- if (is.null(prev)) character(0) else
  .key(prev$p, prev$n, prev$rep, prev$prior, prev$method)
rep_done <- function(p, n, r) {
  all(.key(p, n, r, rep(PRIORS, each = 2), rep(METHODS, 2)) %in% done_keys)
}

# ---- One dataset -----------------------------------------------------------------
run_dataset <- function(p, n, r, n_iter, m_is) {
  # The dataset of the comparative study.
  view <- instance_view(build_instance("ER", p), rho = 0.5, effect_scale = 1.0)
  X <- draw_dataset(view, n = n, rep = r, regime_tag = "mixture")
  truth <- view$dag

  # The start graph of both chains.
  start <- tl_start_dag(X)

  rows <- list()
  for (prior in PRIORS) {
    # The Laplace chain.
    t1 <- proc.time()[["elapsed"]]
    lap <- laplace_chain(X, start, prior, n_iter)
    t_lap <- proc.time()[["elapsed"]] - t1

    # The exact chain.
    t2 <- proc.time()[["elapsed"]]
    ex <- exact_chain(X, start, prior, n_iter, m_is)
    t_ex <- proc.time()[["elapsed"]] - t2

    # SHD between the CPDAGs of the two highest-scoring DAGs.
    agree <- cpdag_shd(lap$dag, ex$dag)
    rows[[length(rows) + 1L]] <- data.frame(
      rep = r, method = "Laplace", prior = prior,
      shd_cpdag = cpdag_shd(lap$dag, truth), shd_between = agree,
      sel_edges = edge_string(lap$dag), runtime_s = t_lap,
      stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      rep = r, method = "Exact", prior = prior,
      shd_cpdag = cpdag_shd(ex$dag, truth), shd_between = agree,
      sel_edges = edge_string(ex$dag), runtime_s = t_ex,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$p <- p
  out$n <- n
  out$n_over_p <- round(n / p, 2)
  out$config <- sprintf("p=%d, n=%d", p, n)
  out$m_is <- m_is
  out$n_iter <- n_iter
  out
}

# ---- Main loop over the pairs (p, n) --------------------------------------------
all_rows <- if (is.null(prev)) list() else list(prev)
n_failed <- 0L
for (k in seq_len(nrow(CONFIGS))) {
  p <- CONFIGS$p[k]
  n <- CONFIGS$n[k]
  n_iter <- n_iter_for_p(p)
  m_is <- m_is_for_p(p)
  to_run <- Filter(function(r) !rep_done(p, n, r), seq_len(N_REPS))
  cat(sprintf("p = %d, n = %d: %d datasets to run (iterations %d, importance draws %d)\n",
              p, n, length(to_run), n_iter, m_is))
  if (length(to_run) == 0L) next

  # Run the datasets in parallel. A failure returns its message.
  res <- mclapply(to_run, function(r) {
    tryCatch(run_dataset(p, n, r, n_iter, m_is),
             error = function(e) conditionMessage(e))
  }, mc.cores = MC_CORES, mc.preschedule = FALSE)

  # Report the failures and keep the complete datasets.
  ok <- vapply(res, is.data.frame, logical(1))
  for (i in which(!ok)) {
    msg <- if (is.character(res[[i]])) res[[i]] else "no result (the worker stopped)"
    cat(sprintf("  FAILED: p = %d, n = %d, dataset %d: %s\n", p, n, to_run[i], msg))
  }
  n_failed <- n_failed + sum(!ok)
  if (any(ok)) all_rows <- c(all_rows, res[ok])

  # Write the file after each pair (p, n), through a temporary file. Nothing
  # is written while no dataset is complete.
  if (length(all_rows) > 0L) {
    tab <- do.call(rbind, all_rows)
    tmp <- paste0(csv_out, ".tmp")
    write.csv(tab, tmp, row.names = FALSE)
    file.rename(tmp, csv_out)
  }
}

# ---- Check that every dataset is complete ---------------------------------------
if (!file.exists(csv_out)) {
  cat("No dataset is complete. Fix the errors above and run the script again.\n")
  quit(status = 1)
}
tab <- read.csv(csv_out, stringsAsFactors = FALSE)
expected <- nrow(CONFIGS) * N_REPS * length(PRIORS) * length(METHODS)
cat(sprintf("%d rows of %d expected. %d datasets failed in this run.\n",
            nrow(tab), expected, n_failed))
if (n_failed > 0L || nrow(tab) != expected) {
  cat("Some datasets are missing. Run the script again to complete them.\n")
  quit(status = 1)
}
cat("All datasets are complete. Next: Rscript reproducibility/06_laplace_vs_exact/make_table_lve.R\n")
