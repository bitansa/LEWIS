# =============================================================================
# measure_lewis_runtime.R
#
# Measures the computing time of LEWIS that Appendix E.1 ("Computing time")
# reports: the median time of a complete run in the settings of Section 7.2,
# and the share of the time taken by the TL initialization at
# (p, n) = (100, 1000).
#
# A complete run is one call of lewis(): the TL initialization, the sampling
# and the post-processing. The script uses the fixed instances and the
# datasets of the comparative study (rho = 0.5, s = 1, regime tag
# "mixture"), the default schedule of lewis_schedule(p), and one dataset per
# core. It does not change the result files of the comparative study.
#
# Output: results/lewis_runtime_<prior>.csv, one row per dataset, with
#     graph, p, n, rep, prior
#     runtime_s        time of the call of lewis(), in seconds
#     warm_start_s     time of LEWIS::tl_lingam on the same dataset
#     warm_start_frac  warm_start_s / runtime_s
#     ms_per_iter      time per iteration net of the TL initialization,
#                      1000 * (runtime_s - warm_start_s) / n_iter
#     n_iter           number of iterations
#     shd              SHD of D_hat from the true RPDAG
#     edges            number of edges of D_hat
#     n_nG_hat         number of nodes in nG_hat
#     n_unique_dags    number of distinct DAGs after the burn-in
#     acceptance_rate  acceptance rate of the chain
#     engine           software that ran LEWIS
#
# TIME OF THE TL INITIALIZATION. lewis() does not report the time of its TL
# initialization separately. The script therefore times one more call of
# LEWIS::tl_lingam on the same dataset, AFTER the call of lewis(). The
# random stream that lewis() uses is then the same as in the comparative
# study, so the SHD column can be checked against it (see the end of the
# script).
#
# Usage, from any folder, for example from the root of the repository:
#     MC_CORES=100 Rscript reproducibility/02_comparison/measure_lewis_runtime.R
# Environment variables:
#     MC_CORES  number of workers (default: all cores but one)
#     N_REPS    number of datasets per setting (default 100)
#     CELLS     settings as "pxn" pairs (default: every setting that the
#               paper quotes, "10x50,10x500,10x1000,20x500,50x500,100x1000")
#     PRIOR     "complexity" (default, the prior of the paper) or "uniform"
#
# Resume: the rows already in the output file are skipped.
#
# The check at the end compares the SHD with results/results_lewis_<prior>.csv
# of the comparative study, if that file is present (from run_one_method.R
# or from the Zenodo archive). The script exits with status 1 if a dataset
# fails or an SHD differs.
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
suppressPackageStartupMessages({
  source(file.path(COMMON_DIR, "scenario_instances.R"))
  source(file.path(COMMON_DIR, "method_wrappers.R"))
  source(file.path(COMMON_DIR, "scoring.R"))
})
require_method("lewis_complexity")
require_method("tl_zhao")

# ---- Options ------------------------------------------------------------------
# Read the environment variables described in the header.
N_REPS <- as.integer(Sys.getenv("N_REPS", "100"))
PRIOR <- Sys.getenv("PRIOR", "complexity")
stopifnot(PRIOR %in% c("complexity", "uniform"))
MC_CORES <- {
  env <- suppressWarnings(as.integer(Sys.getenv("MC_CORES", "")))
  if (is.finite(env) && env >= 1L) env else max(1L, detectCores() - 1L)
}
ENGINE <- paste("LEWIS", as.character(utils::packageVersion("LEWIS")))

# The settings, both graph types. Each entry is "pxn".
CELL_SPEC <- Sys.getenv("CELLS", "10x50,10x500,10x1000,20x500,50x500,100x1000")
cells <- do.call(rbind, lapply(strsplit(CELL_SPEC, ",", fixed = TRUE)[[1]], function(s) {
  pn <- as.integer(strsplit(trimws(s), "x", fixed = TRUE)[[1]])
  expand.grid(graph = c("ER", "SF"), p = pn[1], n = pn[2], stringsAsFactors = FALSE)
}))
cells <- cells[order(cells$p, cells$n, cells$graph), ]
RHO <- 0.5
EFFECT <- 1.0

# ---- Output file and resume -----------------------------------------------------
# Create results/ if needed and read the rows already measured.
OUT_DIR <- file.path(.here, "results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_CSV <- file.path(OUT_DIR, sprintf("lewis_runtime_%s.csv", PRIOR))
existing <- if (file.exists(OUT_CSV)) read.csv(OUT_CSV, stringsAsFactors = FALSE) else NULL
done_key <- if (!is.null(existing)) {
  paste(existing$graph, existing$p, existing$n, existing$rep, sep = "|")
} else {
  character(0)
}
all_rows <- if (is.null(existing)) list() else list(existing)

# Print the plan of the run.
cat(sprintf("Run time of LEWIS (%s), prior %s, %d datasets per setting, %d workers\n",
            ENGINE, PRIOR, N_REPS, MC_CORES))
cat(sprintf("Output: %s\n\n", OUT_CSV))

# ---- One dataset --------------------------------------------------------------
time_one <- function(view, gtype, p, n, r, ref, n_iter) {
  # Generate the dataset. draw_dataset() sets the seed of this dataset.
  X <- draw_dataset(view, n, r)
  # The measured call. Nothing runs between draw_dataset() and lewis(), as
  # in the comparative study.
  t0 <- Sys.time()
  fit <- if (PRIOR == "complexity") {
    lewis(X, d = LEWIS_D, alpha = LEWIS_ALPHA, center = FALSE, ess_min = 0)
  } else {
    lewis(X, d = 0, center = FALSE, ess_min = 0)
  }
  total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # The TL initialization, timed by a separate call after lewis().
  t1 <- Sys.time()
  LEWIS::tl_lingam(as.matrix(X))
  ws <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  D_hat <- unname(fit$D_hat)
  data.frame(graph = gtype, p = p, n = n, rep = r, prior = PRIOR,
             runtime_s = total,
             warm_start_s = ws,
             warm_start_frac = ws / total,
             ms_per_iter = 1000 * (total - ws) / n_iter,
             n_iter = n_iter,
             shd = LEWIS::shd(D_hat, ref),
             edges = n_edges_partial(D_hat),
             n_nG_hat = length(fit$nG_hat),
             n_unique_dags = fit$diagnostics$n_unique_dags,
             acceptance_rate = fit$diagnostics$acceptance_rate,
             engine = ENGINE,
             stringsAsFactors = FALSE)
}

# ---- Main loop -----------------------------------------------------------------
n_errors <- 0L
for (i in seq_len(nrow(cells))) {
  gtype <- cells$graph[i]
  p <- cells$p[i]
  n <- cells$n[i]
  inst <- build_instance(gtype, p)
  view <- instance_view(inst, RHO, EFFECT)
  ref <- true_rpdag_of(view)
  sched <- lewis_schedule(p)

  # The datasets of this setting that are not in the file yet.
  to_run <- setdiff(seq_len(N_REPS),
                  as.integer(sub(".*\\|", "",
                                 grep(sprintf("^%s\\|%d\\|%d\\|", gtype, p, n),
                                      done_key, value = TRUE))))
  cat(sprintf("[%s p = %d n = %d] %d iterations, burn-in %d, thinning %d, %d datasets to run\n",
              gtype, p, n, sched$n_iter, sched$burn_in, sched$thin, length(to_run)))
  if (length(to_run) == 0L) next

  # Run the datasets in parallel. An error is returned as its message.
  t_cell <- Sys.time()
  rows <- mclapply(to_run, function(r) {
    tryCatch(time_one(view, gtype, p, n, r, ref, sched$n_iter),
             error = function(e) sprintf("dataset %d: %s", r, conditionMessage(e)))
  }, mc.cores = MC_CORES, mc.preschedule = FALSE)

  # Report the failed datasets.
  failed <- !vapply(rows, is.data.frame, logical(1))
  for (k in which(failed)) {
    msg <- if (is.character(rows[[k]])) rows[[k]] else
      sprintf("dataset %d: no result (the worker process was killed)", to_run[k])
    message(sprintf("[%s p = %d n = %d] %s", gtype, p, n, msg))
  }
  n_errors <- n_errors + sum(failed)
  rows <- rows[!failed]

  # Save the rows after each setting.
  if (length(rows) > 0L) {
    all_rows[[length(all_rows) + 1L]] <- do.call(rbind, rows)
    write.csv(do.call(rbind, all_rows), OUT_CSV, row.names = FALSE)
    cat(sprintf("    %.1f minutes, median %.2f s per dataset\n",
                as.numeric(difftime(Sys.time(), t_cell, units = "mins")),
                median(do.call(rbind, rows)$runtime_s)))
  }
}
if (length(all_rows) == 0L) {
  message("No row was written.")
  quit(status = 1)
}
res <- do.call(rbind, all_rows)

# ---- The numbers of Appendix E.1 ---------------------------------------------------
cat("\nMedians over the datasets\n")
agg <- aggregate(cbind(runtime_s, ms_per_iter, warm_start_s, warm_start_frac) ~
                   graph + p + n, data = res, FUN = median)
agg <- agg[order(agg$p, agg$n, agg$graph), ]
cat(sprintf("%-3s %-4s %-5s %12s %14s %12s %10s\n",
            "gr", "p", "n", "total (s)", "ms/iteration", "TL (s)", "TL share"))
for (k in seq_len(nrow(agg))) {
  cat(sprintf("%-3s %-4d %-5d %12.2f %14.3f %12.2f %9.0f%%\n",
              agg$graph[k], agg$p[k], agg$n[k], agg$runtime_s[k],
              agg$ms_per_iter[k], agg$warm_start_s[k], 100 * agg$warm_start_frac[k]))
}

# ---- Check against the comparative study --------------------------------------------
# The chains are those of the comparative study, so on the datasets present
# in both files the SHD must agree exactly. The run times may differ.
bench <- file.path(.here, "results", sprintf("results_lewis_%s.csv", PRIOR))
n_mismatch <- 0L
if (file.exists(bench)) {
  b <- read.csv(bench, stringsAsFactors = FALSE)
  b <- b[b$experiment == "core" & b$rho == RHO & b$effect_scale == EFFECT,
         c("graph_type", "p", "n", "rep", "shd_rpdag")]
  names(b)[names(b) == "graph_type"] <- "graph"
  m <- merge(res[, c("graph", "p", "n", "rep", "shd")], b,
             by = c("graph", "p", "n", "rep"))
  n_mismatch <- sum(m$shd != m$shd_rpdag)
  cat(sprintf("\nSHD agrees with the comparative study on %d of %d datasets.\n",
              nrow(m) - n_mismatch, nrow(m)))
} else {
  cat("\nThe file", basename(bench), "is not in results/. The SHD was not checked.\n")
}

# Report the output file. Exit with status 1 if a dataset failed or an SHD
# differs.
cat("Written:", OUT_CSV, "\n")
if (n_errors > 0L || n_mismatch > 0L) {
  message(sprintf("NOT COMPLETE: %d datasets failed and %d SHD values differ from the comparative study.",
                  n_errors, n_mismatch))
  quit(status = 1)
}
