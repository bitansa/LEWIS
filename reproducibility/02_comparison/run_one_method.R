# =============================================================================
# run_one_method.R
#
# Runs one method on all 54 settings and 100 datasets of the comparative
# study (Section 7.2), and writes one row per dataset to
#
#     results/results_<method>.csv
#
# Usage, from any folder, for example from the root of the repository:
#
#     MC_CORES=100 Rscript reproducibility/02_comparison/run_one_method.R <method>
#
# with <method> one of lewis_complexity, lewis_uniform, pc, ges, lingam,
# direct_lingam, rcd, tl_zhao, notears. Pairwise-LiNGAM is run by the
# scripts of pairwise_lingam/.
#
# DESIGN
#   * One file per method, so that the methods run as independent jobs.
#   * The datasets are generated again inside this process by
#     draw_dataset() with the regime tag "mixture". The seed depends only on
#     the setting and the dataset number, so every method sees the same
#     datasets. A method that uses random numbers continues the random
#     stream that draw_dataset() started, so every row can be reproduced on
#     its own.
#   * Resume: the rows already in the file without failure are skipped, so
#     a job that is submitted again continues where it stopped. A dataset
#     whose row has skipped = TRUE is run again, and the new row replaces
#     the old one.
#   * The datasets of one setting run in parallel with parallel::mclapply()
#     on MC_CORES workers (environment variable, default: all cores but one).
#   * There is no time limit on the methods or on the SID.
#
# FAILURES. If a method fails on a dataset, the row is written with
# skipped = TRUE and the message is printed. If the evaluation of a dataset
# fails, no row is written and the message is printed. At the end the script
# checks the file and exits with status 1 if a setting has fewer than 100
# distinct datasets or if a row has skipped = TRUE. Fix the cause and run
# the script again. It then runs only the missing and failed datasets.
#
# Columns of the output:
#     experiment, graph_type, p, n, rho, effect_scale, rep, method,
#     shd_rpdag, sid_lower, sid_upper, n_edges_est, runtime_s, skipped
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

# ---- The method -------------------------------------------------------------
cli <- commandArgs(trailingOnly = TRUE)
METHOD_NAME <- if (length(cli) >= 1L) cli[[1L]] else ""
VALID_METHODS <- c("lewis_complexity", "lewis_uniform", "pc", "ges", "lingam",
                   "direct_lingam", "rcd", "tl_zhao", "notears")
if (!METHOD_NAME %in% VALID_METHODS) {
  stop("Give the method as the first argument. Valid methods: ",
       paste(VALID_METHODS, collapse = ", "), ".")
}

# ---- Shared code --------------------------------------------------------------
source(file.path(COMMON_DIR, "scenario_instances.R"))
source(file.path(COMMON_DIR, "method_wrappers.R"))
source(file.path(COMMON_DIR, "scoring.R"))
source(file.path(.here, "comparison_cells.R"))

# Stop at once if the software of the method is missing.
require_method(METHOD_NAME)

# ---- Number of workers --------------------------------------------------------
# The environment variable MC_CORES, or all cores but one.
MC_CORES <- {
  env <- suppressWarnings(as.integer(Sys.getenv("MC_CORES", "")))
  if (is.finite(env) && env >= 1L) env else max(1L, detectCores() - 1L)
}

# ---- Settings and output file -------------------------------------------------
# Create results/ if needed.
cells <- comparison_cells()
csv_dir <- file.path(.here, "results")
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
csv_out <- file.path(csv_dir, paste0("results_", METHOD_NAME, ".csv"))
cat(sprintf("[%s] %d settings, %d datasets each, %d workers\n",
            METHOD_NAME, nrow(cells), N_REPS, MC_CORES))
cat(sprintf("[%s] output: %s\n", METHOD_NAME, csv_out))

# ---- Keys of the rows already written ------------------------------------------
# Only rows without failure count as done.
row_key <- function(experiment, graph_type, p, n, rho, effect_scale, rep) {
  sprintf("%s|r%d", cell_key(experiment, graph_type, p, n, rho, effect_scale),
          as.integer(rep))
}
done_keys <- character(0)
if (file.exists(csv_out)) {
  prev <- read.csv(csv_out, stringsAsFactors = FALSE)
  if (nrow(prev) > 0L) {
    ok <- !as.logical(prev$skipped)
    done_keys <- row_key(prev$experiment, prev$graph_type, prev$p, prev$n,
                         prev$rho, prev$effect_scale, prev$rep)[ok]
    cat(sprintf("[%s] resume: %d rows are already in the file, %d without failure\n",
                METHOD_NAME, nrow(prev), sum(ok)))
  }
}

# ---- One dataset --------------------------------------------------------------
run_one_rep <- function(cell, r) {
  # Generate dataset r of this setting.
  inst <- build_instance(cell$graph_type, cell$p)
  view <- instance_view(inst, rho = cell$rho, effect_scale = cell$effect_scale)
  X <- draw_dataset(view, n = cell$n, rep = r, regime_tag = "mixture")
  # Run the method and evaluate its estimate.
  out <- run_method(METHOD_NAME, X)
  met <- score_estimate(METHOD_NAME, out, view)
  cbind(data.frame(experiment = cell$experiment, graph_type = cell$graph_type,
                   p = cell$p, n = cell$n, rho = cell$rho,
                   effect_scale = cell$effect_scale, rep = r,
                   method = METHOD_NAME, stringsAsFactors = FALSE),
        met)
}

# ---- Main loop: settings in turn, datasets in parallel --------------------------
t_start <- Sys.time()
n_errors <- 0L
for (ci in seq_len(nrow(cells))) {
  cell <- cells[ci, ]
  # The datasets of this setting that are not in the file yet.
  to_run <- Filter(function(r) {
    !(row_key(cell$experiment, cell$graph_type, cell$p, cell$n, cell$rho,
              cell$effect_scale, r) %in% done_keys)
  }, seq_len(N_REPS))
  if (length(to_run) == 0L) next

  cat(sprintf("[%s] setting %d of %d: %s %s p = %d n = %d rho = %g s = %g (%d datasets)\n",
              METHOD_NAME, ci, nrow(cells), cell$experiment, cell$graph_type,
              cell$p, cell$n, cell$rho, cell$effect_scale, length(to_run)))

  # Run the datasets. An error is returned as its message.
  t0 <- Sys.time()
  rows <- mclapply(to_run, function(r) {
    tryCatch(run_one_rep(cell, r),
             error = function(e) sprintf("dataset %d: %s", r, conditionMessage(e)))
  }, mc.cores = MC_CORES, mc.preschedule = FALSE)

  # Report the datasets that gave no row. A NULL result means that the
  # worker process was killed.
  failed <- !vapply(rows, is.data.frame, logical(1))
  for (k in which(failed)) {
    msg <- if (is.character(rows[[k]])) rows[[k]] else
      sprintf("dataset %d: no result (the worker process was killed)", to_run[k])
    message(sprintf("[%s] setting %d, %s", METHOD_NAME, ci, msg))
  }
  n_errors <- n_errors + sum(failed)
  rows <- rows[!failed]

  # Append the new rows to the file. Failed rows of an earlier run that the
  # new rows replace are removed first.
  if (length(rows) > 0L) {
    block <- do.call(rbind, rows)
    bkey <- row_key(block$experiment, block$graph_type, block$p, block$n,
                    block$rho, block$effect_scale, block$rep)
    if (file.exists(csv_out)) {
      prev <- read.csv(csv_out, stringsAsFactors = FALSE)
      pkey <- row_key(prev$experiment, prev$graph_type, prev$p, prev$n,
                      prev$rho, prev$effect_scale, prev$rep)
      drop <- as.logical(prev$skipped) & pkey %in% bkey
      if (any(drop)) {
        write.csv(prev[!drop, , drop = FALSE], csv_out, row.names = FALSE)
      }
    }
    new_file <- !file.exists(csv_out)
    write.table(block, csv_out, sep = ",", row.names = FALSE,
                col.names = new_file, append = !new_file)
    done_keys <- c(done_keys, bkey[!block$skipped])
  }
  cat(sprintf("[%s]   %.1f minutes, %d rows written, %d datasets failed\n",
              METHOD_NAME, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              length(rows), sum(failed)))
}

# ---- Final check of the whole file ---------------------------------------------
if (!file.exists(csv_out)) {
  message(sprintf("[%s] No row was written. Check the messages above.", METHOD_NAME))
  quit(status = 1)
}
res <- read.csv(csv_out, stringsAsFactors = FALSE)
res_cell <- cell_key(res$experiment, res$graph_type, res$p, res$n, res$rho,
                     res$effect_scale)
grid_cell <- cell_key(cells$experiment, cells$graph_type, cells$p, cells$n,
                      cells$rho, cells$effect_scale)
# Number of distinct datasets of each setting.
n_distinct <- vapply(grid_cell, function(k) length(unique(res$rep[res_cell == k])),
                     integer(1))
short <- grid_cell[n_distinct < N_REPS]
n_skipped <- sum(as.logical(res$skipped))
n_dup <- sum(duplicated(paste(res_cell, res$rep)))

# Report the check. Exit with status 1 if the file is not complete.
cat(sprintf("[%s] %.1f hours. %d rows, %d settings with fewer than %d datasets, %d rows with skipped = TRUE, %d duplicated rows.\n",
            METHOD_NAME, as.numeric(difftime(Sys.time(), t_start, units = "hours")),
            nrow(res), length(short), N_REPS, n_skipped, n_dup))
if (length(short) > 0L || n_skipped > 0L || n_dup > 0L || n_errors > 0L) {
  for (k in head(short, 20L)) message("  incomplete setting: ", k)
  message(sprintf("[%s] NOT COMPLETE. Check the messages above, fix the cause ",
                  METHOD_NAME),
          "and run the script again. It runs only the missing and failed datasets.")
  quit(status = 1)
}
cat(sprintf("[%s] Complete: all %d settings have %d datasets.\n",
            METHOD_NAME, nrow(cells), N_REPS))
