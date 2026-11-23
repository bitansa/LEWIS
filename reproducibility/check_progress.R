# =============================================================================
# check_progress.R
#
# Reports the progress of the simulation studies from the result files on
# disk. It only reads files, so it can run at any time, also while jobs
# are still writing.
#
#     Rscript reproducibility/check_progress.R              # all studies
#     Rscript reproducibility/check_progress.R comparison   # one study
#
# Studies: motivation (01_motivation), comparison (02_comparison),
# calibration (03_calibration), laplace (06_laplace_vs_exact), and health
# (a check of the comparative study and of the study of Figure 1).
#
# For each study it reports the rows written against the rows expected and
# the number of rows with skipped = TRUE. A method wrapper that fails on a
# dataset prints the error and writes a row with skipped = TRUE and missing
# metrics (see common/method_wrappers.R). So a failing method shows up as a
# large count of skipped rows.
#
# The section "health" adds a check that does not rely on that flag. It
# reports, for each method, the median number of estimated edges, the
# median run time and the share of empty estimates. A method that returns
# the empty graph on most datasets is marked for inspection.
# =============================================================================

# ---- Folder of this script ----------------------------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})

# The studies to report. All of them if none is given.
sections <- commandArgs(trailingOnly = TRUE)
if (length(sections) == 0L) {
  sections <- c("motivation", "comparison", "calibration", "laplace", "health")
}

# ---- Helpers --------------------------------------------------------------------
hdr <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")
pct <- function(a, b) if (b > 0) sprintf("%5.1f%%", 100 * a / b) else "   NA "

# Read a result file that another process may be writing. A file that
# cannot be read is reported, not hidden.
read_now <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) {
    cat(sprintf("Could not read %s: %s\n", basename(path), conditionMessage(e)))
    NULL
  })
}

# =============================================================================
# 1. Figure 1 (01_motivation)
# =============================================================================
if ("motivation" %in% sections) {
  hdr("FIGURE 1 (01_motivation)")
  d <- read_now(file.path(.here, "01_motivation", "results", "motivation_results.csv"))
  if (is.null(d) || nrow(d) == 0L) {
    cat("No rows yet.\n")
  } else {
    # 25 regimes, 3 sample sizes, 100 datasets and 10 methods. The nine
    # methods in R are run by motivation_experiment.R. Pairwise-LiNGAM is
    # added by the import of the MATLAB results.
    n_cells <- 25L * 3L
    expected <- n_cells * 100L * 10L
    cells <- unique(d[, c("regime", "n")])
    cat(sprintf("Rows: %d of %d (%s). Settings: %d of %d. Methods: %d of 10.\n",
                nrow(d), expected, pct(nrow(d), expected), nrow(cells), n_cells,
                length(unique(d$method))))
    cat(sprintf("Rows with an SHD: %d. Rows with skipped = TRUE: %d.\n\n",
                sum(!is.na(d$shd_rpdag)), sum(as.logical(d$skipped))))
    tb <- do.call(rbind, lapply(split(d, d$method), function(z) {
      data.frame(method = z$method[1], rows = nrow(z), skipped = sum(as.logical(z$skipped)),
                 median_shd = round(median(z$shd_rpdag, na.rm = TRUE), 1))
    }))
    print(tb, row.names = FALSE)
  }
}

# =============================================================================
# 2. Comparative study (02_comparison)
# =============================================================================
if ("comparison" %in% sections) {
  hdr("COMPARATIVE STUDY (02_comparison)")
  source(file.path(.here, "02_comparison", "comparison_cells.R"))
  cells <- comparison_cells()
  expected <- nrow(cells) * N_REPS
  cat(sprintf("%d settings x %d datasets = %d rows expected per method.\n\n",
              nrow(cells), N_REPS, expected))
  files <- sort(list.files(file.path(.here, "02_comparison", "results"),
                           pattern = "^results_.*\\.csv$", full.names = TRUE))
  if (length(files) == 0L) {
    cat("No result files yet.\n")
  } else {
    cat(sprintf("%-20s %8s %9s %8s %10s %10s %9s\n", "Method", "Rows", "Settings",
                "Done", "SHD rows", "SID rows", "Skipped"))
    for (f in files) {
      d <- read_now(f)
      if (is.null(d) || nrow(d) == 0L) next
      m <- sub("^results_(.*)\\.csv$", "\\1", basename(f))
      n_set <- nrow(unique(d[, c("experiment", "graph_type", "p", "n", "rho",
                                 "effect_scale")]))
      cat(sprintf("%-20s %8d %9d %8s %10d %10d %9d\n", m, nrow(d), n_set,
                  pct(nrow(d), expected), sum(!is.na(d$shd_rpdag)),
                  sum(!is.na(d$sid_upper)), sum(as.logical(d$skipped))))
    }
    cat("\nSettings: settings with at least one row. Done: rows over rows expected.\n",
        "A large count of skipped rows means that the method is failing.\n", sep = "")
  }
}

# =============================================================================
# 3. Calibration study (03_calibration)
# =============================================================================
if ("calibration" %in% sections) {
  hdr("CALIBRATION STUDY (03_calibration)")
  # 84 settings, 10 values of the constant and 100 datasets.
  expected <- 84L * 10L * 100L
  for (nm in c("d_sweep", "alpha_sweep")) {
    par_col <- if (nm == "d_sweep") "d" else "alpha"
    # The result file of the sweep.
    fs <- list.files(file.path(.here, "03_calibration", "results"),
                     pattern = sprintf("^%s_results\\.csv$", nm), full.names = TRUE)
    d <- do.call(rbind, lapply(fs, function(f) {
      z <- read_now(f)
      if (is.null(z)) NULL else z[, c("scenario", par_col, "rep", "shd_rpdag")]
    }))
    if (is.null(d) || nrow(d) == 0L) {
      cat(sprintf("%-12s: no rows yet\n", nm))
      next
    }
    # Count each (setting, value, dataset) once.
    d <- d[!duplicated(d[, c("scenario", par_col, "rep")]), ]
    cat(sprintf("%-12s: %6d of %d distinct rows (%s), settings %d of 84, values %d of 10, median SHD %.1f\n",
                nm, nrow(d), expected, pct(nrow(d), expected),
                length(unique(d$scenario)), length(unique(d[[par_col]])),
                median(d$shd_rpdag, na.rm = TRUE)))
  }
}

# =============================================================================
# 4. Accuracy of the Laplace approximation (06_laplace_vs_exact)
# =============================================================================
if ("laplace" %in% sections) {
  hdr("LAPLACE APPROXIMATION AND EXACT MARGINAL LIKELIHOOD (06_laplace_vs_exact)")
  f <- list.files(file.path(.here, "06_laplace_vs_exact", "results"),
                  pattern = "^exact_vs_laplace\\.csv$", recursive = TRUE, full.names = TRUE)
  d <- if (length(f) > 0L) read_now(f[1]) else NULL
  if (is.null(d) || nrow(d) == 0L) {
    cat("No rows yet.\n")
  } else {
    # 13 pairs (p, n), 100 datasets, two chains and two priors.
    expected <- 13L * 100L * 4L
    cat(sprintf("Rows: %d of %d (%s). Pairs (p, n): %d of 13.\n", nrow(d), expected,
                pct(nrow(d), expected), nrow(unique(d[, c("p", "n")]))))
  }
}

# =============================================================================
# 5. Health check of the comparative study and of Figure 1
# =============================================================================
# A method whose typical estimate is the empty graph, or whose run time is
# close to zero, is marked for inspection. A healthy method estimates about
# p edges on ER graphs and about 2p edges on SF graphs, and its run time
# grows with p.
if ("health" %in% sections) {
  hdr("HEALTH CHECK: empty estimates and very short run times")

  report <- function(label, d) {
    if (is.null(d) || nrow(d) == 0L) return(invisible(NULL))
    cat(sprintf("\n%s\n", label))
    cat(sprintf("  %-20s %10s %11s %14s %10s\n", "Method", "Med. edges",
                "Med. time s", "Empty rows", "Verdict"))
    for (m in sort(unique(d$method))) {
      z <- d[d$method == m & !as.logical(d$skipped), ]
      if (nrow(z) == 0L) next
      me <- median(z$n_edges_est)
      mt <- median(z$runtime_s, na.rm = TRUE)
      n0 <- sum(z$n_edges_est == 0)
      bad <- isTRUE(me == 0) || (n0 > 0.5 * nrow(z))
      cat(sprintf("  %-20s %10.1f %11.3f %6d (%3.0f%%) %10s\n", m, me, mt, n0,
                  100 * n0 / nrow(z), if (bad) "CHECK" else "ok"))
    }
  }

  # Comparative study.
  cfs <- sort(list.files(file.path(.here, "02_comparison", "results"),
                         pattern = "^results_.*\\.csv$", full.names = TRUE))
  if (length(cfs) > 0L) {
    cmp <- do.call(rbind, lapply(cfs, read_now))
    report("COMPARATIVE STUDY", cmp)
    ok <- !as.logical(cmp$skipped)
    cat("\n  Median number of estimated edges by p (ER: about p, SF: about 2p):\n")
    print(round(tapply(cmp$n_edges_est[ok], list(cmp$method[ok], cmp$p[ok]), median), 1))
  }

  # Figure 1 (the true graph has 18 edges).
  mot <- read_now(file.path(.here, "01_motivation", "results", "motivation_results.csv"))
  report("FIGURE 1 (the true graph has 18 edges)", mot)

  cat("\n  CHECK means that the typical estimate of the method is the empty graph.\n",
      "  Inspect the error log of its job, for example\n",
      "     grep -i 'failed' reproducibility/slurm/logs/cmp_<method>_*.err\n", sep = "")
}

cat("\n")
