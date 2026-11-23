# =============================================================================
# run_reported_analyses.R
#
# Runs the two reported analyses of Section 7.4 (RD_REPORTED in
# real_data_config.R) from the raw data to the tables and figures:
#   1. prepare_tcga.R for each reported network.
#   2. real_data_lewis.R for each, with 10^7 iterations per chain: the
#      reported chain, the ten stability runs and the four convergence
#      chains.
#   3. plot_real_data_figures.R, when the tables of both networks exist.
# The tables go to results/tables/ and replace the stored tables. The
# figures go to <LEWIS_FIG_DIR>/real_data/TCGA/. The output of each step is
# written to a log file in results/logs/. A failed step stops the script.
# The raw data must be in data/raw/ (run download_tcga.R first).
#
# Run:  Rscript run_reported_analyses.R
#       The two networks can also run in parallel in two terminals:
#       REPORTED_ONLY=PRAD_ar8 Rscript run_reported_analyses.R
#       REPORTED_ONLY=LIHC_chol8 Rscript run_reported_analyses.R
#       Each run redraws the figures from the tables present in
#       results/tables/. The run that finishes last draws them from the new
#       tables of both networks.
# =============================================================================

# ---- Folder of this script and settings -------------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
OUT_DIR <- file.path(.here, "results")
Sys.setenv(RD_OUT_DIR = OUT_DIR)
source(file.path(.here, "real_data_config.R"))

# The folder of the logs and the Rscript of this R installation.
LOG_DIR <- file.path(OUT_DIR, "logs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")

# Run one script of this folder with the given environment. Stop if it fails.
step <- function(script, env, log) {
  lg <- file.path(LOG_DIR, log)
  t0 <- Sys.time()
  st <- system2(rscript, shQuote(file.path(.here, script)), env = env,
                stdout = lg, stderr = lg)
  if (!identical(st, 0L)) stop(script, " failed. See the log ", lg)
  cat(sprintf("  %-26s %6.1f min  (log: results/logs/%s)\n", script,
              as.numeric(difftime(Sys.time(), t0, units = "mins")), log))
}

# ---- Steps 1 and 2 for each reported network --------------------------------
ONLY <- Sys.getenv("REPORTED_ONLY", "")
if (nzchar(ONLY) && !ONLY %in% RD_REPORTED$network) {
  stop("REPORTED_ONLY must be one of: ", paste(RD_REPORTED$network, collapse = ", "))
}
for (k in seq_len(nrow(RD_REPORTED))) {
  net <- RD_REPORTED$network[k]
  if (nzchar(ONLY) && net != ONLY) next
  cat(sprintf("%s (%s)\n", RD_NETWORKS[[net]]$label, net))
  env <- c(sprintf("RD_NETWORK=%s", net),
           sprintf("RD_N_ITER=%d", RD_REPORTED$n_iter[k]),
           sprintf("RD_OUT_DIR=%s", shQuote(OUT_DIR)))
  step("prepare_tcga.R", env, sprintf("prepare_%s.log", net))
  step("real_data_lewis.R", env, sprintf("analysis_%s.log", net))
}

# ---- Step 3: the figures, once every reported network is complete ----------
tags <- vapply(seq_len(nrow(RD_REPORTED)), function(k) {
  rd_tag(RD_REPORTED$network[k], RD_REPORTED$n_iter[k])
}, character(1))
if (all(file.exists(file.path(RD_TAB_DIR, sprintf("summary_%s.csv", tags))))) {
  step("plot_real_data_figures.R", character(0), "plot_real_data_figures.log")
} else {
  cat("The figures are not drawn yet, because not every reported network is complete.\n")
}
cat("Done.\n")
