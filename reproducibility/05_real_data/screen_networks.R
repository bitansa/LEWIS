# =============================================================================
# screen_networks.R
#
# The screen of the 18 gene sets of Appendix I ("Choice of the gene sets").
# Each gene set is prepared once and analyzed in full by real_data_lewis.R
# (the reported chain, the ten stability runs and the four convergence
# chains) at the chain lengths RD_SCREEN_LENGTHS in turn: the schedule of
# lewis_schedule(), then 10^5, 10^6 and 10^7 iterations. The sets with 14
# genes stop at 10^6 and the sets with 8 genes at 10^7. A gene set stops at
# the first length at which it passes the four checks of real_data_config.R:
#   1. PSRF of the score below 1.05.
#   2. the highest scores of the four convergence chains within 0.5.
#   3. the posterior probabilities of every adjacency and orientation within
#      0.1 between the chains.
#   4. all ten stability runs return the same D_hat.
#
# The analyses of the screen are written to screen/ in this folder (or to
# the folder in the environment variable RD_SCREEN_DIR). The summary is
# written to results/tables/screen_summary.csv.
#
# MODES
#   Rscript screen_networks.R
#       Summary mode. Rebuilds the summary from the analyses that exist in
#       screen/. No chain is run. If screen/ holds no analysis (for example
#       in a fresh copy of the repository), the stored summary is kept and
#       a message explains how to run the screen.
#   SCREEN_RUN=all Rscript screen_networks.R
#       Runs the whole screen: all 18 gene sets, then the summary. This
#       takes several days on a laptop.
#   SCREEN_RUN=PRAD_ar8,LIHC_chol8 Rscript screen_networks.R
#       Runs the named gene sets only, then the summary.
#   SCREEN_FORCE=1 reruns analyses that exist already.
#
# NOTE ON THE STORED SUMMARY. For PRAD_chol8, PRAD_ifn8, LIHC_ppara8 and
# LIHC_myc8 the stored summary has rows up to 10^6 iterations. Their runs
# with 10^7 iterations were stopped once three gene sets had passed. A full
# run of the screen adds their rows at 10^7 iterations.
# =============================================================================

# ---- Folder of this script and settings -------------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
SCREEN_DIR <- Sys.getenv("RD_SCREEN_DIR", file.path(.here, "screen"))
Sys.setenv(RD_OUT_DIR = SCREEN_DIR)
source(file.path(.here, "real_data_config.R"))

# ---- Mode -------------------------------------------------------------------
RUN <- Sys.getenv("SCREEN_RUN", "")
FORCE <- identical(Sys.getenv("SCREEN_FORCE", "0"), "1")
NETWORKS <- character(0)
if (identical(RUN, "all")) {
  NETWORKS <- names(RD_NETWORKS)
} else if (nzchar(RUN)) {
  NETWORKS <- strsplit(RUN, ",", fixed = TRUE)[[1]]
}
unknown <- setdiff(NETWORKS, names(RD_NETWORKS))
if (length(unknown) > 0L) stop("Unknown gene sets in SCREEN_RUN: ", paste(unknown, collapse = ", "))

# Folders of the tables and the logs, and the file of the summary.
TAB_DIR <- RD_TAB_DIR
LOG_DIR <- file.path(SCREEN_DIR, "logs")
OUT_FILE <- file.path(.here, "results", "tables", "screen_summary.csv")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")

# The chain lengths of a gene set.
lengths_of <- function(net) {
  L <- RD_SCREEN_LENGTHS
  L[is.na(L) | L <= rd_screen_max_length(net)]
}

# Run one script of this folder with the given environment. Return TRUE if
# it succeeds. The output goes to a log file.
run_step <- function(script, env, log) {
  st <- system2(rscript, shQuote(file.path(.here, script)), env = env,
                stdout = file.path(LOG_DIR, log), stderr = file.path(LOG_DIR, log))
  identical(st, 0L)
}

# One row of the summary for a gene set and a chain length.
collect <- function(net, L) {
  tag <- rd_tag(net, L)
  rd <- function(stem) {
    read.csv(file.path(TAB_DIR, sprintf("%s_%s.csv", stem, tag)), stringsAsFactors = FALSE)
  }
  s <- rd("summary")
  nd <- rd("node_posterior")
  cv <- rd("convergence")
  sb <- rd("stability")
  # A gene is counted as Gaussian when it is declared non-Gaussian in at
  # most 5 percent of the draws.
  g <- nd$gene[nd$p_non_gaussian <= 0.05]
  data.frame(network = net, cohort = RD_NETWORKS[[net]]$cohort, length = rd_length_label(L),
             n = s$n, p = s$p,
             n_gaussian = length(g), gaussian_genes = paste(g, collapse = ";"),
             n_class_uncertain = sum(nd$p_non_gaussian > 0.05 & nd$p_non_gaussian < 0.95),
             psrf_log_score = cv$psrf_log_score[1],
             best_score_gap = max(cv$max_log_score) - min(cv$max_log_score),
             max_pair_diff = max(cv$max_diff_p_adjacent[1], cv$max_diff_p_directed[1]),
             n_reproduced = sum(sb$shd_to_reported == 0), n_stability = nrow(sb),
             mass_D_hat = s$mass_D_hat, eps_alpha = s$eps_alpha, n_edges = s$n_edges,
             agrees = s$agrees, against = s$against, stringsAsFactors = FALSE)
}

# The four checks.
passes <- function(r) {
  r$psrf_log_score < RD_CONV_MAX_PSRF & r$best_score_gap <= RD_CONV_MAX_GAP &
    r$max_pair_diff < RD_CONV_MAX_PDIFF & r$n_reproduced == r$n_stability
}

# ---- Run the gene sets of SCREEN_RUN ------------------------------------------
failed <- character(0)
for (net in NETWORKS) {
  cat(sprintf("%s (%s)\n", net, RD_NETWORKS[[net]]$label))
  env0 <- c(sprintf("RD_NETWORK=%s", net), sprintf("RD_OUT_DIR=%s", shQuote(SCREEN_DIR)))
  # Prepare the data once.
  if (FORCE || !file.exists(file.path(TAB_DIR, sprintf("expression_%s.csv", net)))) {
    if (!run_step("prepare_tcga.R", env0, sprintf("prepare_%s.log", net))) {
      message("The preparation failed for ", net, ". See ", file.path(LOG_DIR, sprintf("prepare_%s.log", net)))
      failed <- c(failed, net)
      next
    }
  }
  # Analyze at increasing chain lengths until the gene set passes.
  for (L in lengths_of(net)) {
    tag <- rd_tag(net, L)
    if (FORCE || !file.exists(file.path(TAB_DIR, sprintf("summary_%s.csv", tag)))) {
      t0 <- Sys.time()
      env <- c(env0, sprintf("RD_N_ITER=%s", if (is.na(L)) "" else L))
      if (!run_step("real_data_lewis.R", env, sprintf("%s.log", tag))) {
        message("The analysis failed for ", tag, ". See ", file.path(LOG_DIR, sprintf("%s.log", tag)))
        failed <- c(failed, tag)
        break
      }
      cat(sprintf("  %-22s %6.1f min\n", tag, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
    r <- collect(net, L)
    r$pass <- passes(r)
    cat(sprintf("  %-22s PSRF %.3f, score gap %.2f, difference %.3f, reproduced %d of %d: %s\n",
                tag, r$psrf_log_score, r$best_score_gap, r$max_pair_diff,
                r$n_reproduced, r$n_stability, if (r$pass) "pass" else "fail"))
    if (r$pass) break
  }
}

# ---- Summary --------------------------------------------------------------------
# The summary is rebuilt from every completed analysis in the screen folder,
# so that several screens run in parallel all leave the full summary.
rows <- list()
for (net in names(RD_NETWORKS)) for (L in lengths_of(net)) {
  if (!file.exists(file.path(TAB_DIR, sprintf("summary_%s.csv", rd_tag(net, L))))) break
  r <- collect(net, L)
  r$pass <- passes(r)
  rows[[length(rows) + 1L]] <- r
  if (r$pass) break
}
if (length(rows) == 0L) {
  cat("No analysis of the screen was found in ", SCREEN_DIR, ".\n",
      "The stored summary ", OUT_FILE, " is kept.\n",
      "Run the whole screen with:  SCREEN_RUN=all Rscript screen_networks.R\n", sep = "")
} else {
  out <- do.call(rbind, rows)
  write.csv(out, OUT_FILE, row.names = FALSE)
  cat("\n==== Summary of the screen ====\n")
  show <- out[, c("network", "length", "n", "n_gaussian", "psrf_log_score", "best_score_gap",
                  "max_pair_diff", "n_reproduced", "mass_D_hat", "eps_alpha", "agrees",
                  "against", "pass")]
  print(show, row.names = FALSE, digits = 3)
  cat(sprintf("\nWritten: %s\n", OUT_FILE))
}

# A failed step makes the script exit with an error.
if (length(failed) > 0L) {
  stop("These steps failed: ", paste(failed, collapse = ", "))
}
