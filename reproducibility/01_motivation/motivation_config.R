# =============================================================================
# motivation_config.R
#
# The settings of the study of Figure 1 (design in Appendix G.1, "Design of
# Figure 1"). It is sourced by motivation_experiment.R, plot_motivation.R
# and the two Pairwise-LiNGAM scripts, so that all of them use the same
# graph, regimes and datasets.
#
# THE GRAPH is the 12-node DAG of Figure 7(a) of the paper
# (common/paper12_instance.R), with coefficients and error scales drawn
# once from a fixed seed.
#
# REGIMES. Figure 1 has two rows. In both rows, the first k nodes of the
# order NG_ORDER = (1, 2, 3, 5, 10, 4, 6, 7, 8, 9, 11, 12) have Student's
# t5 errors, for k = 12, 11, ..., 0. The rows differ in the errors of the
# other 12 - k nodes:
#
#     top row     Gaussian errors. The true RPDAG goes from the CPDAG at
#                 k = 0 to the DAG at k = 12.
#     bottom row  Student's t15 errors. Every node is non-Gaussian, so the
#                 true RPDAG is the DAG at every k.
#
# The two rows share the regime with k = 12 (all_t5), which is run once.
# This gives 25 regimes:
#
#     all_t5                    k = 12 in both rows
#     count_t5g_k01 .. k11      t5 on k nodes, Gaussian on the others
#     gaussian                  k = 0 of the top row
#     count_t5t15_k01 .. k11    t5 on k nodes, t15 on the others
#     t15                       k = 0 of the bottom row
#
# SEEDS. Dataset r of regime R and sample size n has the
# seed stable_seed("errors", "instance|paper12", "rho0", "s1", "n<n>",
# "motivation12_<R>", "rep<r>"). It depends only on the regime name, n and
# r, not on the order or the number of the regimes. The regime names are
# part of the seeds and must never be renamed.
#
# The fields `family`, `k_ng` and `df` are stored in the results. Their
# values are those of the runs of the paper: `family` is "shared" for
# gaussian and all_t5, "t5_gauss" and "t5_t15" for the two rows, and "df"
# for t15. `k_ng` is the number of nodes with non-Gaussian errors, so it is
# 12 for t15. The position k of a regime in Figure 1 is taken from
# INTRO_ROWS below, where t15 has k = 0.
#
# Requires common/scenario_instances.R and common/paper12_instance.R.
# =============================================================================

stopifnot(exists("VIEW12"))

# Sample sizes (the columns of Figure 1) and number of datasets per setting.
N_VALUES <- c(100L, 250L, 500L)
N_REPS <- 100L

# Error distribution of the non-Gaussian nodes in both rows, and of the
# other nodes in the bottom row.
NG_SHAPE_T5 <- "t5"
PROXY_G_SHAPE <- "t15"

# -----------------------------------------------------------------------------
# make_regimes: the 25 regimes of Figure 1.
# -----------------------------------------------------------------------------
# `nG` is the true set of non-Gaussian nodes. It defines the true RPDAG of
# the regime. In the bottom row it is the set of all 12 nodes at every k,
# because the t15 errors are non-Gaussian.
make_regimes <- function() {
  reg <- list()
  add <- function(...) reg[[length(reg) + 1L]] <<- list(...)

  # All errors Gaussian (k = 0 of the top row) and all errors t5 (k = 12 of
  # both rows).
  add(regime = "gaussian", family = "shared", k_ng = 0L, df = Inf,
      shapes = rep("gaussian", P12), nG = integer(0))
  add(regime = "all_t5", family = "shared", k_ng = 12L, df = 5,
      shapes = rep(NG_SHAPE_T5, P12), nG = seq_len(P12))

  # Top row: t5 on the first k nodes of NG_ORDER, Gaussian on the others.
  for (k in 1:11) {
    ng <- sort(NG_ORDER[seq_len(k)])
    sh <- rep("gaussian", P12)
    sh[ng] <- NG_SHAPE_T5
    add(regime = sprintf("count_t5g_k%02d", k), family = "t5_gauss",
        k_ng = k, df = NA_real_, shapes = sh, nG = ng)
  }

  # Bottom row: t5 on the first k nodes of NG_ORDER, t15 on the others.
  for (k in 1:11) {
    ng <- sort(NG_ORDER[seq_len(k)])
    sh <- rep(PROXY_G_SHAPE, P12)
    sh[ng] <- NG_SHAPE_T5
    add(regime = sprintf("count_t5t15_k%02d", k), family = "t5_t15",
        k_ng = k, df = NA_real_, shapes = sh, nG = seq_len(P12))
  }

  # The end of the bottom row: t15 on all nodes.
  add(regime = "t15", family = "df", k_ng = 12L, df = 15,
      shapes = rep(PROXY_G_SHAPE, P12), nG = seq_len(P12))
  reg
}
REGIMES <- make_regimes()
REGIME_NAMES <- vapply(REGIMES, function(z) z$regime, character(1))
stopifnot(length(REGIMES) == 25L, !anyDuplicated(REGIME_NAMES))

# The two rows of Figure 1 as tables of (k, regime). The plot script reads
# these tables.
INTRO_ROWS <- list(
  t5_gauss = data.frame(
    k_ng = 12:0,
    regime = c("all_t5", sprintf("count_t5g_k%02d", 11:1), "gaussian"),
    stringsAsFactors = FALSE),
  t5_t15 = data.frame(
    k_ng = 12:0,
    regime = c("all_t5", sprintf("count_t5t15_k%02d", 11:1), "t15"),
    stringsAsFactors = FALSE)
)
INTRO_ROW_LABELS <- c(
  t5_gauss = "Non-Gaussian nodes t5, remaining nodes Gaussian",
  t5_t15 = "Non-Gaussian nodes t5, remaining nodes t15")
stopifnot(setequal(unlist(lapply(INTRO_ROWS, function(z) z$regime)), REGIME_NAMES))

# -----------------------------------------------------------------------------
# draw_motivation_dataset: one dataset of a regime.
# -----------------------------------------------------------------------------
# The regime tag "motivation12_<regime>" is part of the seed.
draw_motivation_dataset <- function(rg, n, rep, return_errors = FALSE) {
  draw_dataset(VIEW12, n = n, rep = rep,
               regime_tag = paste0("motivation12_", rg$regime),
               shapes = rg$shapes, return_errors = return_errors)
}
