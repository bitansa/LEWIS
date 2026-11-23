# =============================================================================
# scoring.R
#
# Evaluation of one estimate in the comparative study (Section 7.2 and
# Appendix G.1).
#
# REFERENCE. The reference of a setting is the true RPDAG, the RPDAG of the
# true DAG and its non-Gaussian nodes, LEWIS::rpdag(view$dag, view$nG_star).
# It is the same for all datasets of a setting, so it is computed once per
# setting.
#
# METRICS of one row:
#     shd_rpdag     SHD (LEWIS::shd) between the graph of the method (see
#                   est_partial_dag() in method_wrappers.R) and the true
#                   RPDAG
#     sid_lower,    lower and upper SID bounds of that graph against the
#     sid_upper     true DAG, without any time limit
#     n_edges_est   number of edges of that graph (an undirected edge
#                   counts once)
#     runtime_s     run time of the method in seconds
#     skipped       TRUE when the method failed on this dataset. The metric
#                   columns are then NA, and the calling script reports the
#                   row and exits with a non-zero status.
#
# An error while scoring is not caught. It stops the evaluation of the
# dataset, and the calling script reports it.
#
# The calling script must source method_wrappers.R first.
# =============================================================================

.true_rpdag_memo <- new.env(parent = emptyenv())

# The true RPDAG of a setting, cached by instance and rho. The effect-size
# scale changes neither the graph nor the non-Gaussian nodes.
true_rpdag_of <- function(view) {
  key <- sprintf("%s|rho%s", view$key, format(view$rho, trim = TRUE))
  hit <- .true_rpdag_memo[[key]]
  if (!is.null(hit)) return(hit)
  ref <- unname(LEWIS::rpdag(view$dag, view$nG_star))
  .true_rpdag_memo[[key]] <- ref
  ref
}

# Evaluate the output `out` of run_method() for the setting `view`.
# Returns a one-row data frame with the metric columns.
score_estimate <- function(method_name, out, view) {
  # A failed method gives a flagged row with missing metrics.
  if (isTRUE(out$skipped)) {
    return(data.frame(shd_rpdag = NA_real_, sid_lower = NA_real_,
                      sid_upper = NA_real_, n_edges_est = NA_integer_,
                      runtime_s = out$runtime, skipped = TRUE))
  }
  # The graph of the method and the reference.
  ref <- true_rpdag_of(view)
  est <- est_partial_dag(method_name, out$dag)
  # SHD from the true RPDAG and SID bounds against the true DAG.
  shd_val <- LEWIS::shd(est, ref)
  sids <- sid_bounds(est, view$dag)
  data.frame(shd_rpdag = shd_val,
             sid_lower = sids$sid_lower,
             sid_upper = sids$sid_upper,
             n_edges_est = n_edges_partial(est),
             runtime_s = out$runtime,
             skipped = FALSE)
}
