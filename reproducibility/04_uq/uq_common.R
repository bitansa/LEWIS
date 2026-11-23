# =============================================================================
# uq_common.R
#
# The setting of the uncertainty quantification study (Section 7.3 and
# Appendix H), and the functions that run LEWIS (complexity) on one dataset
# and write its tables. It is sourced by uq_decision_rule_demo.R and by
# choose_seed_for_demo.R, so that the two scripts run exactly the same
# computation for a given seed.
#
# SETTING. The ER graph with p = 10 of the comparative study (Section 7.2),
# which has 9 edges, with rho = 0.5, s = 1 and n = 500. The graph, its
# coefficients, its error scales and its non-Gaussian nodes are those of the
# comparative study. Only the errors depend on the seed. The dataset of seed
# s is draw_dataset(view, n = 500, rep = s, regime_tag = "uq_demo").
#
# SEEDING. The computation of the paper draws the dataset, calls
# set.seed(s), computes the TL start once, and then runs the chain, which
# computes the TL start again. The TL start may draw random numbers, for
# example when it removes edges to break a cycle. So run_uq() makes the same
# first TL computation before it calls lewis() without a seed. The chain
# then uses the same random numbers as in the paper.
#
# FUNCTIONS
#     uq_dataset(seed)          the dataset of a seed and the truth
#     uq_tl_start(X)            the first TL computation
#     run_uq(seed)              LEWIS (complexity) on the dataset of a seed,
#                               and the tables of the study
#     write_uq_tables(res, dir) writes the tables of run_uq() to a folder
#
# The calling script must define COMMON_DIR before it sources this file.
# Dependencies: LEWIS, igraph.
# =============================================================================

stopifnot(exists("COMMON_DIR"))
source(file.path(COMMON_DIR, "settings.R"))
source(file.path(COMMON_DIR, "scenario_instances.R"))
suppressPackageStartupMessages(library(LEWIS))

# ---- The setting ------------------------------------------------------------
UQ_P <- 10L
UQ_N <- 500L
UQ_RHO <- 0.5
UQ_S <- 1.0
UQ_REGIME_TAG <- "uq_demo"
UQ_ALPHA_BALL <- 0.05
# Number of graphs shown in the figures and in posterior_uncertainty_table.csv.
UQ_TOP_K <- 15L

# ---- The dataset of one seed -------------------------------------------------
# Returns the data matrix with the column names "1", ..., "10", so that every
# table names the nodes by their numbers, the true DAG and the true
# non-Gaussian nodes. The seed is the dataset number (argument rep) of
# draw_dataset().
uq_dataset <- function(seed) {
  stopifnot(length(seed) == 1L, seed >= 1, seed == round(seed))
  # The fixed instance of the setting and its view.
  inst <- build_instance("ER", UQ_P)
  view <- instance_view(inst, rho = UQ_RHO, effect_scale = UQ_S)
  # The errors of this seed. draw_dataset() sets its own seed.
  X <- draw_dataset(view, n = UQ_N, rep = as.integer(seed),
                    regime_tag = UQ_REGIME_TAG)
  colnames(X) <- as.character(seq_len(UQ_P))
  list(X = X, dag = view$dag, nG_star = view$nG_star)
}

# ---- The first TL computation of the runs of the paper -----------------------
# The same steps as the start of lewis(): the TL estimate, then removal of
# random edges while the graph has a cycle. Returns the start graph.
uq_tl_start <- function(X) {
  p <- ncol(X)
  # Run TL. An error is reported, and the start is then the empty graph, as
  # in lewis().
  A <- tryCatch({
    Bmat <- as.matrix(LEWIS::tl_lingam(unname(X))$B)
    Bmat[!is.finite(Bmat)] <- 0
    out <- matrix(0L, p, p)
    out[abs(Bmat) > 1e-8] <- 1L
    diag(out) <- 0L
    out
  }, error = function(e) {
    message("The TL method failed: ", conditionMessage(e))
    matrix(0L, p, p)
  })
  # Remove random edges while the graph has a cycle.
  is_acyclic <- function(M) igraph::is_dag(igraph::graph_from_adjacency_matrix(M, mode = "directed"))
  while (!is_acyclic(A)) {
    ij <- which(A == 1L, arr.ind = TRUE)
    drop <- ij[sample.int(nrow(ij), 1L), ]
    A[drop[1L], drop[2L]] <- 0L
  }
  A
}

# ---- Helpers for the edit lists ----------------------------------------------
# Number of entries of an edit list such as "5->1 | 4 - 5".
.n_entries <- function(s) {
  ifelse(is.na(s) | !nzchar(s), 0L, lengths(strsplit(s, " | ", fixed = TRUE)))
}
# Split a list of nodes such as "1, 2, 3" into a character vector.
.split_nodes <- function(s) {
  if (is.na(s) || !nzchar(s)) character(0) else strsplit(s, ", ", fixed = TRUE)[[1]]
}

# ---- Run LEWIS (complexity) on the dataset of one seed ----------------------
# Returns a list with the fit and the tables of the study.
run_uq <- function(seed) {
  # Draw the dataset. draw_dataset() sets its own seed.
  dat <- uq_dataset(seed)
  X <- dat$X

  # Set the seed of the chain and repeat the first TL computation.
  set.seed(seed)
  start0 <- uq_tl_start(X)

  # Run LEWIS (complexity) with the settings of the paper. The data are not
  # centered, as in the paper.
  fit <- lewis(X, d = LEWIS_D, alpha = LEWIS_ALPHA, alpha_sw = LEWIS_ALPHA_SW,
               alpha_ball = UQ_ALPHA_BALL, center = FALSE, ess_min = 0)
  if (!identical(unname(start0), unname(fit$diagnostics$start_dag))) {
    message("Note: the two TL computations gave different start graphs.")
  }

  # The true RPDAG and the SHD of the point estimate from it.
  truth <- unname(rpdag(dat$dag, dat$nG_star))
  shd_truth <- LEWIS::shd(unname(fit$D_hat), truth)
  radius <- fit$ball$radius
  M <- fit$diagnostics$n_draws

  # Table of all distinct RPDAGs, in decreasing order of posterior mass.
  g <- fit$graphs
  g_nG <- lapply(g$nongaussian, .split_nodes)
  top <- data.frame(
    rank = g$rank,
    count = as.integer(round(g$mass * M)),
    mass = g$mass,
    cum_mass = cumsum(g$mass),
    is_D_hat = g$class == "identical",
    shd_to_Dhat = g$shd,
    n_edges = g$n_edges,
    category = g$class,
    inside_ball = g$in_ball,
    stringsAsFactors = FALSE)

  # Table of the most probable graphs with their edits relative to D_hat.
  # The non-Gaussian set of a graph is the most frequent set among its draws.
  k <- seq_len(min(UQ_TOP_K, nrow(g)))
  uq <- data.frame(
    rank = g$rank[k],
    mass = g$mass[k],
    cum_mass = cumsum(g$mass)[k],
    is_D_hat = g$class[k] == "identical",
    shd_to_Dhat = g$shd[k],
    inside_ball = g$in_ball[k],
    category = g$class[k],
    n_edges = g$n_edges[k],
    n_edges_added = .n_entries(g$added[k]),
    n_edges_removed = .n_entries(g$removed[k]),
    n_edges_reversed = .n_entries(g$reversed[k]),
    added_edges_str = g$added[k],
    removed_edges_str = g$removed[k],
    reversed_edges_str = g$reversed[k],
    nongaussian = g$nongaussian[k],
    nG_added_str = vapply(g_nG[k], function(v) paste(setdiff(v, fit$nG_hat), collapse = ","),
                          character(1)),
    nG_removed_str = vapply(g_nG[k], function(v) paste(setdiff(fit$nG_hat, v), collapse = ","),
                            character(1)),
    stringsAsFactors = FALSE)

  # One-row summary of the run.
  summary_row <- data.frame(
    seed = seed,
    p = UQ_P,
    n = UQ_N,
    nG_star = paste(dat$nG_star, collapse = ","),
    n_iter = fit$diagnostics$n_iter,
    burn_in = fit$diagnostics$burn_in,
    thin = fit$diagnostics$thin,
    M = M,
    acceptance_rate = round(fit$diagnostics$acceptance_rate, 4),
    distinct_dags = fit$diagnostics$n_unique_dags,
    distinct_rpdags = fit$diagnostics$n_unique_rpdags,
    D_hat_edges = nrow(fit$edges),
    D_hat_mean_shd = round(sum(g$mass * g$shd), 3),
    nG_hat = paste(fit$nG_hat, collapse = ","),
    eps_alpha = radius,
    mass_outside_ball = fit$ball$mass_outside,
    shd_Dhat_truth = shd_truth,
    truth_in_ball_shd = shd_truth <= radius,
    alpha_ball = UQ_ALPHA_BALL,
    m_max = round(max(g$mass), 4),
    ess_log_score = round(fit$diagnostics$ess_log_score, 1),
    runtime_s = round(fit$diagnostics$runtime, 2),
    stringsAsFactors = FALSE)

  # Row of the seed sweep (choose_seed_for_demo.R).
  sweep_row <- data.frame(
    seed = seed,
    accept_rate = round(fit$diagnostics$acceptance_rate, 4),
    n_dags = fit$diagnostics$n_unique_dags,
    n_rpdags = fit$diagnostics$n_unique_rpdags,
    m_max = round(max(g$mass), 4),
    eps_alpha = radius,
    shd_to_truth = shd_truth,
    truth_in_ball = shd_truth <= radius,
    n_buckets = length(unique(uq$category)),
    nG_correct = identical(sort(as.integer(fit$nG_hat)), sort(as.integer(dat$nG_star))),
    stringsAsFactors = FALSE)

  list(fit = fit, summary = summary_row, top = top, uq = uq,
       edges = fit$edges, nodes = fit$nodes, sweep_row = sweep_row)
}

# ---- Write the tables of one run ---------------------------------------------
write_uq_tables <- function(res, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(res$summary, file.path(out_dir, "decision_rule_summary.csv"), row.names = FALSE)
  write.csv(res$uq, file.path(out_dir, "posterior_uncertainty_table.csv"), row.names = FALSE)
  write.csv(res$top, file.path(out_dir, "top_rpdag_posterior.csv"), row.names = FALSE)
  write.csv(res$edges, file.path(out_dir, "posterior_edges.csv"), row.names = FALSE)
  write.csv(res$nodes, file.path(out_dir, "posterior_nodes.csv"), row.names = FALSE)
  cat("Wrote the tables to", out_dir, "\n")
}
