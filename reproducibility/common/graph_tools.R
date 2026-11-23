# =============================================================================
# graph_tools.R
#
# Graph functions used to evaluate the methods that are not part of the
# package LEWIS.
#
#     to_cpdag(dag)              CPDAG of a DAG, with pcalg::dag2cpdag
#     sid_bounds(est, true_dag)  lower and upper SID bounds, with the
#                                package SID
#     n_edges_partial(A)         number of edges of a partially directed
#                                graph
#     edge_list(A)               the edges of a graph as text
#
# The RPDAG of a DAG and the SHD between two graphs come from the package
# LEWIS: LEWIS::rpdag() and LEWIS::shd().
#
# ADJACENCY CONVENTION. Every matrix in this folder has A[k, j] = 1 for the
# directed edge k -> j (the row is the parent). An undirected edge k - j has
# A[k, j] = A[j, k] = 1. This is also the convention of the package SID and
# of the graphNEL to matrix coercion of the package graph.
#
# Dependencies: pcalg (to_cpdag), SID (sid_bounds).
# =============================================================================

# -----------------------------------------------------------------------------
# to_cpdag: the CPDAG of a DAG.
# -----------------------------------------------------------------------------
# The input must be a DAG. The matrix is first stripped of any class and
# names, so that the coercion to graphNEL works for every input.
to_cpdag <- function(dag) {
  if (!requireNamespace("pcalg", quietly = TRUE)) {
    stop("The package pcalg is needed for to_cpdag().")
  }
  # Plain integer matrix without names.
  dag <- matrix(as.integer(dag), nrow = nrow(dag), ncol = ncol(dag))
  # DAG to graphNEL, CPDAG, and back to a matrix.
  g <- methods::as(dag, "graphNEL")
  M <- methods::as(pcalg::dag2cpdag(g), "matrix")
  M <- matrix(as.integer(M != 0), nrow(M), ncol(M))
  M
}

# -----------------------------------------------------------------------------
# sid_bounds: structural intervention distance (Peters and Buhlmann, 2015).
# -----------------------------------------------------------------------------
# Returns list(sid_lower, sid_upper) for the estimate `est_graph` against
# the true DAG. When the estimate has undirected edges, SID returns the
# lower and upper bounds over the DAGs consistent with it. For a DAG the two
# bounds are equal. There is no time limit. An error of SID is not caught:
# it stops the caller, which reports it.
sid_bounds <- function(est_graph, true_dag) {
  if (!requireNamespace("SID", quietly = TRUE)) {
    stop("The package SID is needed for sid_bounds().")
  }
  if (nrow(est_graph) != nrow(true_dag)) {
    stop("The estimate and the true DAG have different numbers of nodes.")
  }
  # Plain integer matrices without names.
  est <- matrix(as.integer(est_graph), nrow(est_graph), ncol(est_graph))
  tru <- matrix(as.integer(true_dag), nrow(true_dag), ncol(true_dag))
  res <- SID::structIntervDist(trueGraph = tru, estGraph = est)
  # SID returns sid for a DAG, and sidLowerBound and sidUpperBound for a
  # graph with undirected edges.
  lower <- if (!is.null(res$sidLowerBound)) res$sidLowerBound else res$sid
  upper <- if (!is.null(res$sidUpperBound)) res$sidUpperBound else res$sid
  if (is.null(lower) || is.null(upper)) {
    stop("SID::structIntervDist returned no SID value.")
  }
  list(sid_lower = as.numeric(lower), sid_upper = as.numeric(upper))
}

# -----------------------------------------------------------------------------
# n_edges_partial: number of edges of a partially directed graph.
# -----------------------------------------------------------------------------
# A directed edge counts once. An undirected edge is stored twice and also
# counts once.
n_edges_partial <- function(A) {
  A <- matrix(as.integer(A != 0), nrow(A), ncol(A))
  und <- sum(A == 1L & t(A) == 1L) / 2L
  dir <- sum(A == 1L & t(A) == 0L)
  as.integer(und + dir)
}

# -----------------------------------------------------------------------------
# edge_list: the edges of a graph as text, for example "1->3  2--4".
# -----------------------------------------------------------------------------
edge_list <- function(A) {
  p <- nrow(A)
  out <- character(0)
  for (i in seq_len(p)) for (j in seq_len(p)) if (i < j) {
    if (A[i, j] != 0 && A[j, i] != 0) {
      out <- c(out, sprintf("%d--%d", i, j))
    } else if (A[i, j] != 0) {
      out <- c(out, sprintf("%d->%d", i, j))
    } else if (A[j, i] != 0) {
      out <- c(out, sprintf("%d->%d", j, i))
    }
  }
  if (length(out) == 0L) "(empty)" else paste(out, collapse = "  ")
}
