# The Metropolis-Hastings sampler and its starting graph.

#' Warm start from the TL method
#'
#' Runs `tl_lingam()`, the TL method. The coefficient matrix `B` has
#' `B[i, j] != 0` for the edge `i -> j`. The row is the parent. If TL
#' fails, the empty graph is returned together with the error message. If
#' the TL graph has a cycle, edges are removed at random until it is
#' acyclic. This makes sure that the start graph is always a DAG.
#' @noRd
.warm_start_tl <- function(X) {
  p <- ncol(X)
  empty <- matrix(0L, p, p)
  failed <- NULL

  # Run TL and read off the edges.
  A <- tryCatch({
    fit <- tl_lingam(unname(X))
    if (is.null(fit) || is.null(fit$B)) stop("The TL method returned no coefficient matrix.")
    Bmat <- as.matrix(fit$B)
    Bmat[!is.finite(Bmat)] <- 0
    out <- matrix(0L, p, p)
    for (i in seq_len(p)) for (j in seq_len(p)) {
      if (i != j && abs(Bmat[i, j]) > 1e-8) out[i, j] <- 1L
    }
    out
  }, error = function(e) {
    failed <<- conditionMessage(e)
    empty
  })
  storage.mode(A) <- "integer"

  # Remove random edges until the graph is acyclic.
  while (!isTRUE(lewis_is_acyclic_cpp(A))) {
    ij <- which(A == 1L, arr.ind = TRUE)
    if (nrow(ij) == 0L) break
    drop <- ij[sample.int(nrow(ij), 1L), ]
    A[drop[1L], drop[2L]] <- 0L
  }
  list(dag = A, failed = failed)
}

#' Run one chain
#'
#' `start` is "TL", "empty" or an integer adjacency matrix. Returns the
#' output of the compiled sampler together with the start graph and a
#' label of the start actually used.
#' @noRd
.run_chain <- function(X, start, sched, d, alpha, verbose = FALSE) {
  p <- ncol(X)

  # The starting graph.
  failed <- NULL
  if (is.matrix(start)) {
    start_dag <- start
    used <- "user"
  } else if (identical(start, "empty")) {
    start_dag <- matrix(0L, p, p)
    used <- "empty"
  } else {
    ws <- .warm_start_tl(X)
    start_dag <- ws$dag
    failed <- ws$failed
    used <- if (is.null(failed)) "TL" else "empty (TL failed)"
  }

  # The chain. alpha is not used when d = 0.
  alpha_used <- if (d == 0) 0.75 else alpha
  res <- lewis_mcmc_cpp(X, start_dag, sched$n_iter, sched$burn_in, sched$thin,
                        alpha_used, d, isTRUE(verbose))
  res$start_dag <- start_dag
  res$start_used <- used
  res$start_failed <- failed
  res
}
