# =============================================================================
# laplace_exact.R
#
# Tools for Appendix F of the paper (accuracy of the Laplace approximation).
# They are used only by 06_laplace_vs_exact and are not part of the package
# LEWIS.
#
#     load_exact_engine()        compile laplace_exact.cpp and connect it to
#                                the Fortran routines of quantreg
#     tl_start_dag(X)            the start graph of both chains: the estimate
#                                of the TL method
#     laplace_chain(X, start, prior, n_iter)
#                                the chain of LEWIS with the Laplace
#                                approximation. Returns the highest-scoring
#                                DAG.
#     exact_chain(X, start, prior, n_iter, M)
#                                the exact (pseudo-marginal) chain with M
#                                importance draws per node. Returns the
#                                highest-scoring DAG.
#     cpdag_shd(A, B)            SHD between the CPDAGs of two DAGs
#     edge_string(A)             the edges of a DAG as "i>j;k>l"
#
# THE LAPLACE CHAIN uses the sampler of the package LEWIS through its
# unexported function LEWIS:::lewis_mcmc_cpp(). The exported function
# lewis() does not return the highest-scoring DAG of the chain, which
# Appendix F compares. The sampler is the one of lewis(), with the same
# random numbers.
#
# THE EXACT CHAIN replaces the Laplace approximation of each node by an
# unbiased importance-sampling estimate of the marginal likelihood of the
# node. A move draws new estimates for the nodes whose parent sets change
# and keeps the other estimates. The stationary distribution is the exact
# posterior (Andrieu and Roberts, 2009). The chain is written in C++
# (laplace_exact.cpp). A version in plain R gives the same draws but takes
# 50 to 160 times as long, which would also change the time ratios of
# Table 1.
#
# PRIORS of the importance sampler: b_j ~ N(0, 100 I) and
# theta_j ~ Inv-Gamma(1, 1), as in the model of the paper. The importance
# distributions are a multivariate t with 5 degrees of freedom, centered at
# the L1 fit b_hat, for b_j, and the log-normal distribution
# LogNormal(log theta_hat - c_n / 2, c_n) with c_n = log(1 + 1/n) for
# theta_j.
#
# The calling script must define COMMON_DIR and source settings.R and
# graph_tools.R first.
# Dependencies: LEWIS, Rcpp, quantreg, igraph, pcalg.
# =============================================================================

stopifnot(exists("COMMON_DIR"), exists("LEWIS_ALPHA"), exists("to_cpdag"))

# -----------------------------------------------------------------------------
# load_exact_engine: compile the C++ code once per R session.
# -----------------------------------------------------------------------------
# Call it once, before any parallel work, so that forked workers share the
# compiled code. The compiled library is cached in the folder given by the
# environment variable LEWIS_EXACT_CACHE, or in a temporary folder.
load_exact_engine <- function() {
  if (exists("exact_mcmc_cpp", mode = "function")) return(invisible(TRUE))
  for (pk in c("Rcpp", "quantreg")) {
    if (!requireNamespace(pk, quietly = TRUE)) stop("The package ", pk, " is needed.")
  }
  # quantreg must be loaded, because the C++ code calls its Fortran routines.
  loadNamespace("quantreg")
  cache <- Sys.getenv("LEWIS_EXACT_CACHE", "")
  if (!nzchar(cache)) cache <- file.path(tempdir(), "laplace_exact")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(COMMON_DIR, "laplace_exact.cpp"), env = globalenv(),
                  cacheDir = cache, showOutput = FALSE, verbose = FALSE)
  # Hand over the routines rqbr and rqfnb of quantreg.
  sym <- function(nm) getNativeSymbolInfo(nm, PACKAGE = "quantreg")$address
  exact_set_quantreg_symbols(sym("rqbr"), sym("rqfnb"))
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Small helpers.
# -----------------------------------------------------------------------------

# TRUE if the 0/1 matrix A has no directed cycle.
.is_acyclic <- function(A) {
  igraph::is_dag(igraph::graph_from_adjacency_matrix(A, mode = "directed"))
}

# A plain double matrix without names, as the compiled code expects.
.as_plain_double <- function(X) {
  X <- matrix(as.numeric(X), nrow(X), ncol(X))
  X
}

# The constant d of a prior. The uniform prior is d = 0.
.prior_d <- function(prior) if (prior == "complexity") LEWIS_D else 0

# -----------------------------------------------------------------------------
# tl_start_dag: the start graph of both chains.
# -----------------------------------------------------------------------------
# The estimate of LEWIS::tl_lingam (the TL method). Its matrix B has
# B[i, j] != 0 for the edge i -> j. If TL fails, the message is printed
# and the empty graph is used. If the graph has a cycle, edges are removed
# at random until it is acyclic. This is the start of the chains of the
# paper, with the same random numbers.
tl_start_dag <- function(X) {
  p <- ncol(X)
  empty <- matrix(0L, p, p)
  # Run TL and read off the edges.
  A <- tryCatch({
    fit <- LEWIS::tl_lingam(as.matrix(X))
    if (is.null(fit) || is.null(fit$B)) stop("tl_lingam returned no matrix B.")
    Bmat <- as.matrix(fit$B)
    Bmat[!is.finite(Bmat)] <- 0
    out <- matrix(0L, p, p)
    for (i in seq_len(p)) for (j in seq_len(p)) {
      if (i != j && abs(Bmat[i, j]) > 1e-8) out[i, j] <- 1L
    }
    out
  }, error = function(e) {
    message("TL failed and the chains start at the empty graph: ", conditionMessage(e))
    empty
  })
  storage.mode(A) <- "integer"
  # Remove random edges until the graph is acyclic.
  while (!.is_acyclic(A)) {
    ij <- which(A == 1L, arr.ind = TRUE)
    if (nrow(ij) == 0L) break
    drop <- ij[sample.int(nrow(ij), 1L), ]
    A[drop[1L], drop[2L]] <- 0L
  }
  A
}

# -----------------------------------------------------------------------------
# laplace_chain: the chain of LEWIS with the Laplace approximation.
# -----------------------------------------------------------------------------
# The schedule is lewis_schedule(p, n_iter). The highest-scoring DAG is
# taken over all iterations, including the start graph.
laplace_chain <- function(X, start, prior = c("complexity", "uniform"), n_iter) {
  prior <- match.arg(prior)
  sched <- lewis_schedule(ncol(X), n_iter = n_iter)
  st <- matrix(as.integer(start), ncol(X), ncol(X))
  res <- LEWIS:::lewis_mcmc_cpp(.as_plain_double(X), st, sched$n_iter,
                                sched$burn_in, sched$thin, LEWIS_ALPHA,
                                .prior_d(prior), FALSE)
  list(dag = res$best_dag, score = res$best_score,
       acceptance_rate = res$acceptance_rate)
}

# -----------------------------------------------------------------------------
# exact_chain: the pseudo-marginal chain.
# -----------------------------------------------------------------------------
# The highest-scoring DAG is taken over all iterations, including the start
# graph. The score is the importance-sampling estimate.
exact_chain <- function(X, start, prior = c("complexity", "uniform"), n_iter, M) {
  prior <- match.arg(prior)
  load_exact_engine()
  st <- matrix(as.integer(start), ncol(X), ncol(X))
  res <- exact_mcmc_cpp(.as_plain_double(X), st, as.integer(n_iter), prior,
                        LEWIS_ALPHA, .prior_d(prior), as.integer(M), FALSE)
  list(dag = res$best_dag, score = res$best_score,
       acceptance_rate = res$acceptance_rate)
}

# -----------------------------------------------------------------------------
# cpdag_shd and edge_string.
# -----------------------------------------------------------------------------

# SHD between the CPDAG of the DAG A and the CPDAG of the DAG B.
cpdag_shd <- function(A, B) {
  LEWIS::shd(to_cpdag(A), to_cpdag(B))
}

# The edges of a DAG as "i>j;k>l;...", ordered by parent and then child.
edge_string <- function(A) {
  idx <- which(A == 1L, arr.ind = TRUE)
  if (nrow(idx) == 0L) return("")
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  paste(sprintf("%d>%d", idx[, 1], idx[, 2]), collapse = ";")
}
