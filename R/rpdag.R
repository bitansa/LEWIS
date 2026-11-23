#' Non-Gaussian nodes of a DAG
#'
#' Step 1 of the post-processing of LEWIS (estimating the non-Gaussian
#' nodes). Estimates which nodes of a DAG have non-Gaussian errors.
#'
#' @details
#' For each node, the least absolute deviation regression of the node on
#' its parents in `dag` gives residuals. The residuals of a node without
#' parents are its data (centered if `center = TRUE`), because the model
#' has no intercept. The Shapiro-Wilk test of normality is
#' applied to the residuals of each node. A node is declared non-Gaussian
#' when the p-value is below `alpha_sw / p`, a Bonferroni correction over
#' the `p` nodes. With more than 5000 observations, the test uses the 5000
#' residuals at the positions `round(seq(1, n, length.out = 5000))`,
#' because [stats::shapiro.test()] accepts at most 5000 values. The test is
#' the algorithm of [stats::shapiro.test()].
#'
#' @inheritParams lewis_score
#' @param alpha_sw Number strictly between 0 and 1. The level of the
#'   Shapiro-Wilk tests before the Bonferroni division by `p`. The default
#'   is 0.05.
#'
#' @return A character vector with the names of the nodes declared
#'   non-Gaussian, in the order of the columns of `X`. The attribute
#'   `p_values` holds the Shapiro-Wilk p-value of every node. The attribute
#'   `level` holds the level `alpha_sw / p`.
#'
#' @references
#' Chaudhuri, A., Bhattacharya, A., Ni, Y. and Sarkar, B. (2026). Consistent
#' DAG selection for Bayesian Causal Discovery under general error
#' distributions. Manuscript.
#'
#' Shapiro, S. S. and Wilk, M. B. (1965). An analysis of variance test for
#' normality (complete samples). Biometrika, 52, 591-611.
#'
#' @seealso [rpdag()], [lewis()]
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' x1 <- rexp(n) - rexp(n)          # Laplace error
#' x2 <- 0.8 * x1 + rnorm(n)        # Gaussian error
#' X <- cbind(x1, x2)
#' A <- matrix(c(0, 0, 1, 0), 2, 2) # the edge x1 -> x2
#' nongaussian_nodes(X, A)
#'
#' @export
nongaussian_nodes <- function(X, dag, alpha_sw = 0.05, center = TRUE) {
  # Check the arguments.
  .check_flag(center, "center")
  .check_level(alpha_sw, "alpha_sw")
  dat <- .check_data(X, center = center)
  A <- .check_dag(dag, dat$p, dat$names)

  # Test every node.
  res <- lewis_nongaussian_cpp(dat$X, A, alpha_sw)
  out <- dat$names[res$nongaussian]
  attr(out, "p_values") <- stats::setNames(res$p_value, dat$names)
  attr(out, "level") <- alpha_sw / dat$p
  out
}

#' Restricted partially directed acyclic graph (RPDAG)
#'
#' Step 2 of the post-processing of LEWIS (projection onto the space of
#' RPDAGs). Projects a DAG and a set of non-Gaussian nodes onto the RPDAG
#' that represents their distribution equivalence class.
#'
#' @details
#' The RPDAG is built in three steps, as in the paper.
#' (A) Take the CPDAG of `dag`, which represents its Markov equivalence
#' class.
#' (B) For every undirected edge of the CPDAG with at least one endpoint in
#' `nG`, restore the direction it has in `dag`.
#' (C) Orient further undirected edges with Meek's rules R1 to R4 until no
#' rule applies.
#'
#' With an empty `nG` the result is the CPDAG of `dag`. When every node is
#' in `nG` the result is `dag` itself.
#'
#' @param dag A `p` x `p` adjacency matrix of a DAG with entries 0 and 1.
#'   `dag[i, j] = 1` means the edge `i -> j`.
#' @param nG The non-Gaussian nodes. Node indices, node names (matched to
#'   the column names of `dag`), or a logical vector of length `p`. The
#'   default is the empty set.
#'
#' @return The RPDAG as a `p` x `p` integer matrix with row and column
#'   names. `D[i, j] = 1` and `D[j, i] = 0` means the directed edge
#'   `i -> j`. `D[i, j] = D[j, i] = 1` means the undirected edge `i - j`.
#'
#' @references
#' Chaudhuri, A., Bhattacharya, A., Ni, Y. and Sarkar, B. (2026). Consistent
#' DAG selection for Bayesian Causal Discovery under general error
#' distributions. Manuscript.
#'
#' Hoyer, P. O., Hyvarinen, A., Scheines, R., Spirtes, P., Ramsey, J.,
#' Lacerda, G. and Shimizu, S. (2008). Causal discovery of linear acyclic
#' models with arbitrary distributions. Proceedings of the 24th Conference
#' on Uncertainty in Artificial Intelligence.
#'
#' Meek, C. (1995). Causal inference and causal explanation with background
#' knowledge. Proceedings of the 11th Conference on Uncertainty in
#' Artificial Intelligence.
#'
#' @seealso [nongaussian_nodes()], [shd()], [lewis()]
#'
#' @examples
#' # The chain 1 -> 2 -> 3.
#' A <- matrix(0, 3, 3)
#' A[1, 2] <- 1
#' A[2, 3] <- 1
#' rpdag(A)          # all Gaussian: the CPDAG, 1 - 2 - 3
#' rpdag(A, nG = 1)  # node 1 non-Gaussian: 1 -> 2 -> 3
#'
#' @export
rpdag <- function(dag, nG = integer(0)) {
  # Check the graph and the node set.
  if (is.data.frame(dag)) dag <- as.matrix(dag)
  if (!is.matrix(dag)) stop("dag must be a square adjacency matrix.", call. = FALSE)
  p <- nrow(dag)
  if (p < 2L) stop("dag must have at least 2 nodes.", call. = FALSE)
  nm <- colnames(dag) %||% rownames(dag) %||% paste0("V", seq_len(p))
  A <- .check_dag(dag, p)
  idx <- .check_nodes(nG, p, nm)

  # Build the RPDAG.
  D <- lewis_rpdag_cpp(A, idx)
  dimnames(D) <- list(nm, nm)
  D
}

#' Structural Hamming distance
#'
#' Step 3 of the post-processing of LEWIS (point estimate and credible
#' ball) uses this distance. Computes the structural Hamming distance (SHD)
#' between two partially directed graphs.
#'
#' @details
#' The SHD is the number of unordered pairs of nodes whose status differs
#' between the two graphs. The status of a pair \eqn{\{i, j\}} is one of: no
#' edge, \eqn{i \to j}{i -> j}, \eqn{j \to i}{j -> i}, or undirected. A
#' missing, extra, reversed or wrongly undirected edge therefore costs one.
#' The SHD is symmetric and is zero only for identical graphs.
#'
#' @param D1,D2 Square adjacency matrices of the same size with entries 0
#'   and 1 and a zero diagonal. `D[i, j] = 1` and `D[j, i] = 0` means
#'   `i -> j`. `D[i, j] = D[j, i] = 1` means the undirected edge `i - j`.
#'   DAGs, CPDAGs and RPDAGs are all allowed.
#'
#' @return A single integer, the SHD.
#'
#' @references
#' Tsamardinos, I., Brown, L. E. and Aliferis, C. F. (2006). The max-min
#' hill-climbing Bayesian network structure learning algorithm. Machine
#' Learning, 65, 31-78.
#'
#' @seealso [rpdag()], [lewis()]
#'
#' @examples
#' A <- matrix(0, 3, 3)
#' A[1, 2] <- 1
#' A[2, 3] <- 1
#' B <- t(A)          # both edges reversed
#' shd(A, B)          # 2
#' shd(A, rpdag(A))   # 2, both edges are undirected in the CPDAG
#' shd(A, A)          # 0
#'
#' @export
shd <- function(D1, D2) {
  # Check the two graphs.
  A <- .check_pdag(D1, "D1")
  B <- .check_pdag(D2, "D2")
  if (nrow(A) != nrow(B)) {
    stop("D1 and D2 must have the same number of nodes.", call. = FALSE)
  }
  # Count the pairs whose status differs.
  lewis_shd_cpp(A, B)
}
