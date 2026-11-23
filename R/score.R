#' Score of a DAG under LEWIS
#'
#' Computes the score \eqn{s_n(\gamma)} of one DAG, as used by the sampler
#' of [lewis()]. The score is the log of the approximate posterior of the
#' DAG up to an additive constant.
#'
#' @details
#' LEWIS uses a Laplace error working model. Each variable is a linear
#' function of its parents plus an error with a Laplace distribution. The
#' coefficients have the prior \eqn{N(0, 100 I)} and each error scale has
#' the prior Inv-Gamma(1, 1). These priors are fixed, as in the paper.
#'
#' The marginal likelihood of a DAG is replaced by its Laplace
#' approximation (the Laplace approximation theorem of the paper). For node
#' \eqn{j} with parent set \eqn{S}, the score is
#' \deqn{s_{n,j}(S) = -n(1 + \log 2) - n \log \hat\theta_j - \frac{|S| + 1}{2} \log n + \hat c_j - n^\alpha d |S|,}{s_nj(S) = -n (1 + log 2) - n log(theta_j) - (|S| + 1) / 2 log(n) + c_j - n^alpha d |S|,}
#' where \eqn{\hat\theta_j} is the mean absolute residual of the least
#' absolute deviation regression of node \eqn{j} on its parents and
#' \eqn{\hat c_j} estimates the constant of the Laplace approximation. The
#' estimate \eqn{\hat c_j} uses a Gaussian kernel estimate of the residual
#' density at zero with Silverman's bandwidth. With fewer than 5 residuals
#' this density is set to 1. The score of the DAG is the sum of the node
#' scores. The last term comes from the complexity prior. With `d = 0` it
#' is zero. This gives the uniform prior.
#'
#' The model has no intercept. With `center = TRUE` (the default) the
#' column means of `X` are subtracted first.
#'
#' @param X Numeric matrix or data frame with `n` rows (observations) and
#'   `p` columns (variables). It must contain only finite values, and
#'   `n > p` is required.
#' @param dag A `p` x `p` adjacency matrix of a DAG with entries 0 and 1.
#'   `dag[i, j] = 1` means the edge `i -> j`.
#' @param d Number, at least 0. The constant of the complexity prior. The
#'   default 0.03 is the value of the paper. `d = 0` gives the uniform prior.
#' @param alpha Number strictly between 0.5 and 1. The exponent of the
#'   complexity prior. The default 0.75 is the value of the paper. It is
#'   not used when `d = 0`.
#' @param center Logical. If `TRUE` (the default), the columns of `X` are
#'   centered before the score is computed.
#'
#' @return A single number, the score \eqn{s_n} of `dag`. The attribute
#'   `per_node` holds the score of each node, named by the variables.
#'
#' @references
#' Chaudhuri, A., Bhattacharya, A., Ni, Y. and Sarkar, B. (2026). Consistent
#' DAG selection for Bayesian Causal Discovery under general error
#' distributions. Manuscript.
#'
#' Barrodale, I. and Roberts, F. D. K. (1974). Solution of an overdetermined
#' system of equations in the l1 norm. Communications of the ACM, 17,
#' 319-320.
#'
#' Koenker, R. (2005). Quantile Regression. Cambridge University Press.
#'
#' @seealso [lewis()], [nongaussian_nodes()]
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- 0.8 * x1 + rnorm(n)
#' X <- cbind(x1, x2)
#' A <- matrix(c(0, 0, 1, 0), 2, 2)  # the edge x1 -> x2
#' lewis_score(X, A)
#' lewis_score(X, matrix(0, 2, 2))   # the empty graph
#'
#' @export
lewis_score <- function(X, dag, d = 0.03, alpha = 0.75, center = TRUE) {
  # Check the arguments.
  .check_flag(center, "center")
  dat <- .check_data(X, center = center)
  A <- .check_dag(dag, dat$p, dat$names)
  .check_prior(d, alpha)
  alpha_used <- if (d == 0) 0.75 else alpha

  # Score every node and sum.
  s <- lewis_score_dag_cpp(dat$X, A, alpha_used, d)
  out <- s$total
  attr(out, "per_node") <- stats::setNames(s$per_node, dat$names)
  out
}
