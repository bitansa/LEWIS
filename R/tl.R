# The TL method of Zhao, He and Wang (2022).
#
# The functions in this file are copies of the functions TLLiNGAM,
# Compute.precision, Compute.stdPreciMat.refit.high, Compute.dcor and
# Parent.set.refit of the R package TransGraph, version 1.1.0, by Mingyang
# Ren, Ruixuan Zhao, Xin He and Junhui Wang, licensed under the GNU General
# Public License, version 2. TransGraph is no longer on CRAN. The code is
# unchanged except for the function names, the comments, a check of the
# input and the removal of the CLIME option. The default settings do not
# use this option. The results are the same as those of
# TransGraph::TLLiNGAM with its default settings.

#' Learn a linear non-Gaussian DAG with the TL method
#'
#' Estimates a DAG with the topological layer (TL) method of Zhao, He and
#' Wang (2022). LEWIS uses this estimate as the default start of its chain.
#' The function is a copy of `TLLiNGAM()` of the package TransGraph, version
#' 1.1.0, with its default settings.
#'
#' @details The method finds the layers of the DAG one at a time. At each
#'   step, it estimates a precision matrix by the graphical lasso, refits its
#'   coefficients by least squares, and tests with distance covariance which
#'   remaining nodes are independent of their regression residuals. The
#'   distance covariance test uses 100 random permutations. The result
#'   therefore depends on the state of the random number generator. Set a
#'   seed for a reproducible result.
#'
#' @param X Numeric matrix with `n` rows (observations) and `p` columns
#'   (variables).
#' @param hardth Hard threshold for the regression coefficients. The
#'   default is 0.3.
#' @param criti.val Critical value of the distance covariance test. The
#'   default is 0.01.
#' @param precision.refit Logical. If `TRUE` (the default), the coefficients
#'   of the precision matrix are refitted by least squares.
#' @param B.refit Logical. If `TRUE` (the default), the coefficients of the
#'   structural equations are refitted by least squares once the parents
#'   of every node are found.
#'
#' @return A list with three elements.
#' \describe{
#'   \item{A}{A `p` x 2 matrix. The first column is the node and the second
#'     column is its layer. The nodes found in the first step have layer 0.}
#'   \item{B}{The `p` x `p` matrix of coefficients. `B[i, j] != 0` means
#'     the edge `i -> j`. The row is the parent.}
#'   \item{check}{`TRUE` if no node passed the first test. Then `B` is
#'     zero.}
#' }
#'
#' @references
#' Zhao, R., He, X. and Wang, J. (2022). Learning linear non-Gaussian
#' directed acyclic graph with diverging number of nodes. Journal of
#' Machine Learning Research, 23.
#'
#' Ren, M., Zhao, R., He, X. and Wang, J. (2025). TransGraph: Transfer
#' Graph Learning. R package version 1.1.0.
#'
#' @seealso [lewis()]
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x1 <- rexp(n) - rexp(n)
#' x2 <- 0.8 * x1 + rexp(n) - rexp(n)
#' x3 <- 0.7 * x2 + rexp(n) - rexp(n)
#' fit <- tl_lingam(cbind(x1, x2, x3))
#' fit$B != 0
#'
#' @export
tl_lingam <- function(X, hardth = 0.3, criti.val = 0.01, precision.refit = TRUE,
                      B.refit = TRUE) {
  # Check the input.
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a numeric matrix of finite values.", call. = FALSE)
  }
  if (ncol(X) < 2L) stop("X must have at least 2 columns.", call. = FALSE)

  check <- FALSE
  n <- nrow(X)
  p <- ncol(X)

  # A stores the layer of every node. B stores the coefficients.
  A <- matrix(0, p, 2)
  A[, 1] <- 1:p
  colnames(A) <- c("Node", "Layer")
  B <- matrix(0, p, p)

  # Find the first layer: the nodes independent of their residuals.
  S <- 1:p
  stdPreciMat <- .tl_precision(X, hardth, precision.refit)
  dcor_value <- sapply(1:p, function(k) .tl_dcor(k, stdPreciMat[, k], X))
  node_S_select <- S[dcor_value > criti.val]
  if (length(node_S_select) == 0) {
    check <- TRUE
    return(list(A = A, B = B, check = check))
  }
  S <- S[-node_S_select]
  B[S, node_S_select] <- stdPreciMat[S, node_S_select]
  t <- 1

  # Find the next layers among the remaining nodes.
  while (length(S) > 1) {
    len_S <- length(S)
    X_wait <- X[, S]
    stdPreciMat <- .tl_precision(X_wait, hardth, precision.refit)
    dcor_value <- sapply(1:len_S, function(k) .tl_dcor(k, stdPreciMat[, k], X_wait))
    node_S_select <- S[dcor_value > criti.val]
    index_node_S_select <- which(dcor_value > criti.val)
    if (length(node_S_select) == 0) {
      node_S_select <- S
    }
    A[node_S_select, 2] <- t
    S <- setdiff(S, node_S_select)
    index_S <- setdiff(c(1:len_S), index_node_S_select)
    B[S, node_S_select] <- stdPreciMat[index_S, index_node_S_select]
    t <- t + 1
  }
  if (length(S) == 1) {
    A[S, 2] <- t
  }

  # Refit the coefficients given the parent sets.
  if (B.refit) {
    B <- .tl_refit_parents(X, B)
  }
  list(A = A, B = B, check = check)
}

# Standardized precision matrix by the graphical lasso, with an optional
# least-squares refit. Copy of TransGraph's Compute.precision.
.tl_precision <- function(comp.X, comp.hardth, precision.refit = TRUE) {
  n <- nrow(comp.X)
  p <- ncol(comp.X)
  # Graphical lasso estimate of the precision matrix.
  PreciMat <- huge::huge(comp.X, lambda = 0.5 * sqrt(log(max(p, n)) / n),
                         method = "glasso", verbose = FALSE)$icov[[1]]
  # Divide each column by its diagonal entry.
  stdPreciMat <- t(PreciMat / diag(PreciMat))
  # Refit by least squares and set small coefficients to zero.
  if (precision.refit) {
    stdPreciMat <- .tl_refit_precision(comp.X, stdPreciMat, comp.hardth)
    stdPreciMat[abs(stdPreciMat) < comp.hardth] <- 0
  }
  stdPreciMat
}

# Least-squares refit of the precision coefficients with a hard threshold.
# Copy of TransGraph's Compute.stdPreciMat.refit.high.
.tl_refit_precision <- function(comp.X, pre_stdPreciMat, hardth = 0.3) {
  p <- ncol(comp.X)
  stdPreciMat <- matrix(0, p, p)
  diag(stdPreciMat) <- 1
  diag(pre_stdPreciMat) <- 0
  for (i in 1:p) {
    # Regress node i on the nodes with a nonzero coefficient.
    index_regres <- which(pre_stdPreciMat[, i] != 0)
    if (length(index_regres) < nrow(comp.X)) {
      if (length(index_regres) != 0) {
        pre.coeff <- MASS::ginv(t(comp.X[, index_regres]) %*% comp.X[, index_regres]) %*%
          t(comp.X[, index_regres]) %*% comp.X[, i]
        # Refit on the coefficients above the threshold.
        index_select <- which(abs(pre.coeff) > hardth)
        if (length(index_select) != 0) {
          X_regres <- comp.X[, index_regres]
          if (is.matrix(X_regres) == TRUE) {
            select.coeff <- MASS::ginv(t(X_regres[, index_select]) %*% X_regres[, index_select]) %*%
              t(X_regres[, index_select]) %*% comp.X[, i]
            comp.coeff <- rep(0, length(index_regres))
            comp.coeff[index_select] <- select.coeff
            stdPreciMat[index_regres, i] <- -comp.coeff
          } else {
            select.coeff <- MASS::ginv(t(X_regres) %*% X_regres) %*% t(X_regres) %*% comp.X[, i]
            stdPreciMat[index_regres, i] <- -select.coeff
          }
        }
      }
    }
  }
  stdPreciMat
}

# P-values of the distance covariance tests between the residual of one
# node and every other node. Copy of TransGraph's Compute.dcor.
.tl_dcor <- function(index.Y, preci.coeff, comp.X) {
  p <- ncol(comp.X)
  # Residual of the node and the data of the other nodes.
  comp.resid <- comp.X %*% preci.coeff
  other.X <- comp.X[, -index.Y]
  # Distance covariance tests of independence.
  if (p > 2) {
    comp.dcor.index <- 1:ncol(other.X)
    comp.dcor.test <- dcov::dcor.test(comp.resid, other.X[, comp.dcor.index], R = 100)$p.values
  } else {
    comp.dcor.test <- dcov::dcor.test(comp.resid, other.X, R = 100)$p.values
  }
  comp.dcor.test
}

# Least-squares refit of the structural coefficients given the parent sets.
# Copy of TransGraph's Parent.set.refit.
.tl_refit_parents <- function(comp.X, pre_B) {
  p <- ncol(comp.X)
  B <- matrix(0, p, p)
  for (i in 1:p) {
    # Least-squares regression of node i on its parents.
    index <- which(pre_B[, i] != 0)
    if (length(index) > 0) {
      comp.coeff <- MASS::ginv(t(comp.X[, index]) %*% comp.X[, index]) %*%
        t(comp.X[, index]) %*% comp.X[, i]
      B[index, i] <- comp.coeff
    }
  }
  B
}
