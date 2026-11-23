# Convergence diagnostics of the chain.

#' Effective sample size
#'
#' Geyer's initial positive sequence estimator of the integrated
#' autocorrelation time. The autocorrelations are summed in pairs of
#' consecutive lags until a pair sum is not positive. The effective sample
#' size is the number of draws divided by the integrated autocorrelation
#' time. A constant series gets the number of draws.
#' @param x Numeric vector, one value per retained draw.
#' @return A single number.
#' @noRd
.ess_geyer <- function(x) {
  N <- length(x)
  # Degenerate series.
  if (N < 3L || stats::var(x) == 0) return(N)
  # Autocorrelations up to lag 2000. r[1] is lag 0.
  r <- stats::acf(x, lag.max = min(N - 1L, 2000L), plot = FALSE)$acf[, 1, 1]
  # Sum the pairs of autocorrelations while they are positive.
  tau <- -1
  m <- 0L
  while (2L * m + 2L <= length(r)) {
    G <- r[2L * m + 1L] + r[2L * m + 2L]
    if (G <= 0) break
    tau <- tau + 2 * G
    m <- m + 1L
  }
  N / max(tau, 1)
}

#' Potential scale reduction factor
#'
#' The potential scale reduction factor of Gelman and Rubin for a list of
#' chains. All chains are cut to the length of the shortest one.
#' @param chains List of numeric vectors, one per chain.
#' @return A single number. It is 1 when every chain is constant at the
#'   same value.
#' @noRd
.psrf <- function(chains) {
  # Put the chains side by side.
  N <- min(lengths(chains))
  if (N < 2L) return(NA_real_)
  x <- vapply(chains, function(v) as.numeric(v[seq_len(N)]), numeric(N))
  # Within-chain and between-chain variances.
  W <- mean(apply(x, 2L, stats::var))
  B <- N * stats::var(colMeans(x))
  if (W == 0) return(if (B == 0) 1 else Inf)
  sqrt(((N - 1) / N * W + B / N) / W)
}
