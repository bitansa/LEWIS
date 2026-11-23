#' MCMC schedule of LEWIS
#'
#' Computes the number of iterations, the burn-in and the thinning interval
#' that [lewis()] uses for `p` variables.
#'
#' @details
#' By default the chain runs for
#' \eqn{10000 \max(1, (p/20)^2)}{10000 * max(1, (p/20)^2)} iterations,
#' rounded to a whole number, as in the paper. The first
#' `floor(burn_in * n_iter)` iterations are the burn-in. After the burn-in,
#' every `thin`-th state is kept, with
#' `thin = max(1, ceiling((n_iter - burn_in_iterations) / max_draws))`, so
#' that at most `max_draws` draws are kept. For example, the schedule is
#' 10,000 iterations with thinning 1 for p up to 20, 62,500 iterations
#' with burn-in 12,500 and thinning 7 for p = 50, and 250,000 iterations
#' with thinning 25 for p = 100.
#'
#' With `max_draws = Inf` there is no thinning. Every state after the
#' burn-in is kept. A long chain then needs more memory and a longer
#' post-processing.
#'
#' @param p Whole number, at least 2. The number of variables.
#' @param n_iter `NULL` or a positive whole number. The number of
#'   iterations. `NULL` (the default) uses the schedule of the paper.
#' @param burn_in Number in \eqn{[0, 1)}. The fraction of the iterations
#'   discarded as burn-in. The default is 0.2.
#' @param max_draws Positive whole number, or `Inf`. The largest number of
#'   draws kept after the burn-in. The default 8000 is the value of the
#'   paper. `Inf` keeps every draw after the burn-in.
#'
#' @return A list with elements
#' \describe{
#'   \item{p}{The number of variables.}
#'   \item{n_iter}{The number of iterations.}
#'   \item{burn_in}{The number of burn-in iterations.}
#'   \item{thin}{The thinning interval.}
#'   \item{n_draws}{The number of draws kept. It is the whole part of
#'     `(n_iter - burn_in) / thin`.}
#' }
#'
#' @references
#' Chaudhuri, A., Bhattacharya, A., Ni, Y. and Sarkar, B. (2026). Consistent
#' DAG selection for Bayesian Causal Discovery under general error
#' distributions. Manuscript.
#'
#' @seealso [lewis()]
#'
#' @examples
#' lewis_schedule(10)
#' lewis_schedule(50)
#' lewis_schedule(100)
#' lewis_schedule(5, n_iter = 2000)
#' lewis_schedule(5, n_iter = 20000, max_draws = Inf)
#'
#' @export
lewis_schedule <- function(p, n_iter = NULL, burn_in = 0.2, max_draws = 8000) {
  # Check the number of variables and the number of iterations.
  .check_count(p, "p", min = 2)
  if (!is.null(n_iter)) .check_count(n_iter, "n_iter", min = 1)

  # Check the burn-in fraction.
  .check_number(burn_in, "burn_in")
  if (burn_in < 0 || burn_in >= 1) {
    stop("burn_in must be a fraction in [0, 1).", call. = FALSE)
  }

  # Check max_draws. Inf is allowed and keeps every draw after the burn-in.
  keep_all <- is.numeric(max_draws) && length(max_draws) == 1L &&
    !is.na(max_draws) && max_draws == Inf
  if (!keep_all) {
    if (!is.numeric(max_draws) || length(max_draws) != 1L || !is.finite(max_draws) ||
        max_draws != round(max_draws) || max_draws < 1) {
      stop("max_draws must be a whole number of at least 1, or Inf to keep every draw.",
           call. = FALSE)
    }
  }

  # Number of iterations: the schedule of the paper unless given.
  if (is.null(n_iter)) {
    total <- round(10000 * max(1, (as.numeric(p) / 20)^2))
    if (total > .Machine$integer.max) {
      stop("p is too large for the default schedule. Give n_iter explicitly.",
           call. = FALSE)
    }
    n_iter <- as.integer(total)
  } else {
    n_iter <- as.integer(n_iter)
  }

  # Burn-in: the first floor(n_iter * burn_in) iterations.
  burn <- as.integer(floor(n_iter * burn_in))
  if (n_iter - burn < 1L) {
    stop("No draws are left after the burn-in. Increase n_iter or reduce burn_in.",
         call. = FALSE)
  }

  # Thinning: keep at most max_draws draws, or every draw if max_draws = Inf.
  thin <- if (keep_all) {
    1L
  } else {
    max(1L, as.integer(ceiling((n_iter - burn) / max_draws)))
  }
  list(p = as.integer(p), n_iter = n_iter, burn_in = burn, thin = thin,
       n_draws = as.integer((n_iter - burn) %/% thin))
}
