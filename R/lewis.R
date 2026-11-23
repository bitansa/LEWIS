#' Bayesian causal discovery with LEWIS
#'
#' Fits LEWIS (Laplace Error Working model for Inference and Selection of
#' DAG) to data. LEWIS samples DAGs from an approximate posterior and
#' summarizes the draws by an estimate of the distribution equivalence
#' class, a credible ball and posterior probabilities.
#'
#' @details
#' **Model.** For a DAG, each variable is a linear function of its parents
#' plus an error. The working model takes the errors to be Laplace. The
#' true errors need not be Laplace. Each true error may be Gaussian or a
#' Gaussian scale mixture, such as a Student's t error. The errors of
#' different variables may have different distributions. The coefficients
#' have the prior \eqn{N(0, 100 I)} and the error scales have the prior
#' Inv-Gamma(1, 1). The DAG has the complexity
#' prior \eqn{\pi(\gamma) \propto \exp(-n^\alpha d |\gamma|)}{pi(gamma)
#' proportional to exp(-n^alpha d |gamma|)}, where \eqn{|\gamma|}{|gamma|}
#' is the number of edges. With `d = 0` this is the uniform prior. The
#' model has no intercept. With `center = TRUE` (the default) the column
#' means of `X` are subtracted first.
#'
#' **Sampling.** The marginal likelihood is replaced by its Laplace
#' approximation, which gives the score \eqn{s_n} of a DAG (see
#' [lewis_score()]). A Metropolis-Hastings sampler explores the DAGs. Each
#' iteration proposes an addition, a deletion or a reversal of one edge,
#' each with probability 1/3. An impossible move or a move that gives a
#' cycle is rejected. Otherwise the move is accepted with probability
#' \eqn{\min\{1, \exp(s_n(\gamma') - s_n(\gamma)) r\}}{min(1, exp(s_n(new) -
#' s_n(old)) r)}, where \eqn{r} is the Hastings ratio of the paper. The
#' number of iterations, the burn-in and the thinning follow
#' [lewis_schedule()].
#'
#' **Post-processing.** The retained draws are summarized in four steps,
#' as in the paper.
#' Step 1 (estimating the non-Gaussian nodes) declares the non-Gaussian
#' nodes of each draw by Shapiro-Wilk tests at the level `alpha_sw / p`
#' (see [nongaussian_nodes()]).
#' Step 2 (projection onto the space of RPDAGs) maps each draw to the RPDAG
#' of its distribution equivalence class (see [rpdag()]).
#' Step 3 (point estimate and credible ball) takes as point estimate
#' `D_hat` the distinct RPDAG with the smallest mean structural Hamming
#' distance (SHD, see [shd()]) to the projected draws. Ties go to fewer
#' edges and then to the graph that the chain visited first. The estimate
#' `nG_hat` of the non-Gaussian nodes is the most frequent set among the
#' draws whose RPDAG is `D_hat`. Ties go to the smaller set and then to the
#' set that the chain visited first. The credible ball is the set of graphs
#' within SHD `radius` of `D_hat`. Its radius is the smallest value for
#' which the ball holds at least a fraction `1 - alpha_ball` of the
#' projected draws.
#' Step 4 (the shape of the posterior) lists the distinct RPDAGs in
#' decreasing order of posterior mass. Each graph is compared with `D_hat`.
#' It is identical (the graph is `D_hat`), upper (every edge of `D_hat` is
#' an edge of the graph with the same status, and the graph is not
#' `D_hat`), lower (every edge of the graph is an edge of `D_hat` with the
#' same status, and the graph is not `D_hat`) or rearranged (otherwise).
#' The posterior probability of an edge of `D_hat` is the proportion of
#' projected draws in which the edge has the same status as in `D_hat`. Its
#' adjacency probability is the proportion of projected draws in which the
#' two nodes are adjacent. All probabilities are proportions of the
#' retained draws.
#'
#' **Convergence.** The effective sample size (ESS) of the score of the
#' retained draws is computed with Geyer's initial positive sequence
#' estimator. If the score is the same in all retained draws, the ESS is
#' set to the number of draws. A warning is given when the ESS is below
#' `ess_min`. With `n_chains > 1`, further chains are run. Chains 2, 4, ...
#' start at the empty graph and chains 3, 5, ... start at the TL estimate.
#' If `seed` is given, chain `k` uses the seed `seed + k - 1`. All estimates
#' come from chain 1. LEWIS reports the potential scale reduction factor
#' (PSRF) of the score and of the number of edges. It also reports the
#' largest difference between chains in the posterior probability of an
#' adjacency or an orientation. A warning is given when a PSRF is above
#' 1.05 or a difference is above 0.1. If a warning appears, increase
#' `n_iter`. The real-data analysis of the paper uses \eqn{10^7}{10^7}
#' iterations.
#'
#' **Small samples.** LEWIS needs `n > p`. The kernel estimate of the
#' residual density in the score needs at least 5 residuals. With fewer,
#' the density is set to 1.
#'
#' @param X Numeric matrix or data frame with `n` rows (observations) and
#'   `p` columns (variables). It must contain only finite values. Constant,
#'   duplicated or linearly dependent columns are not allowed. The column
#'   names are used in all outputs. Without names, the variables are called
#'   `V1`, ..., `Vp`.
#' @param d Number, at least 0. The constant of the complexity prior. The
#'   default 0.03 is the value of the paper. `d = 0` gives the uniform
#'   prior.
#' @param alpha Number strictly between 0.5 and 1. The exponent of the
#'   complexity prior. The default 0.75 is the value of the paper. It is
#'   not used when `d = 0`.
#' @param start The starting graph of the chain. `"TL"` (the default)
#'   starts at the estimate of the TL method, computed by [tl_lingam()],
#'   as in the paper. If TL fails, a warning is given and the chain starts
#'   at the empty graph. `"empty"` starts at the empty graph. A `p` x `p`
#'   adjacency matrix of a DAG with entries 0 and 1 (row = parent, column =
#'   child) starts at that DAG. If this matrix has column names, they must
#'   equal the variable names.
#' @param n_iter `NULL` or a positive whole number. The number of
#'   iterations. `NULL` (the default) uses the schedule of the paper,
#'   \eqn{10000 \max(1, (p/20)^2)}{10000 * max(1, (p/20)^2)}.
#' @param burn_in Number in \eqn{[0, 1)}. The fraction of the iterations
#'   discarded as burn-in. The default is 0.2.
#' @param max_draws Positive whole number, or `Inf`. The chain is thinned
#'   so that at most `max_draws` draws are kept after the burn-in. The
#'   default 8000 is the value of the paper. `Inf` keeps every draw after
#'   the burn-in.
#' @param alpha_sw Number strictly between 0 and 1. The level of the
#'   Shapiro-Wilk tests of Step 1 before the Bonferroni division by `p`.
#'   The default is 0.05.
#' @param alpha_ball Number strictly between 0 and 1. One minus the
#'   credible level of the ball. The default 0.05 gives the 95 percent
#'   credible ball.
#' @param center Logical. If `TRUE` (the default), the columns of `X` are
#'   centered.
#' @param n_chains Positive whole number. The number of chains. The default
#'   is 1, as in the paper. Use 2 or more to check convergence.
#' @param ess_min Number, at least 0. A warning is given when the effective
#'   sample size of the score is below `ess_min`. The default is 100.
#'   Use 0 to switch the check off.
#' @param seed `NULL` or a whole number. If given, [set.seed()] is called
#'   with it before chain 1 starts. The results can then be reproduced.
#'   Chain `k` uses the seed `seed + k - 1`. This changes the state of the
#'   random number generator of the R session. If `NULL` (the default), the
#'   current state of the random number generator is used and no seed is
#'   set.
#' @param verbose Logical. If `TRUE`, the progress of the chain is printed.
#'   The default is `FALSE`.
#'
#' @return An object of class `"lewis"`, a list with elements
#' \describe{
#'   \item{D_hat}{The point estimate, an RPDAG as a `p` x `p` integer
#'     matrix with row and column names. `D_hat[i, j] = 1` and
#'     `D_hat[j, i] = 0` means `i -> j`. `D_hat[i, j] = D_hat[j, i] = 1`
#'     means the undirected edge `i - j`.}
#'   \item{nG_hat}{Character vector. The names of the nodes declared
#'     non-Gaussian.}
#'   \item{edges}{Data frame with one row per edge of `D_hat` and columns
#'     `from`, `to`, `type` (`"directed"` or `"undirected"`), `prob` (the
#'     posterior probability of the edge with the same status as in
#'     `D_hat`) and `prob_adjacent` (the posterior probability that the two
#'     nodes are adjacent). An undirected edge is listed once.}
#'   \item{nodes}{Data frame with one row per node and columns `node`,
#'     `prob_nongaussian` (the proportion of draws in which the node is
#'     declared non-Gaussian) and `in_nG_hat` (logical).}
#'   \item{graphs}{Data frame with one row per distinct RPDAG among the
#'     projected draws, in decreasing order of posterior mass. Graphs with
#'     the same mass are in the order of the first visit. Columns: `rank`,
#'     `mass` (posterior mass), `class` (`"identical"`, `"upper"`,
#'     `"lower"` or `"rearranged"`), `shd` (SHD to `D_hat`), `n_edges`,
#'     `added`, `removed` and `reversed` (the edges that differ from
#'     `D_hat`, separated by `" | "`, with a reversed edge written as its
#'     orientation in `D_hat`, then `"=>"`, then its orientation in the
#'     graph), `nongaussian` (the most frequent set of non-Gaussian nodes
#'     among the draws of the graph, with ties broken as for `nG_hat`) and
#'     `in_ball` (logical, `TRUE` if the graph is in the credible ball).}
#'   \item{ball}{List with the credible ball: `radius` (the SHD radius),
#'     `alpha_ball`, and `mass_outside` (the proportion of draws outside
#'     the ball).}
#'   \item{trace}{Data frame with one row per retained draw and columns
#'     `iteration`, `log_score` (the score \eqn{s_n} of the DAG, that is,
#'     the log of the approximate posterior up to a constant) and `n_edges`
#'     (the number of edges of the DAG).}
#'   \item{diagnostics}{List with `acceptance_rate`, `invalid_rate` (the
#'     proportion of impossible or cyclic proposals), `moves` (proposed and
#'     accepted moves by type), `n_iter`, `burn_in` (iterations), `thin`,
#'     `n_draws`, `n_unique_dags`, `n_unique_rpdags`, `ess_log_score` (the
#'     ESS of the score), `best_score` (the highest score visited),
#'     `start` (the start used:
#'     `"TL"`, `"empty"`, `"user"` or `"empty (TL failed)"`), `start_dag`,
#'     `runtime` (seconds) and `warnings`. With `n_chains > 1` also
#'     `psrf_log_score`, `psrf_n_edges`, `max_diff_adjacency`,
#'     `max_diff_orientation` and `chains` (a data frame with one row per
#'     chain).}
#'   \item{settings}{List with the arguments used, and `n`, `p` and
#'     `variables`.}
#'   \item{call}{The matched call.}
#' }
#'
#' @references
#' Chaudhuri, A., Bhattacharya, A., Ni, Y. and Sarkar, B. (2026). Consistent
#' DAG selection for Bayesian Causal Discovery under general error
#' distributions. Manuscript.
#'
#' Geyer, C. J. (1992). Practical Markov chain Monte Carlo. Statistical
#' Science, 7, 473-483.
#'
#' Gelman, A. and Rubin, D. B. (1992). Inference from iterative simulation
#' using multiple sequences. Statistical Science, 7, 457-472.
#'
#' Wade, S. and Ghahramani, Z. (2018). Bayesian cluster analysis: point
#' estimation and credible balls. Bayesian Analysis, 13, 559-626.
#'
#' Zhao, R., He, X. and Wang, J. (2022). Learning linear non-Gaussian
#' directed acyclic graph with diverging number of nodes. Journal of
#' Machine Learning Research, 23.
#'
#' @seealso [lewis_schedule()], [lewis_score()], [nongaussian_nodes()],
#'   [rpdag()], [shd()], [tl_lingam()], [plot.lewis()], [summary.lewis()],
#'   [tcga_prostate]
#'
#' @examples
#' # Simulated data from the chain x1 -> x2 -> x3. The errors of x1 and x3
#' # are Laplace and the error of x2 is Gaussian.
#' set.seed(1)
#' n <- 300
#' x1 <- rexp(n) - rexp(n)
#' x2 <- 0.8 * x1 + rnorm(n)
#' x3 <- 0.8 * x2 + rexp(n) - rexp(n)
#' X <- cbind(x1, x2, x3)
#'
#' # Fit LEWIS with the default settings of the paper.
#' fit <- lewis(X, seed = 1)
#' fit
#' fit$D_hat
#' fit$edges
#'
#' @export
lewis <- function(X, d = 0.03, alpha = 0.75, start = "TL", n_iter = NULL,
                  burn_in = 0.2, max_draws = 8000, alpha_sw = 0.05,
                  alpha_ball = 0.05, center = TRUE, n_chains = 1,
                  ess_min = 100, seed = NULL, verbose = FALSE) {
  cl <- match.call()
  t0 <- proc.time()[["elapsed"]]

  # Check the scalar arguments.
  .check_flag(center, "center")
  .check_flag(verbose, "verbose")
  .check_prior(d, alpha)
  if (d == 0 && !missing(alpha)) {
    message("alpha is not used because d = 0 gives the uniform prior.")
  }
  .check_level(alpha_sw, "alpha_sw")
  .check_level(alpha_ball, "alpha_ball")
  n_chains <- .check_count(n_chains, "n_chains", min = 1)
  .check_number(ess_min, "ess_min")
  if (ess_min < 0) stop("ess_min must be zero or positive.", call. = FALSE)
  .check_seed(seed)

  # Check the data and the starting graph, and set the schedule.
  dat <- .check_data(X, center = center)
  st <- .check_start(start, dat$p, dat$names)
  sched <- lewis_schedule(dat$p, n_iter = n_iter, burn_in = burn_in,
                          max_draws = max_draws)
  Xc <- dat$X
  p <- dat$p
  warn <- character(0)

  # Set the seed once, before any random number is drawn.
  if (!is.null(seed)) set.seed(seed)

  # Chain 1.
  if (verbose && n_chains > 1L) message("Chain 1 of ", n_chains)
  ch <- .run_chain(Xc, if (st$type == "user") st$dag else st$type, sched,
                   d, alpha, verbose)
  if (!is.null(ch$start_failed)) {
    warn <- c(warn, paste0("The TL warm start failed (", ch$start_failed, "). ",
                           "Chain 1 starts at the empty graph."))
  }

  # Steps 1 to 4 for chain 1.
  proj <- .project(Xc, ch, alpha_sw)
  dec <- .decide(ch, proj, p)
  out <- .summarize(ch, proj, dec, dat$names, alpha_ball)

  # Trace and effective sample size of the score.
  trace <- data.frame(iteration = ch$iteration, log_score = ch$log_score,
                      n_edges = ch$n_edges)
  ess <- .ess_geyer(ch$log_score)
  if (ess < ess_min) {
    warn <- c(warn, sprintf(paste0("The effective sample size of the score is %.0f, ",
                                   "below ess_min = %g. Increase n_iter. Use ",
                                   "n_chains = 2 to check convergence."), ess, ess_min))
  }

  diagnostics <- list(acceptance_rate = ch$acceptance_rate,
                      invalid_rate = ch$invalid_rate, moves = ch$moves,
                      n_iter = sched$n_iter, burn_in = sched$burn_in,
                      thin = sched$thin, n_draws = length(ch$draw_dag),
                      n_unique_dags = length(ch$dags),
                      n_unique_rpdags = out$n_unique_rpdags,
                      ess_log_score = ess, best_score = ch$best_score,
                      start = ch$start_used,
                      start_dag = .positions_to_matrix(which(ch$start_dag == 1L), p,
                                                       dat$names))

  # Further chains for the convergence check.
  if (n_chains > 1L) {
    scores <- list(ch$log_score)
    edges_n <- list(ch$n_edges)
    adj <- list(out$freq$adjacent)
    dir <- list(out$freq$directed)
    info <- list(data.frame(chain = 1L, start = ch$start_used,
                            acceptance_rate = ch$acceptance_rate,
                            mean_log_score = mean(ch$log_score),
                            best_score = ch$best_score,
                            stringsAsFactors = FALSE))
    for (k in 2:n_chains) {
      if (verbose) message("Chain ", k, " of ", n_chains)
      # Seed and start of chain k.
      if (!is.null(seed)) set.seed(.chain_seed(seed, k))
      start_k <- if (k %% 2L == 0L) "empty" else "TL"
      chk <- .run_chain(Xc, start_k, sched, d, alpha, verbose)
      if (!is.null(chk$start_failed)) {
        warn <- c(warn, paste0("The TL warm start failed for chain ", k, " (",
                               chk$start_failed, "). It starts at the empty graph."))
      }
      # Edge probabilities of chain k from its RPDAGs.
      projk <- .project(Xc, chk, alpha_sw)
      fk <- .edge_freq(projk$rpdag, tabulate(chk$draw_dag, nbins = length(chk$dags)), p)
      scores[[k]] <- chk$log_score
      edges_n[[k]] <- chk$n_edges
      adj[[k]] <- fk$adjacent
      dir[[k]] <- fk$directed
      info[[k]] <- data.frame(chain = k, start = chk$start_used,
                              acceptance_rate = chk$acceptance_rate,
                              mean_log_score = mean(chk$log_score),
                              best_score = chk$best_score,
                              stringsAsFactors = FALSE)
    }

    # PSRF and the largest differences between chains.
    psrf_score <- .psrf(scores)
    psrf_edges <- .psrf(edges_n)
    spread <- function(mats, sel) {
      v <- matrix(vapply(mats, function(m) m[sel], numeric(sum(sel))), nrow = sum(sel))
      max(apply(v, 1L, function(r) diff(range(r))))
    }
    upper <- upper.tri(diag(p))
    offdiag <- row(diag(p)) != col(diag(p))
    max_adj <- spread(adj, upper)
    max_dir <- spread(dir, offdiag)

    diagnostics$psrf_log_score <- psrf_score
    diagnostics$psrf_n_edges <- psrf_edges
    diagnostics$max_diff_adjacency <- max_adj
    diagnostics$max_diff_orientation <- max_dir
    diagnostics$chains <- do.call(rbind, info)

    # Convergence warnings.
    if (isTRUE(psrf_score > 1.05) || isTRUE(psrf_edges > 1.05)) {
      warn <- c(warn, sprintf(paste0("The potential scale reduction factor is %.3f for the ",
                                     "score and %.3f for the number of edges. ",
                                     "Increase n_iter."), psrf_score, psrf_edges))
    }
    if (max(max_adj, max_dir) > 0.1) {
      warn <- c(warn, sprintf(paste0("The posterior probabilities of an adjacency or an ",
                                     "orientation differ by up to %.2f between chains. ",
                                     "Increase n_iter."), max(max_adj, max_dir)))
    }
  }

  diagnostics$runtime <- proc.time()[["elapsed"]] - t0
  diagnostics$warnings <- warn

  # The settings used.
  settings <- list(d = d, alpha = alpha, start = start, n_iter = sched$n_iter,
                   burn_in = burn_in, max_draws = max_draws, alpha_sw = alpha_sw,
                   alpha_ball = alpha_ball, center = center, n_chains = n_chains,
                   ess_min = ess_min, seed = seed, verbose = verbose,
                   n = dat$n, p = p, variables = dat$names)

  res <- structure(list(D_hat = out$D_hat, nG_hat = out$nG_hat, edges = out$edges,
                        nodes = out$nodes, graphs = out$graphs, ball = out$ball,
                        trace = trace, diagnostics = diagnostics,
                        settings = settings, call = cl),
                   class = "lewis")

  # Report the warnings.
  for (w in warn) warning(w, call. = FALSE)
  res
}
