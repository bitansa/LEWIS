# A pure R implementation of LEWIS, used only by the tests.
#
# It computes the same quantities as the compiled engine of the package with
# R code: quantreg::rq.fit() for the L1 regressions, stats::shapiro.test()
# for Step 1, pcalg::dag2cpdag() and pcalg::addBgKnowledge() for Step 2, and
# an R loop for the sampler. The test test-reference.R checks that the
# package gives the same draws and the same estimates from the same seed.
#
# Source this file into an environment. It needs the packages quantreg,
# igraph and pcalg.

# ---- Settings of the paper --------------------------------------------------

LEWIS_DEFAULT_ALPHA <- 0.75
LEWIS_DEFAULT_D <- 0.03
LEWIS_DEFAULT_ALPHA_SW <- 0.05
LEWIS_ITER_BASE <- 10000L
LEWIS_ITER_P0 <- 20L
LEWIS_BURN_FRAC <- 0.2
LEWIS_MAX_RETAINED <- 8000L

# The default number of iterations, 10000 max(1, (p/20)^2).
lewis_n_iter <- function(p) {
  as.integer(round(LEWIS_ITER_BASE * max(1, (as.numeric(p) / LEWIS_ITER_P0)^2)))
}
# The burn-in: the first 20 percent of the iterations.
lewis_burn_in <- function(n_iter) as.integer(floor(n_iter * LEWIS_BURN_FRAC))
# The thinning interval that keeps at most 8000 draws.
lewis_thin <- function(n_iter, burn_in = lewis_burn_in(n_iter)) {
  max(1L, as.integer(ceiling((n_iter - burn_in) / LEWIS_MAX_RETAINED)))
}
# The default schedule for p variables.
lewis_mcmc_schedule <- function(p) {
  n_iter <- lewis_n_iter(p)
  burn_in <- lewis_burn_in(n_iter)
  thin <- lewis_thin(n_iter, burn_in)
  list(p = as.integer(p), n_iter = n_iter, burn_in = burn_in, thin = thin,
       n_retained = as.integer((n_iter - burn_in) %/% thin))
}
# The schedule with a given number of iterations.
resolve_mcmc_schedule <- function(p, n_iter = NULL, burn_in = NULL, thin = NULL) {
  s <- lewis_mcmc_schedule(p)
  if (!is.null(n_iter)) s$n_iter <- as.integer(n_iter)
  if (!is.null(burn_in)) s$burn_in <- as.integer(burn_in)
  else if (!is.null(n_iter)) s$burn_in <- lewis_burn_in(s$n_iter)
  if (!is.null(thin)) s$thin <- max(1L, as.integer(thin))
  else if (!is.null(n_iter)) s$thin <- lewis_thin(s$n_iter, s$burn_in)
  s$n_retained <- as.integer((s$n_iter - s$burn_in) %/% s$thin)
  s
}

# ---- DAG utilities ----------------------------------------------------------

# TRUE if the directed graph A has no directed cycle.
is_acyclic <- function(A) {
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
  igraph::is_dag(g)
}

# ---- L1 fit of one node -----------------------------------------------------

# L1 fit of node j on its parents. A node without parents keeps its data
# as residuals. The fallbacks are the Frisch-Newton method and ridge least
# squares.
fit_node_l1 <- function(X, j, parents) {
  n <- nrow(X)
  if (length(parents) == 0L) {
    res <- X[, j]
    return(list(b = numeric(0),
                theta = max(mean(abs(res)), 1e-12),
                residuals = res))
  }
  y <- as.numeric(X[, j])
  Xpa <- X[, parents, drop = FALSE]
  fit_coefs <- tryCatch(
    suppressWarnings(quantreg::rq.fit(Xpa, y, tau = 0.5, method = "br")$coefficients),
    error = function(e) tryCatch(
      suppressWarnings(quantreg::rq.fit(Xpa, y, tau = 0.5, method = "fn")$coefficients),
      error = function(e2) NULL
    )
  )
  if (is.null(fit_coefs) || any(!is.finite(fit_coefs))) {
    G <- crossprod(Xpa) + 1e-8 * diag(ncol(Xpa))
    fit_coefs <- as.numeric(solve(G, crossprod(Xpa, y)))
  }
  fit_coefs <- as.numeric(unname(fit_coefs))
  if (any(!is.finite(fit_coefs))) fit_coefs[!is.finite(fit_coefs)] <- 0
  res <- as.numeric(y - as.numeric(Xpa %*% fit_coefs))
  res <- res[is.finite(res)]
  theta <- if (length(res) == 0L) 1e-12 else max(mean(abs(res)), 1e-12)
  list(b = fit_coefs, theta = theta, residuals = res)
}

# ---- The Laplace-approximated score ------------------------------------------

SIGMA_B2_REFINED <- 100
IG_ALPHA_REFINED <- 1
IG_BETA_REFINED <- 1

# Gaussian kernel estimate of the residual density at zero, with
# Silverman's bandwidth.
.resid_density_at_zero <- function(resid) {
  r <- resid[is.finite(resid)]
  n <- length(r)
  if (n < 5L) return(1)
  s <- stats::sd(r)
  iqr <- stats::IQR(r) / 1.349
  scale_r <- min(s, if (is.finite(iqr) && iqr > 0) iqr else s)
  if (!is.finite(scale_r) || scale_r <= 0) return(1)
  h <- 0.9 * scale_r * n^(-1 / 5)
  if (!is.finite(h) || h <= 0) return(1)
  max(mean(stats::dnorm(r / h)) / h, 1e-8)
}

# The information matrix J_j of node j.
.build_J_node <- function(b_hat, theta_hat, X, parents, f0) {
  n <- nrow(X)
  d_j <- length(parents)
  eta <- 1 / theta_hat
  if (d_j == 0L) {
    return(matrix(theta_hat^2, 1L, 1L))
  }
  Xpa <- X[, parents, drop = FALSE]
  w_hat <- b_hat / theta_hat
  M <- crossprod(Xpa) / n
  v <- as.numeric(M %*% w_hat)
  S <- sum(w_hat * v)
  J <- matrix(0, d_j + 1L, d_j + 1L)
  J[1, 1] <- 2 * f0 / eta^3 * S + 1 / eta^2
  J[1, 2:(d_j + 1L)] <- -2 * f0 / eta^2 * v
  J[2:(d_j + 1L), 1] <- -2 * f0 / eta^2 * v
  J[2:(d_j + 1L), 2:(d_j + 1L)] <- 2 * f0 / eta * M
  J
}

# The estimate c_hat_j of the constant of the Laplace approximation.
.node_c_gamma <- function(b_hat, theta_hat, X, parents, resid = NULL,
                          sigma_b2 = SIGMA_B2_REFINED,
                          ig_alpha = IG_ALPHA_REFINED,
                          ig_beta = IG_BETA_REFINED) {
  n <- nrow(X)
  d_j <- length(parents)
  d_total <- d_j + 1L
  gauss_norm <- d_total / 2 * log(2 * pi)
  log_prior_b <- if (d_j > 0L) {
    -d_j / 2 * log(2 * pi * sigma_b2) -
      sum(b_hat^2) / (2 * sigma_b2)
  } else 0
  log_prior_theta <- ig_alpha * log(ig_beta) -
    lgamma(ig_alpha) -
    (ig_alpha + 1) * log(theta_hat) -
    ig_beta / theta_hat
  log_jacobian <- (d_j + 2L) * log(theta_hat)
  log_prior_etaw <- log_prior_b + log_prior_theta + log_jacobian
  f0 <- if (length(parents) == 0L) 1 else .resid_density_at_zero(resid)
  J <- .build_J_node(b_hat, theta_hat, X, parents, f0 = f0)
  log_det_J <- as.numeric(determinant(J, logarithm = TRUE)$modulus)
  gauss_norm + log_prior_etaw - 0.5 * log_det_J
}

# The score s_{n,j}(S) of node j with parent set S.
node_score <- function(X, j, parents, alpha = LEWIS_DEFAULT_ALPHA,
                       d = LEWIS_DEFAULT_D,
                       prior = c("uniform", "complexity")) {
  prior <- match.arg(prior)
  n <- nrow(X)
  fit <- fit_node_l1(X, j, parents)
  d_j <- length(parents)
  base <- -n * log(2) - n * log(fit$theta) - n - (1 + d_j) / 2 * log(n)
  if (prior == "complexity") base <- base - n^alpha * d * d_j
  c_gamma_j <- .node_c_gamma(fit$b, fit$theta, X, parents,
                             resid = fit$residuals)
  base + c_gamma_j
}

# The score of a DAG and the score of each node.
score_dag_cached <- function(X, dag, prior = "uniform",
                             alpha = LEWIS_DEFAULT_ALPHA,
                             d = LEWIS_DEFAULT_D) {
  p <- ncol(X)
  per <- numeric(p)
  for (j in seq_len(p)) {
    pa <- which(dag[, j] == 1L)
    per[j] <- node_score(X, j, pa, alpha, d, prior)
  }
  list(total = sum(per), per_node = per)
}

# ---- The sampler ------------------------------------------------------------

# The Metropolis-Hastings sampler over DAGs.
mcmc_dag_core <- function(start, n_iter, burn_in, thin,
                          score_init, score_update,
                          keep_trace = TRUE) {
  current <- start
  storage.mode(current) <- "integer"
  p <- nrow(current)
  stopifnot(ncol(current) == p, is_acyclic(current))
  n_iter <- as.integer(n_iter)
  burn_in <- as.integer(burn_in)
  thin <- max(1L, as.integer(thin))
  stopifnot(n_iter >= 0L, burn_in >= 0L, burn_in <= n_iter)

  n_pairs <- p * (p - 1L)
  cache <- score_init(current)
  n_edges <- sum(current)

  best <- current
  best_score <- cache$total
  n_accept <- 0L
  n_invalid <- 0L
  moves <- matrix(0L, 3L, 2L,
                  dimnames = list(c("add", "delete", "reverse"),
                                  c("proposed", "accepted")))

  n_ret <- if (n_iter > burn_in) (n_iter - burn_in) %/% thin else 0L
  trace <- if (isTRUE(keep_trace)) vector("list", n_ret) else NULL
  scores <- numeric(n_ret)
  slot <- 0L

  draw_absent_pair <- function(A, n_absent) {
    if (n_absent <= 0L) return(NULL)
    for (tries in seq_len(64L)) {
      ij <- sample.int(p, 2L)
      if (A[ij[1L], ij[2L]] == 0L && A[ij[2L], ij[1L]] == 0L) return(ij)
    }
    cand <- which(A == 0L & t(A) == 0L & !diag(p), arr.ind = TRUE)
    cand[sample.int(nrow(cand), 1L), ]
  }
  draw_edge <- function(A, n_e) {
    if (n_e <= 0L) return(NULL)
    e <- which(A == 1L, arr.ind = TRUE)
    e[sample.int(n_e, 1L), ]
  }

  for (iter in seq_len(n_iter)) {
    move <- sample.int(3L, 1L)
    moves[move, 1L] <- moves[move, 1L] + 1L
    n_absent <- n_pairs - 2L * n_edges
    new_dag <- current
    changed <- integer(0L)
    valid <- FALSE
    log_q <- 0

    if (move == 1L) {
      ij <- draw_absent_pair(current, n_absent)
      if (!is.null(ij)) {
        new_dag[ij[1L], ij[2L]] <- 1L
        if (is_acyclic(new_dag)) {
          changed <- ij[2L]
          valid <- TRUE
          log_q <- log(n_absent) - log(n_edges + 1L)
        }
      }
    } else if (move == 2L) {
      ij <- draw_edge(current, n_edges)
      if (!is.null(ij)) {
        new_dag[ij[1L], ij[2L]] <- 0L
        changed <- ij[2L]
        valid <- TRUE
        log_q <- log(n_edges) - log(n_absent + 2L)
      }
    } else {
      ij <- draw_edge(current, n_edges)
      if (!is.null(ij)) {
        new_dag[ij[1L], ij[2L]] <- 0L
        new_dag[ij[2L], ij[1L]] <- 1L
        if (is_acyclic(new_dag)) {
          changed <- c(ij[1L], ij[2L])
          valid <- TRUE
          log_q <- 0
        }
      }
    }

    if (valid) {
      new_cache <- score_update(new_dag, changed, cache)
      log_acc <- (new_cache$total - cache$total) + log_q
      if (!is.nan(log_acc) && log(runif(1)) < log_acc) {
        current <- new_dag
        cache <- new_cache
        n_edges <- if (move == 1L) n_edges + 1L else if (move == 2L) n_edges - 1L else n_edges
        n_accept <- n_accept + 1L
        moves[move, 2L] <- moves[move, 2L] + 1L
        if (cache$total > best_score) {
          best <- current
          best_score <- cache$total
        }
      }
    } else {
      n_invalid <- n_invalid + 1L
    }

    if (iter > burn_in && ((iter - burn_in) %% thin) == 0L) {
      slot <- slot + 1L
      if (isTRUE(keep_trace) && slot <= n_ret) {
        trace[[slot]] <- current
        scores[slot] <- cache$total
      }
    }
  }

  list(final = current, final_score = cache$total,
       best = best, best_score = best_score,
       trace = trace, scores = scores, n_retained = n_ret,
       acceptance_rate = if (n_iter > 0L) n_accept / n_iter else NA_real_,
       invalid_rate = if (n_iter > 0L) n_invalid / n_iter else NA_real_,
       moves = moves, n_iter = n_iter, burn_in = burn_in, thin = thin)
}

# Score functions for the sampler. Only the changed nodes are scored again.
laplace_score_fns <- function(X, prior, alpha, d) {
  list(
    init = function(dag) score_dag_cached(X, dag, prior, alpha, d),
    update = function(new_dag, changed, cache) {
      per <- cache$per_node
      for (k in changed) per[k] <- node_score(X, k, which(new_dag[, k] == 1L),
                                              alpha, d, prior)
      list(total = sum(per), per_node = per)
    })
}

# ---- Step 1: non-Gaussian nodes ---------------------------------------------

# L1 residuals of node j, with the same fallbacks as fit_node_l1().
.insample_resid <- function(X, j, pa) {
  if (length(pa) == 0L) {
    return(as.numeric(X[, j]))
  }
  Xpa <- X[, pa, drop = FALSE]
  y <- as.numeric(X[, j])
  b <- tryCatch(
    suppressWarnings(quantreg::rq.fit(Xpa, y, tau = 0.5, method = "br")$coefficients),
    error = function(e) tryCatch(
      suppressWarnings(quantreg::rq.fit(Xpa, y, tau = 0.5, method = "fn")$coefficients),
      error = function(e2) NULL))
  if (is.null(b) || any(!is.finite(b))) {
    G <- crossprod(Xpa) + 1e-8 * diag(ncol(Xpa))
    b <- as.numeric(solve(G, crossprod(Xpa, y)))
  }
  b <- as.numeric(unname(b))
  if (any(!is.finite(b))) b[!is.finite(b)] <- 0
  as.numeric(y - as.numeric(Xpa %*% b))
}

# Step 1: the nodes whose Shapiro-Wilk p-value is below alpha / p.
estimate_nG_hat <- function(X, dag, alpha = LEWIS_DEFAULT_ALPHA_SW) {
  p <- ncol(X)
  pvals <- numeric(p)
  for (j in seq_len(p)) {
    pa <- which(dag[, j] == 1L)
    res <- .insample_resid(X, j, pa)
    res <- res[is.finite(res)]
    if (length(res) > 5000L) {
      res <- res[round(seq(1, length(res), length.out = 5000L))]
    }
    if (length(res) < 3L) {
      pvals[j] <- 1
    } else {
      sw <- tryCatch(stats::shapiro.test(res), error = function(e) NULL)
      pvals[j] <- if (is.null(sw)) 1 else sw$p.value
    }
  }
  thr <- alpha / p
  sort(which(pvals < thr))
}

# ---- Step 2: the RPDAG with pcalg --------------------------------------------

# Conversion to and from the adjacency convention of pcalg, which has the
# parent in the column.
.to_pcalg_amat <- function(A) {
  B <- t(A)
  storage.mode(B) <- "integer"
  B
}
.from_pcalg_amat <- function(A) {
  B <- t(A)
  storage.mode(B) <- "integer"
  B
}

# Step 2: the CPDAG by pcalg::dag2cpdag(), then the restored orientations
# and Meek's rules by pcalg::addBgKnowledge().
reference_rpdag <- function(dag, nG_hat) {
  if (!requireNamespace("pcalg", quietly = TRUE))
    stop("pcalg is required for reference_rpdag()")
  p <- nrow(dag)
  if (p == 0L) return(matrix(0L, 0, 0))
  M <- dag
  storage.mode(M) <- "integer"
  gNEL <- as(M, "graphNEL")
  cpdag <- as(pcalg::dag2cpdag(gNEL), "matrix")
  storage.mode(cpdag) <- "integer"
  x_orient <- integer(0)
  y_orient <- integer(0)
  if (length(nG_hat) > 0L) {
    for (i in seq_len(p - 1L)) {
      for (j in seq(i + 1L, p)) {
        if (cpdag[i, j] == 1L && cpdag[j, i] == 1L) {
          if ((i %in% nG_hat) || (j %in% nG_hat)) {
            if (dag[i, j] == 1L && dag[j, i] == 0L) {
              x_orient <- c(x_orient, i)
              y_orient <- c(y_orient, j)
            } else if (dag[j, i] == 1L && dag[i, j] == 0L) {
              x_orient <- c(x_orient, j)
              y_orient <- c(y_orient, i)
            }
          }
        }
      }
    }
  }
  if (length(x_orient) == 0L) return(cpdag)
  res <- tryCatch(
    pcalg::addBgKnowledge(gInput = .to_pcalg_amat(cpdag),
                          x = x_orient, y = y_orient,
                          checkInput = FALSE),
    error = function(e) NULL)
  if (is.null(res)) stop("addBgKnowledge() returned NULL")
  out <- .from_pcalg_amat(as.matrix(res))
  dimnames(out) <- NULL
  out
}

# The status of every unordered pair of nodes.
pair_state <- function(A) {
  ut <- upper.tri(A)
  2L * as.integer(A[ut]) + as.integer(t(A)[ut])
}

# The SHD between two partially directed graphs.
reference_shd <- function(A, B) {
  if (nrow(A) != nrow(B)) stop("graphs differ in dimension")
  sum(pair_state(A) != pair_state(B))
}

# ---- Steps 1 to 3 over the draws ----------------------------------------------

# Turn a key such as "1,4" into an integer vector.
.parse_nG_key <- function(k) {
  if (is.na(k) || identical(k, "")) integer(0)
  else as.integer(strsplit(k, ",", fixed = TRUE)[[1]])
}

# Steps 1 and 2 for every distinct DAG of the draws.
per_sample_summaries <- function(X, trace, alpha_sw = LEWIS_DEFAULT_ALPHA_SW) {
  M <- length(trace)
  p <- ncol(X)
  .key_of <- function(A) paste(which(A == 1L), collapse = ",")
  dag_key <- vapply(trace, .key_of, character(1))
  uniq_keys <- unique(dag_key)
  which_u <- match(dag_key, uniq_keys)
  U <- length(uniq_keys)
  first_of <- match(seq_len(U), which_u)
  u_nG_key <- character(U)
  u_D <- vector("list", U)
  u_D_key <- character(U)
  for (u in seq_len(U)) {
    gamma_u <- trace[[first_of[u]]]
    nG_u <- estimate_nG_hat(X, gamma_u, alpha = alpha_sw)
    D_u <- reference_rpdag(gamma_u, nG_u)
    u_nG_key[u] <- if (length(nG_u) == 0L) "" else paste(nG_u, collapse = ",")
    u_D[[u]] <- D_u
    u_D_key[u] <- .key_of(D_u)
  }
  list(per_m = data.frame(m = seq_len(M),
                          dag_key = dag_key,
                          nG_hat_key = u_nG_key[which_u],
                          D_key = u_D_key[which_u],
                          stringsAsFactors = FALSE),
       D_matrices = u_D[which_u],
       n_unique_dags = U)
}

# The number of edges. An undirected edge counts once.
.n_edges_pdag <- function(D) {
  und <- sum(D == 1L & t(D) == 1L) / 2L
  dir <- sum(D == 1L & t(D) == 0L)
  as.numeric(und + dir)
}

# Step 3: the SHD medoid and its most frequent non-Gaussian set.
choose_shd_medoid <- function(per_m, D_matrices) {
  M <- nrow(per_m)
  unique_Dk <- unique(per_m$D_key)
  idx <- match(per_m$D_key, unique_Dk)
  counts <- as.numeric(tabulate(idx, nbins = length(unique_Dk)))
  rep_of <- match(seq_along(unique_Dk), idx)
  U <- length(unique_Dk)
  D_unique <- lapply(rep_of, function(i) D_matrices[[i]])
  L <- as.integer(nrow(D_unique[[1L]]) * (nrow(D_unique[[1L]]) - 1L) / 2L)
  S <- vapply(D_unique, pair_state, integer(L))
  if (is.null(dim(S))) S <- matrix(S, ncol = U)
  L <- nrow(S)
  W <- matrix(0, 4L, L)
  for (k in 0:3) W[k + 1L, ] <- as.numeric((S == k) %*% counts)
  idx2 <- (rep(seq_len(L), times = U) - 1L) * 4L + as.vector(S) + 1L
  agree <- colSums(matrix(W[idx2], nrow = L))
  mean_shd <- (L * M - agree) / M
  edge_cnts <- vapply(D_unique, .n_edges_pdag, numeric(1))
  best <- order(mean_shd, edge_cnts)[1L]
  bDkey <- unique_Dk[best]
  D_hat <- D_matrices[[rep_of[best]]]
  matches <- per_m$nG_hat_key[per_m$D_key == bDkey]
  tab <- table(matches)
  cand <- names(tab)[tab == max(tab)]
  if (length(cand) > 1L) {
    sizes <- vapply(cand, function(k) length(.parse_nG_key(k)), integer(1))
    cand <- cand[sizes == min(sizes)]
    cand <- cand[order(match(cand, matches))]
  }
  list(D_hat = D_hat, nG_hat = .parse_nG_key(cand[1L]), mean_shd = mean_shd[best])
}

# Step 3: the radius of the credible ball.
credible_ball_radius <- function(D_hat, D_matrices, alpha) {
  M <- length(D_matrices)
  dists <- vapply(D_matrices, function(D) reference_shd(D_hat, D), numeric(1))
  k <- ceiling((1 - alpha) * M)
  list(eps_alpha = sort(dists)[k], distances = dists, k = k, M = M)
}

# ---- The whole procedure from the empty graph ---------------------------------

# The whole procedure, started at the empty graph.
reference_lewis <- function(X, n_iter, prior = "complexity", alpha = LEWIS_DEFAULT_ALPHA,
                            d = LEWIS_DEFAULT_D, alpha_sw = LEWIS_DEFAULT_ALPHA_SW,
                            alpha_ball = 0.05) {
  p <- ncol(X)
  sched <- resolve_mcmc_schedule(p, n_iter = n_iter)
  fns <- laplace_score_fns(X, prior, alpha, d)
  run <- mcmc_dag_core(matrix(0L, p, p), n_iter = sched$n_iter, burn_in = sched$burn_in,
                       thin = sched$thin, score_init = fns$init,
                       score_update = fns$update, keep_trace = TRUE)
  ps <- per_sample_summaries(X, run$trace, alpha_sw = alpha_sw)
  pair <- choose_shd_medoid(ps$per_m, ps$D_matrices)
  ball <- credible_ball_radius(pair$D_hat, ps$D_matrices, alpha_ball)
  list(run = run, ps = ps, D_hat = pair$D_hat, nG_hat = pair$nG_hat,
       radius = ball$eps_alpha)
}
