test_that("a small run returns a valid fit", {
  dat <- sim_data(300, seed = 1)
  fit <- lewis(dat$X, start = "empty", n_iter = 3000, ess_min = 0, seed = 1)
  p <- 4L
  nm <- colnames(dat$X)

  # Class and point estimate.
  expect_s3_class(fit, "lewis")
  expect_identical(dim(fit$D_hat), c(p, p))
  expect_identical(dimnames(fit$D_hat), list(nm, nm))
  expect_true(all(fit$D_hat %in% c(0L, 1L)))
  expect_true(all(diag(fit$D_hat) == 0L))
  expect_true(all(fit$nG_hat %in% nm))

  # D_hat is the RPDAG of some DAG with the non-Gaussian set nG_hat.
  idx <- match(fit$nG_hat, nm)
  is_rpdag <- vapply(all_dags(p), function(g) {
    identical(unname(rpdag(g, idx)), unname(fit$D_hat))
  }, logical(1))
  expect_true(any(is_rpdag))

  # Tables.
  expect_identical(nrow(fit$nodes), p)
  expect_true(all(fit$nodes$prob_nongaussian >= 0 & fit$nodes$prob_nongaussian <= 1))
  expect_identical(fit$nodes$in_nG_hat, nm %in% fit$nG_hat)
  n_edges <- sum(fit$D_hat == 1L & t(fit$D_hat) == 0L) +
    sum(fit$D_hat == 1L & t(fit$D_hat) == 1L) / 2
  expect_identical(nrow(fit$edges), as.integer(n_edges))
  expect_true(all(fit$edges$prob > 0 & fit$edges$prob <= fit$edges$prob_adjacent + 1e-12))
  expect_equal(sum(fit$graphs$mass), 1)
  expect_identical(fit$graphs$class[fit$graphs$shd == 0], "identical")
  expect_true(all(fit$graphs$class %in% c("identical", "upper", "lower", "rearranged")))

  # Trace, ball and diagnostics.
  s <- lewis_schedule(p, n_iter = 3000)
  expect_identical(nrow(fit$trace), s$n_draws)
  expect_identical(fit$diagnostics$n_draws, s$n_draws)
  expect_true(fit$ball$radius >= 0L)
  expect_true(fit$ball$mass_outside <= fit$ball$alpha_ball)
  expect_identical(fit$diagnostics$start, "empty")
  expect_true(fit$diagnostics$acceptance_rate > 0)
})

test_that("the same seed gives the same results", {
  X <- sim_data(200, seed = 2)$X
  a <- lewis(X, start = "empty", n_iter = 1500, ess_min = 0, seed = 7)
  b <- lewis(X, start = "empty", n_iter = 1500, ess_min = 0, seed = 7)
  expect_identical(a$D_hat, b$D_hat)
  expect_identical(a$nG_hat, b$nG_hat)
  expect_identical(a$trace, b$trace)
  expect_identical(a$graphs, b$graphs)
})

test_that("the best score is at least the score of the start graph", {
  dat <- sim_data(200, seed = 3)
  fit <- lewis(dat$X, start = dat$dag, n_iter = 500, ess_min = 0, seed = 1)
  expect_identical(fit$diagnostics$start, "user")
  # The start graph counts as visited. The best score is therefore at least
  # its score. No retained draw has a higher score than the best score.
  expect_true(fit$diagnostics$best_score >= lewis_score(dat$X, dat$dag))
  expect_true(all(fit$trace$log_score <= fit$diagnostics$best_score))
})

test_that("lewis_score() follows the complexity prior", {
  dat <- sim_data(200, seed = 4)
  s0 <- lewis_score(dat$X, dat$dag, d = 0)
  s1 <- lewis_score(dat$X, dat$dag, d = 0.03, alpha = 0.75)
  n <- nrow(dat$X)
  expect_equal(as.numeric(s0 - s1), n^0.75 * 0.03 * sum(dat$dag), tolerance = 1e-8)
  expect_equal(sum(attr(s1, "per_node")), as.numeric(s1), tolerance = 1e-10)
  expect_named(attr(s1, "per_node"), colnames(dat$X))
})

test_that("nongaussian_nodes() finds the Laplace errors", {
  dat <- sim_data(1000, seed = 5)
  ng <- nongaussian_nodes(dat$X, dat$dag)
  expect_true(all(c("x1", "x3") %in% ng))
  expect_length(attr(ng, "p_values"), 4L)
  expect_equal(attr(ng, "level"), 0.05 / 4)
})

test_that("several chains give convergence diagnostics", {
  X <- sim_data(200, seed = 6)$X
  fit <- suppressWarnings(lewis(X, start = "empty", n_iter = 1000, n_chains = 2,
                                ess_min = 0, seed = 3))
  expect_true(is.numeric(fit$diagnostics$psrf_log_score))
  expect_true(is.numeric(fit$diagnostics$max_diff_adjacency))
  expect_identical(nrow(fit$diagnostics$chains), 2L)
  expect_identical(fit$diagnostics$chains$start, c("empty", "empty"))
})

test_that("the TL warm start runs", {
  X <- sim_data(200, seed = 8)$X
  fit <- suppressWarnings(lewis(X, n_iter = 500, ess_min = 0, seed = 1))
  expect_true(fit$diagnostics$start %in% c("TL", "empty (TL failed)"))
})

test_that("print, summary and plot work", {
  X <- sim_data(200, seed = 9)$X
  fit <- lewis(X, start = "empty", n_iter = 1000, ess_min = 0, seed = 1)
  expect_output(print(fit), "LEWIS fit")
  expect_output(print(summary(fit)), "Edges of the point estimate")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  for (type in c("graph", "ball", "trace", "edges")) {
    expect_error(plot(fit, type = type), NA)
  }
})

test_that("max_draws = Inf keeps every draw after the burn-in", {
  X <- sim_data(200, seed = 10)$X
  fit <- lewis(X, start = "empty", n_iter = 1000, max_draws = Inf, ess_min = 0,
               seed = 1)
  expect_identical(fit$diagnostics$thin, 1L)
  expect_identical(nrow(fit$trace), 800L)
  expect_identical(fit$trace$iteration, 201:1000)
})

test_that("plot() rejects a second positional argument", {
  X <- sim_data(200, seed = 9)$X
  fit <- lewis(X, start = "empty", n_iter = 500, ess_min = 0, seed = 1)
  expect_error(plot(fit, "ball"), "type")
})
