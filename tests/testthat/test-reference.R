# The compiled engine against the pure R reference engine in
# tests/testthat/reference. Both must give the same draws and estimates from
# the same seed when the data are not centered.

load_reference <- function() {
  env <- new.env()
  sys.source(test_path("reference", "reference_engine.R"), envir = env)
  env
}

test_that("scores and Step 1 agree with the reference engine", {
  skip_on_cran()
  skip_if_not_installed("quantreg")
  skip_if_not_installed("igraph")
  ref <- load_reference()
  dat <- sim_data(200, seed = 11)
  X <- dat$X
  set.seed(12)
  for (t in 1:10) {
    A <- random_dag(4, 0.5)
    for (pr in c("uniform", "complexity")) {
      d <- if (pr == "uniform") 0 else 0.03
      s_ref <- ref$score_dag_cached(X, A, pr, 0.75, 0.03)
      s_pkg <- lewis_score(X, A, d = d, alpha = 0.75, center = FALSE)
      expect_equal(as.numeric(s_pkg), s_ref$total, tolerance = 1e-10)
      expect_equal(unname(attr(s_pkg, "per_node")), s_ref$per_node, tolerance = 1e-10)
    }
    ng_pkg <- match(nongaussian_nodes(X, A, center = FALSE), colnames(X))
    expect_identical(ng_pkg, as.integer(ref$estimate_nG_hat(X, A)))
  }
})

test_that("Step 2 agrees with pcalg", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  skip_if_not_installed("pcalg")
  ref <- load_reference()
  set.seed(13)
  for (t in 1:100) {
    p <- sample(3:8, 1)
    A <- random_dag(p, runif(1, 0.2, 0.6))
    nG <- sort(sample.int(p, sample(0:p, 1)))
    expect_identical(unname(rpdag(A, nG)), unname(ref$reference_rpdag(A, nG)))
  }
})

test_that("the whole procedure agrees with the reference engine", {
  skip_on_cran()
  skip_if_not_installed("quantreg")
  skip_if_not_installed("igraph")
  skip_if_not_installed("pcalg")
  ref <- load_reference()
  X <- sim_data(200, seed = 14)$X
  for (d in c(0.03, 0)) {
    set.seed(15)
    r <- ref$reference_lewis(X, n_iter = 1500,
                             prior = if (d == 0) "uniform" else "complexity", d = 0.03)
    fit <- lewis(X, d = d, start = "empty", n_iter = 1500, center = FALSE,
                 ess_min = 0, seed = 15)
    # The draws.
    expect_identical(fit$trace$n_edges, vapply(r$run$trace, function(A) sum(A), integer(1)))
    expect_equal(fit$trace$log_score, r$run$scores, tolerance = 1e-10)
    expect_identical(fit$diagnostics$moves, r$run$moves)
    expect_equal(fit$diagnostics$acceptance_rate, r$run$acceptance_rate)
    expect_identical(fit$diagnostics$n_unique_dags, r$ps$n_unique_dags)
    # The estimates.
    expect_identical(unname(fit$D_hat), unname(r$D_hat))
    expect_identical(match(fit$nG_hat, colnames(X)), r$nG_hat)
    expect_identical(fit$ball$radius, as.integer(r$radius))
  }
})
