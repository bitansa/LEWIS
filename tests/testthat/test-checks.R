test_that("lewis() rejects invalid data", {
  X <- sim_data(100)$X
  expect_error(lewis(X[, 1, drop = FALSE]), "at least 2 columns")
  expect_error(lewis(X[1:2, ]), "at least 3 rows")
  expect_error(lewis(X[1:4, ]), "n > p")
  expect_error(lewis("a"), "numeric matrix")
  expect_error(lewis(data.frame(a = letters[1:10], b = 1:10)), "must be numeric")

  Xna <- X
  Xna[1, 1] <- NA
  expect_error(lewis(Xna), "missing or infinite")
  Xinf <- X
  Xinf[2, 2] <- Inf
  expect_error(lewis(Xinf), "missing or infinite")

  expect_error(lewis(cbind(X, 1)), "constant")
  expect_error(lewis(cbind(X, X[, 1])), "duplicate")
  expect_error(lewis(cbind(X, X[, 1] + X[, 2])), "linearly dependent")

  Xdup <- X
  colnames(Xdup) <- c("a", "a", "b", "c")
  expect_error(lewis(Xdup), "unique")
})

test_that("lewis() rejects invalid arguments", {
  X <- sim_data(100)$X
  expect_error(lewis(X, d = -1), "zero or positive")
  expect_error(lewis(X, d = NA), "finite number")
  expect_error(lewis(X, alpha = 0.5), "strictly between 0.5 and 1")
  expect_error(lewis(X, alpha = 1), "strictly between 0.5 and 1")
  expect_error(lewis(X, n_iter = 0), "whole number")
  expect_error(lewis(X, n_iter = 10.5), "whole number")
  expect_error(lewis(X, burn_in = 1), "fraction")
  expect_error(lewis(X, max_draws = 0), "whole number")
  expect_error(lewis(X, alpha_sw = 0), "strictly between 0 and 1")
  expect_error(lewis(X, alpha_ball = 1), "strictly between 0 and 1")
  expect_error(lewis(X, n_chains = 0), "whole number")
  expect_error(lewis(X, ess_min = -1), "zero or positive")
  expect_error(lewis(X, seed = "a"), "seed")
  expect_error(lewis(X, seed = 1.5), "seed")
  expect_error(lewis(X, center = NA), "TRUE or FALSE")
  expect_error(lewis(X, verbose = "yes"), "TRUE or FALSE")
})

test_that("lewis() rejects invalid start graphs", {
  X <- sim_data(100)$X
  expect_error(lewis(X, start = "foo"), "start must be")
  expect_error(lewis(X, start = matrix(0, 3, 3)), "4 x 4")
  bad <- matrix(0, 4, 4)
  bad[1, 2] <- 2
  expect_error(lewis(X, start = bad), "only 0 and 1")
  bad <- matrix(0, 4, 4)
  diag(bad) <- 1
  expect_error(lewis(X, start = bad), "zero diagonal")
  bad <- matrix(0, 4, 4)
  bad[1, 2] <- 1
  bad[2, 1] <- 1
  expect_error(lewis(X, start = bad), "cycle")
  bad <- matrix(0, 4, 4)
  bad[1, 2] <- NA
  expect_error(lewis(X, start = bad), "only 0 and 1")
})

test_that("alpha is ignored with a message when d = 0", {
  X <- sim_data(100)$X
  expect_message(lewis(X, d = 0, alpha = 2, start = "empty", n_iter = 200,
                       ess_min = 0, seed = 1), "not used")
})

test_that("the building blocks check their input", {
  dat <- sim_data(100)
  expect_error(lewis_score(dat$X, matrix(0, 3, 3)), "4 x 4")
  expect_error(nongaussian_nodes(dat$X, dat$dag, alpha_sw = 2), "strictly between")
  cyc <- matrix(0, 3, 3)
  cyc[1, 2] <- 1
  cyc[2, 3] <- 1
  cyc[3, 1] <- 1
  expect_error(rpdag(cyc), "cycle")
  expect_error(rpdag(dat$dag, nG = 7), "between 1 and 4")
  expect_error(rpdag(dat$dag, nG = "zz"), "Unknown node names")
  expect_error(shd(matrix(0, 2, 2), matrix(0, 3, 3)), "same number")
  expect_error(shd(matrix(2, 2, 2), matrix(0, 2, 2)), "only 0 and 1")
  expect_error(lewis_schedule(1), "whole number")
})
