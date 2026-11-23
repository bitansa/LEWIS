test_that("the default schedule follows the paper", {
  s10 <- lewis_schedule(10)
  expect_identical(s10$n_iter, 10000L)
  expect_identical(s10$burn_in, 2000L)
  expect_identical(s10$thin, 1L)
  expect_identical(s10$n_draws, 8000L)

  s20 <- lewis_schedule(20)
  expect_identical(s20$n_iter, 10000L)
  expect_identical(s20$thin, 1L)

  s50 <- lewis_schedule(50)
  expect_identical(s50$n_iter, 62500L)
  expect_identical(s50$burn_in, 12500L)
  expect_identical(s50$thin, 7L)
  expect_identical(s50$n_draws, 7142L)

  s100 <- lewis_schedule(100)
  expect_identical(s100$n_iter, 250000L)
  expect_identical(s100$burn_in, 50000L)
  expect_identical(s100$thin, 25L)
  expect_identical(s100$n_draws, 8000L)
})

test_that("a given number of iterations is thinned to at most max_draws", {
  s <- lewis_schedule(5, n_iter = 2000)
  expect_identical(s$burn_in, 400L)
  expect_identical(s$thin, 1L)
  expect_identical(s$n_draws, 1600L)

  s <- lewis_schedule(8, n_iter = 1e7)
  expect_identical(s$thin, 1000L)
  expect_identical(s$n_draws, 8000L)

  s <- lewis_schedule(5, n_iter = 1000, burn_in = 0, max_draws = 300)
  expect_identical(s$thin, 4L)
  expect_identical(s$n_draws, 250L)
})

test_that("max_draws = Inf keeps every draw after the burn-in", {
  s <- lewis_schedule(5, n_iter = 20000, max_draws = Inf)
  expect_identical(s$burn_in, 4000L)
  expect_identical(s$thin, 1L)
  expect_identical(s$n_draws, 16000L)

  s <- lewis_schedule(50, max_draws = Inf)
  expect_identical(s$n_iter, 62500L)
  expect_identical(s$thin, 1L)
  expect_identical(s$n_draws, 50000L)
})

test_that("invalid values of max_draws are rejected", {
  expect_error(lewis_schedule(5, max_draws = 0), "or Inf")
  expect_error(lewis_schedule(5, max_draws = -Inf), "or Inf")
  expect_error(lewis_schedule(5, max_draws = 10.5), "or Inf")
  expect_error(lewis_schedule(5, max_draws = NA), "or Inf")
  expect_error(lewis_schedule(5, max_draws = "all"), "or Inf")
  expect_error(lewis_schedule(5, max_draws = c(10, 20)), "or Inf")
})
