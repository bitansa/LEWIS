test_that("the SHD is zero for identical graphs and symmetric", {
  set.seed(2)
  for (t in 1:20) {
    p <- sample(3:7, 1)
    A <- random_dag(p)
    B <- random_dag(p)
    expect_identical(shd(A, A), 0L)
    expect_identical(shd(A, B), shd(B, A))
    expect_identical(shd(rpdag(A), rpdag(A)), 0L)
  }
})

test_that("the SHD counts the four states of a pair", {
  none <- matrix(0, 2, 2)
  fwd <- matrix(c(0, 0, 1, 0), 2, 2)   # 1 -> 2
  bwd <- t(fwd)                        # 2 -> 1
  und <- fwd + bwd                     # 1 - 2
  states <- list(none, fwd, bwd, und)
  for (a in seq_along(states)) for (b in seq_along(states)) {
    expect_identical(shd(states[[a]], states[[b]]), if (a == b) 0L else 1L)
  }
})

test_that("the SHD counts each differing pair once", {
  A <- matrix(0, 4, 4)
  A[1, 2] <- 1
  A[2, 3] <- 1
  B <- A
  B[1, 2] <- 0
  B[2, 1] <- 1   # reversed
  B[3, 4] <- 1   # added
  B[4, 3] <- 1   # as an undirected edge
  expect_identical(shd(A, B), 2L)
})
