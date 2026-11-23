test_that("all Gaussian nodes give the CPDAG", {
  # The chain 1 -> 2 -> 3 has the CPDAG 1 - 2 - 3.
  A <- matrix(0L, 3, 3)
  A[1, 2] <- 1L
  A[2, 3] <- 1L
  D <- rpdag(A)
  expect_identical(unname(D), A + t(A))
  expect_identical(colnames(D), c("V1", "V2", "V3"))
})

test_that("all non-Gaussian nodes give the DAG itself", {
  set.seed(4)
  for (t in 1:20) {
    p <- sample(3:7, 1)
    A <- random_dag(p)
    expect_identical(unname(rpdag(A, seq_len(p))), A)
  }
})

test_that("a v-structure stays directed", {
  # 1 -> 3 <- 2 is its own CPDAG.
  A <- matrix(0L, 3, 3)
  A[1, 3] <- 1L
  A[2, 3] <- 1L
  expect_identical(unname(rpdag(A)), A)
})

test_that("restored orientations propagate by Meek's rules", {
  A <- matrix(0L, 3, 3)
  A[1, 2] <- 1L
  A[2, 3] <- 1L
  # Node 1 restores 1 -> 2. Rule R1 then gives 2 -> 3.
  expect_identical(unname(rpdag(A, nG = 1)), A)
  # Node 3 restores 2 -> 3 only. The edge 1 - 2 stays undirected.
  D <- unname(rpdag(A, nG = 3))
  expect_identical(D[2, 3], 1L)
  expect_identical(D[3, 2], 0L)
  expect_identical(D[1, 2], 1L)
  expect_identical(D[2, 1], 1L)
  # Names and logical vectors select the same nodes.
  dimnames(A) <- list(c("a", "b", "c"), c("a", "b", "c"))
  expect_identical(rpdag(A, nG = "c"), rpdag(A, nG = c(FALSE, FALSE, TRUE)))
})

test_that("the example of the paper gives the stated CPDAG and RPDAG", {
  p <- 12
  edges <- rbind(c(8, 4), c(4, 3), c(3, 2), c(2, 1), c(8, 7), c(2, 6), c(1, 5),
                 c(2, 5), c(8, 9), c(6, 7), c(7, 9), c(6, 9), c(5, 12), c(5, 11),
                 c(11, 10), c(9, 10), c(12, 11), c(8, 6))
  A <- matrix(0L, p, p)
  A[edges] <- 1L
  make <- function(directed, undirected) {
    D <- matrix(0L, p, p)
    D[directed] <- 1L
    D[undirected] <- 1L
    D[undirected[, 2:1, drop = FALSE]] <- 1L
    D
  }
  # The CPDAG.
  cp <- make(rbind(c(8, 7), c(2, 6), c(8, 9), c(6, 7), c(6, 9), c(11, 10),
                   c(9, 10), c(8, 6)),
             rbind(c(8, 4), c(4, 3), c(3, 2), c(2, 1), c(1, 5), c(2, 5), c(7, 9),
                   c(5, 12), c(5, 11), c(12, 11)))
  expect_identical(unname(rpdag(A)), cp)
  # The RPDAG with the non-Gaussian nodes 1, 2, 3 and 10.
  rp <- make(rbind(c(4, 3), c(3, 2), c(2, 1), c(1, 5), c(2, 5), c(5, 12), c(5, 11),
                   c(2, 6), c(8, 6), c(8, 7), c(8, 9), c(6, 7), c(6, 9), c(11, 10),
                   c(9, 10)),
             rbind(c(8, 4), c(7, 9), c(12, 11)))
  expect_identical(unname(rpdag(A, c(1, 2, 3, 10))), rp)
})
