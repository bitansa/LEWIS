# Helpers shared by the tests.

# Simulated data from the DAG x1 -> x2 -> x3, x1 -> x4, x3 -> x4. The errors
# of x1 and x3 are Laplace and the errors of x2 and x4 are Gaussian.
sim_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  e1 <- rexp(n) - rexp(n)
  e2 <- rnorm(n)
  e3 <- rexp(n) - rexp(n)
  e4 <- rnorm(n)
  x1 <- e1
  x2 <- 0.8 * x1 + e2
  x3 <- 0.8 * x2 + e3
  x4 <- 0.7 * x1 + 0.7 * x3 + e4
  dag <- matrix(0L, 4, 4)
  dag[1, 2] <- 1L
  dag[2, 3] <- 1L
  dag[1, 4] <- 1L
  dag[3, 4] <- 1L
  list(X = cbind(x1, x2, x3, x4), dag = dag, nG = c(1L, 3L))
}

# A random DAG on p nodes. Each edge of a random order is present with
# probability prob.
random_dag <- function(p, prob = 0.4) {
  ord <- sample.int(p)
  A <- matrix(0L, p, p)
  for (a in seq_len(p - 1L)) for (b in seq(a + 1L, p)) {
    if (runif(1) < prob) A[ord[a], ord[b]] <- 1L
  }
  A
}

# All DAGs on p nodes (p <= 4).
all_dags <- function(p) {
  ij <- which(matrix(1, p, p) - diag(p) == 1, arr.ind = TRUE)
  m <- nrow(ij)
  out <- list()
  for (k in 0:(2^m - 1)) {
    bits <- as.integer(intToBits(k))[seq_len(m)]
    A <- matrix(0L, p, p)
    A[ij[bits == 1L, , drop = FALSE]] <- 1L
    if (lewis_is_acyclic_cpp(A)) out[[length(out) + 1L]] <- A
  }
  out
}
