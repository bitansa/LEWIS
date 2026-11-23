# =============================================================================
# paper12_instance.R
#
# The 12-node DAG of Figure 7(a) of the paper, used by the study of
# Figure 1 (01_motivation). Its coefficients and error scales are drawn once
# from stable_seed("instance|paper12") and never change.
#
# Requires scenario_instances.R (stable_seed, SIGMA_LO, SIGMA_HI).
# =============================================================================

P12 <- 12L

# The 18 directed edges of Figure 7(a), parent -> child. The adjacency
# convention is A[k, j] = 1 for the edge k -> j.
PAPER_EDGES <- rbind(
  c(8, 4), c(4, 3), c(3, 2), c(2, 1), c(8, 7), c(2, 6),
  c(1, 5), c(2, 5), c(8, 9), c(6, 7), c(7, 9), c(6, 9),
  c(5, 12), c(5, 11), c(11, 10), c(9, 10), c(12, 11), c(8, 6)
)

# The first five nodes of the order in which the nodes receive t5 errors
# in Figure 1. Nodes 1, 2, 3 and 10 are the black nodes of Figure 7(a).
PAPER12_NG_STAR <- c(1L, 2L, 3L, 5L, 10L)

# The order in which the nodes receive t5 errors in Figure 1 (Appendix
# G.1, "Design of Figure 1"): 1, 2, 3, 5, 10, then the other nodes in index
# order.
NG_ORDER <- c(PAPER12_NG_STAR, 4L, 6L, 7L, 8L, 9L, 11L, 12L)

# -----------------------------------------------------------------------------
# make_paper12_view: the graph, its coefficients and its error scales.
# -----------------------------------------------------------------------------
# Returns a list in the format of instance_view(), so that it can be passed
# to draw_dataset(). The fields rho and effect_scale are placeholders. The
# error distributions of a study travel in the regime tag and the `shapes`
# argument of draw_dataset().
make_paper12_view <- function() {
  # The graph.
  dag <- matrix(0L, P12, P12)
  for (e in seq_len(nrow(PAPER_EDGES))) {
    dag[PAPER_EDGES[e, 1], PAPER_EDGES[e, 2]] <- 1L
  }
  # Save the random number state of the caller and restore it on exit.
  if (exists(".Random.seed", globalenv())) {
    old_seed <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  }
  # The signed coefficients and the error scales.
  set.seed(stable_seed("instance|paper12"))
  betas <- matrix(0, P12, P12)
  ij <- which(dag == 1L, arr.ind = TRUE)
  for (e in seq_len(nrow(ij))) {
    betas[ij[e, 1], ij[e, 2]] <- sample(c(-1, 1), 1) * runif(1, 0.25, 1.0)
  }
  sigmas <- runif(P12, SIGMA_LO, SIGMA_HI)
  list(dag = dag, betas = betas, sigmas = sigmas, p = P12,
       rho = 0, effect_scale = 1.0, key = "instance|paper12")
}

# Build the view once and check it.
VIEW12 <- make_paper12_view()
stopifnot(sum(VIEW12$dag) == 18L,
          all(PAPER12_NG_STAR %in% seq_len(P12)),
          setequal(NG_ORDER, seq_len(P12)))
