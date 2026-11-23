# =============================================================================
# scenario_instances.R
#
# Data generation shared by the simulation studies (comparative study,
# calibration study, accuracy of the Laplace approximation).
#
# FIXED INSTANCES. For each graph type and each number of nodes p there is
# exactly one instance: one DAG, one set of signed edge coefficients, one
# vector of error scales, and one ordering of the nodes that decides which
# nodes are non-Gaussian. The instance is drawn once from a seed that
# depends only on the graph type and p. It is the same for every dataset,
# sample size, fraction of non-Gaussian errors, effect-size scale and
# method. The datasets of a setting differ only in the errors.
#
# Settings are views of an instance:
#   * effect-size scale s: the coefficients are s times the fixed signed
#     coefficients of the instance.
#   * fraction rho of non-Gaussian errors: the non-Gaussian nodes are the
#     first round(rho * p) nodes of the fixed ordering. The sets are nested
#     in rho.
#   * error distributions: the non-Gaussian distributions are assigned in
#     rotation along the fixed ordering (Laplace, t3, t5, contaminated
#     Gaussian). A node keeps its distribution at every rho at which it is
#     non-Gaussian.
#
# SEEDS. Every seed is computed by stable_seed() from a readable string:
#
#     instance : "instance|<graph type>|p<p>"
#     dataset  : "errors|<instance key>|rho<rho>|s<s>|n<n>|<regime tag>|rep<r>"
#
# The regime tag is part of the seed. It must never be renamed, or the
# datasets change.
#
# ERRORS. Node j has the error sigma_j * Z_j. The scale sigma_j is drawn once
# per instance from Unif[0.5, 1.5]. Z_j has mean zero and unit variance:
#
#     "gaussian"  N(0, 1)
#     "laplace"   Laplace scaled to unit variance
#     "t<df>"     Student's t with df > 2 degrees of freedom, scaled to unit
#                 variance
#     "contam"    0.9 N(0, 1) + 0.1 N(0, 9), scaled to unit variance
#
# All non-Gaussian distributions are scale mixtures of Gaussians, as the
# model of the paper requires.
#
# EFFECT SIZES. The edge k -> j has the coefficient
#     beta_kj = s * sign_kj * u_kj,  sign_kj uniform on {-1, +1},
#     u_kj ~ Unif[0.25, 1].
#
# FUNCTIONS
#     stable_seed(...)                  string to integer seed
#     sample_dag_er(p, mean_degree)     Erdos-Renyi DAG
#     sample_dag_sf(p, mean_degree)     scale-free DAG
#     build_instance(graph_type, p)     the fixed instance
#     instance_view(inst, rho, s)       the view of a setting
#     draw_dataset(view, n, rep, regime_tag, shapes)
#                                       one dataset
#
# Dependency: igraph (topological order).
# =============================================================================

suppressPackageStartupMessages(library(igraph))

# -----------------------------------------------------------------------------
# stable_seed: a seed computed from a string.
# -----------------------------------------------------------------------------
# A polynomial hash of the UTF-8 code points of the pasted arguments, modulo
# the prime 2^31 - 1. All intermediate values stay far below 2^53, so the
# arithmetic is exact and the seed is the same on every platform.
stable_seed <- function(...) {
  s <- paste(..., sep = "|")
  h <- 0
  for (code in utf8ToInt(s)) h <- (h * 131 + code) %% 2147483647
  as.integer(h)
}

# -----------------------------------------------------------------------------
# Graph samplers.
# -----------------------------------------------------------------------------

# Erdos-Renyi DAG with expected degree mean_degree. The nodes are put in a
# random order, and each pair is joined with probability
# mean_degree / (p - 1) by an edge from the earlier to the later node.
sample_dag_er <- function(p, mean_degree = 2) {
  prob_edge <- min(1, mean_degree / (p - 1))
  ord <- sample.int(p)
  A <- matrix(0L, p, p)
  for (i in seq_len(p - 1)) for (j in (i + 1):p) {
    if (runif(1) < prob_edge) A[ord[i], ord[j]] <- 1L
  }
  A
}

# Scale-free DAG grown by preferential attachment along a random order.
# Each new node receives m = round(mean_degree / 2) parents among the
# earlier nodes, chosen with probability proportional to 1 plus their
# current degree. The first m + 1 nodes form a complete DAG.
sample_dag_sf <- function(p, mean_degree = 4) {
  m <- max(1L, round(mean_degree / 2))
  if (m >= p) m <- p - 1L
  ord <- sample.int(p)
  A <- matrix(0L, p, p)
  init_n <- min(m + 1L, p)
  if (init_n >= 2) {
    for (i in seq_len(init_n - 1L)) for (j in (i + 1L):init_n) {
      A[ord[i], ord[j]] <- 1L
    }
  }
  if (p > init_n) {
    for (jpos in (init_n + 1L):p) {
      existing <- ord[seq_len(jpos - 1L)]
      deg <- rowSums(A[existing, , drop = FALSE]) +
        colSums(A[, existing, drop = FALSE])
      w <- deg + 1
      parents <- sample(existing, size = min(m, length(existing)),
                        replace = FALSE, prob = w / sum(w))
      A[parents, ord[jpos]] <- 1L
    }
  }
  A
}

# -----------------------------------------------------------------------------
# Error distributions with unit variance.
# -----------------------------------------------------------------------------

# Laplace draws with the given scale.
.rlaplace <- function(n, scale = 1) scale * (rexp(n) - rexp(n))

# n draws with unit variance from the named distribution.
.draw_shape <- function(shape, n) {
  if (shape == "gaussian") return(rnorm(n))
  if (shape == "laplace") return(.rlaplace(n, 1 / sqrt(2)))
  if (shape == "contam") {
    # 0.9 N(0, 1) + 0.1 N(0, 9) has variance 1.8.
    lam <- ifelse(runif(n) < 0.1, 3, 1)
    return(rnorm(n) * lam / sqrt(1.8))
  }
  if (grepl("^t[0-9]+(\\.[0-9]+)?$", shape)) {
    # Student's t with df degrees of freedom has variance df / (df - 2).
    df <- as.numeric(sub("^t", "", shape))
    stopifnot(df > 2)
    return(rt(n, df = df) / sqrt(df / (df - 2)))
  }
  stop("Unknown error distribution: ", shape)
}

# The non-Gaussian distributions, assigned in rotation along the fixed
# ordering of the nodes.
NG_SHAPE_POOL <- c("laplace", "t3", "t5", "contam")

# Range of the error scales: sigma_j ~ Unif[SIGMA_LO, SIGMA_HI].
SIGMA_LO <- 0.5
SIGMA_HI <- 1.5

# -----------------------------------------------------------------------------
# build_instance: the fixed instance of a graph type and p.
# -----------------------------------------------------------------------------
# Returns a list with
#     dag            p x p adjacency matrix (dag[k, j] = 1 means k -> j)
#     base_betas     p x p signed coefficients (multiply by s for a setting)
#     sigmas         the p error scales
#     ng_perm        the fixed ordering of the nodes
#     shape_by_perm  the non-Gaussian distribution of each position of the
#                    ordering
#     graph_type, p, key, seed
#
# The random number state of the caller is saved and restored, so building
# an instance does not change the random numbers of the caller. The result
# is cached within the R session.
.instance_memo <- new.env(parent = emptyenv())

build_instance <- function(graph_type = c("ER", "SF"), p) {
  graph_type <- match.arg(graph_type)
  key <- sprintf("instance|%s|p%d", graph_type, p)
  if (!is.null(.instance_memo[[key]])) return(.instance_memo[[key]])

  # Save the random number state of the caller.
  if (exists(".Random.seed", globalenv())) {
    old_seed <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  }
  set.seed(stable_seed(key))

  # The graph.
  dag <- if (graph_type == "ER") sample_dag_er(p, mean_degree = 2)
         else sample_dag_sf(p, mean_degree = 4)

  # The signed coefficients.
  base_betas <- matrix(0, p, p)
  ij <- which(dag == 1, arr.ind = TRUE)
  for (e in seq_len(nrow(ij))) {
    sgn <- sample(c(-1, 1), 1)
    base_betas[ij[e, 1], ij[e, 2]] <- sgn * runif(1, 0.25, 1.0)
  }

  # The error scales, the ordering of the nodes and the distributions.
  sigmas <- runif(p, SIGMA_LO, SIGMA_HI)
  ng_perm <- sample.int(p)
  shape_by_perm <- NG_SHAPE_POOL[((seq_len(p) - 1L) %%
                                  length(NG_SHAPE_POOL)) + 1L]

  inst <- list(dag = dag, base_betas = base_betas, sigmas = sigmas,
               ng_perm = ng_perm, shape_by_perm = shape_by_perm,
               graph_type = graph_type, p = p, key = key,
               seed = stable_seed(key))
  .instance_memo[[key]] <- inst
  inst
}

# -----------------------------------------------------------------------------
# instance_view: the view of a setting (rho, s) of an instance.
# -----------------------------------------------------------------------------
# Adds to the instance
#     rho, effect_scale
#     betas    effect_scale * base_betas
#     nG_star  the sorted first round(rho * p) nodes of the ordering
#     shapes   the error distribution of each node ("gaussian" outside
#              nG_star)
instance_view <- function(inst, rho, effect_scale = 1.0) {
  p <- inst$p
  n_ng <- max(0L, round(rho * p))
  sel <- if (n_ng == 0L) integer(0) else inst$ng_perm[seq_len(n_ng)]
  nG_star <- sort(sel)

  # The distribution of each node.
  shapes <- rep("gaussian", p)
  if (n_ng > 0L) {
    shapes[inst$ng_perm[seq_len(n_ng)]] <- inst$shape_by_perm[seq_len(n_ng)]
  }

  out <- inst
  out$rho <- rho
  out$effect_scale <- effect_scale
  out$betas <- effect_scale * inst$base_betas
  out$nG_star <- nG_star
  out$shapes <- shapes
  out
}

# -----------------------------------------------------------------------------
# draw_dataset: one dataset of a setting.
# -----------------------------------------------------------------------------
# Generates n observations of the linear SEM of the view,
#
#     X_j = sum_{k in pa(j)} beta_kj X_k + sigma_j Z_j,
#
# visiting the nodes in topological order. The seed is computed from the
# instance, rho, s, n, the regime tag and the dataset number. So a dataset
# is the same for every method, machine and run.
#
# `shapes` gives the error distribution of each node. The default is the
# assignment of the view. A study that sets the distributions itself passes
# `shapes` together with its own regime tag.
#
# With return_errors = TRUE the n x p matrix of the unit-variance draws Z is
# attached as attr(X, "Z"). The random numbers drawn are the same either way.
draw_dataset <- function(view, n, rep,
                         regime_tag = "mixture",
                         shapes = view$shapes,
                         return_errors = FALSE) {
  p <- view$p
  stopifnot(length(shapes) == p, !is.null(view$betas))

  # Seed of this dataset.
  set.seed(stable_seed("errors", view$key,
                       paste0("rho", format(view$rho, trim = TRUE)),
                       paste0("s", format(view$effect_scale, trim = TRUE)),
                       paste0("n", n), regime_tag, paste0("rep", rep)))

  # Visit the nodes in topological order.
  g <- graph_from_adjacency_matrix(view$dag, mode = "directed")
  order_top <- as.integer(topo_sort(g))
  X <- matrix(0, n, p)
  Z <- if (isTRUE(return_errors)) matrix(0, n, p) else NULL
  for (j in order_top) {
    # Draw the errors of node j and add the effects of its parents.
    z_j <- .draw_shape(shapes[j], n)
    eps_j <- view$sigmas[j] * z_j
    if (!is.null(Z)) Z[, j] <- z_j
    pa <- which(view$dag[, j] == 1)
    X[, j] <- if (length(pa) > 0)
      as.numeric(X[, pa, drop = FALSE] %*% view$betas[pa, j]) + eps_j
    else eps_j
  }
  if (!is.null(Z)) attr(X, "Z") <- Z
  X
}
