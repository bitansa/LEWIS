# Post-processing of the draws (Steps 1 to 4 of the paper).
#
# The chain returns each distinct DAG once, as edge positions. For each
# retained draw it returns the index of its DAG. Steps 1 and 2 depend only
# on the DAG. They are therefore computed once per distinct DAG. Steps 3
# and 4 work with the distinct RPDAGs and their counts.

#' Steps 1 and 2 for every distinct DAG of a chain
#' @param X Numeric data matrix (already centered if requested).
#' @param ch Output of the compiled sampler.
#' @param alpha_sw Shapiro-Wilk level before the Bonferroni division.
#' @return A list with `nongaussian` (one index vector per distinct DAG)
#'   and `rpdag` (one edge position vector per distinct DAG).
#' @noRd
.project <- function(X, ch, alpha_sw) {
  # Steps 1 and 2 in compiled code, once per distinct DAG.
  lewis_project_cpp(X, ch$dags, alpha_sw)
}

#' Step 3: the point estimate and its non-Gaussian set
#'
#' Groups the draws by their RPDAG and by their non-Gaussian set, in the
#' order in which the chain visits them, and calls the compiled SHD medoid.
#' @return A list with the distinct RPDAGs (`D_unique`), the RPDAG index of
#'   each draw (`idx`), the non-Gaussian key of each draw (`nG_key`), the
#'   index of the point estimate (`best`), its non-Gaussian set (`nG_hat`),
#'   and the mean SHD and edge count of every distinct RPDAG.
#' @noRd
.decide <- function(ch, proj, p) {
  # Keys of the RPDAG and of the non-Gaussian set of each distinct DAG.
  u_D_key <- vapply(proj$rpdag, paste, character(1), collapse = ",")
  u_nG_key <- vapply(proj$nongaussian, paste, character(1), collapse = ",")

  # Keys of each draw.
  D_key <- u_D_key[ch$draw_dag]
  nG_key <- u_nG_key[ch$draw_dag]

  # Distinct RPDAGs in the order of the first visit.
  unique_Dk <- unique(D_key)
  idx <- match(D_key, unique_Dk)
  rep_of <- match(seq_along(unique_Dk), idx)
  D_unique <- proj$rpdag[ch$draw_dag[rep_of]]

  # Distinct non-Gaussian sets in the order of the first visit.
  unique_nGk <- unique(nG_key)
  nG_idx <- match(nG_key, unique_nGk)
  nG_size <- vapply(unique_nGk, function(k) length(.parse_key(k)), integer(1),
                    USE.NAMES = FALSE)

  # SHD medoid and the most frequent non-Gaussian set among its draws.
  med <- lewis_medoid_cpp(idx, D_unique, nG_idx, nG_size, p)
  list(D_unique = D_unique, idx = idx, nG_key = nG_key, best = med$best_D,
       nG_hat = .parse_key(unique_nGk[med$best_nG]),
       mean_shd = med$mean_shd, edge_counts = med$edge_counts)
}

#' Edge frequencies over the draws
#'
#' @param graphs List of edge position vectors.
#' @param counts Number of draws of each graph.
#' @param p Number of nodes.
#' @return A list of p x p matrices: `directed` (proportion of draws with
#'   i -> j), `undirected` (proportion with i - j) and `adjacent`
#'   (proportion with i and j adjacent).
#' @noRd
.edge_freq <- function(graphs, counts, p) {
  # Weighted counts of the directed and undirected edges.
  M <- sum(counts)
  f <- lewis_edge_freq_cpp(graphs, as.numeric(counts), p)
  # Proportions. A pair is adjacent if it has an edge of any kind.
  dir <- f$directed / M
  und <- f$undirected / M
  list(directed = dir, undirected = und, adjacent = dir + t(dir) + und)
}

#' Credible ball
#'
#' The radius is the smallest SHD from the point estimate for which the
#' ball holds at least a fraction 1 - alpha_ball of the draws. It is the
#' ceiling((1 - alpha_ball) M)-th smallest distance.
#' @param dists SHD from the point estimate to each projected draw.
#' @param alpha_ball One minus the credible level.
#' @return A list with `radius`, `alpha_ball` and `mass_outside`.
#' @noRd
.credible_ball <- function(dists, alpha_ball) {
  # The k-th smallest distance, with k = ceiling((1 - alpha_ball) M).
  M <- length(dists)
  k <- ceiling((1 - alpha_ball) * M)
  radius <- sort(dists)[k]
  # The radius and the proportion of draws outside the ball.
  list(radius = as.integer(radius), alpha_ball = alpha_ball,
       mass_outside = mean(dists > radius))
}

#' The most frequent key, with ties broken in favor of the smaller set and
#' then of the key that appears first.
#' @noRd
.mode_key <- function(keys) {
  # Count the keys in the order of their first appearance.
  tab <- table(factor(keys, levels = unique(keys)))
  cand <- names(tab)[tab == max(tab)]
  # Break ties in favor of the smaller set, then of the first key.
  if (length(cand) > 1L) {
    sizes <- vapply(cand, function(k) length(.parse_key(k)), integer(1))
    cand <- cand[sizes == min(sizes)]
  }
  cand[1L]
}

#' Steps 3 and 4: point estimate, credible ball, edge and node
#' probabilities, and the table of distinct RPDAGs
#' @noRd
.summarize <- function(ch, proj, dec, names, alpha_ball) {
  p <- length(names)
  M <- length(ch$draw_dag)
  U <- length(dec$D_unique)

  # The point estimate.
  D_hat_pos <- dec$D_unique[[dec$best]]
  D_hat <- .positions_to_matrix(D_hat_pos, p, names)
  nG_hat <- dec$nG_hat

  # Credible ball: the smallest radius that holds a fraction 1 - alpha_ball
  # of the projected draws.
  shd_u <- lewis_shd_list_cpp(D_hat_pos, dec$D_unique, p)
  ball <- .credible_ball(as.numeric(shd_u[dec$idx]), alpha_ball)
  radius <- ball$radius

  # Edge and adjacency probabilities.
  counts_D <- tabulate(dec$idx, nbins = U)
  freq <- .edge_freq(dec$D_unique, counts_D, p)

  # The edges of the point estimate, one row per adjacent pair.
  ij <- which(upper.tri(D_hat) & (D_hat == 1L | t(D_hat) == 1L), arr.ind = TRUE)
  ij <- ij[order(ij[, 1L], ij[, 2L]), , drop = FALSE]
  i <- unname(ij[, 1L])
  j <- unname(ij[, 2L])
  und <- D_hat[cbind(i, j)] == 1L & D_hat[cbind(j, i)] == 1L
  fwd <- D_hat[cbind(i, j)] == 1L & !und
  from <- ifelse(fwd | und, i, j)
  to <- ifelse(fwd | und, j, i)
  prob <- ifelse(und, freq$undirected[cbind(i, j)], freq$directed[cbind(from, to)])
  edges <- data.frame(from = names[from], to = names[to],
                      type = ifelse(und, "undirected", "directed"),
                      prob = as.numeric(prob),
                      prob_adjacent = as.numeric(freq$adjacent[cbind(i, j)]),
                      stringsAsFactors = FALSE)

  # Proportion of draws in which each node is declared non-Gaussian.
  counts_dag <- tabulate(ch$draw_dag, nbins = length(ch$dags))
  ng_count <- numeric(p)
  for (u in seq_along(proj$nongaussian)) {
    v <- proj$nongaussian[[u]]
    if (length(v) > 0L) ng_count[v] <- ng_count[v] + counts_dag[u]
  }
  nodes <- data.frame(node = names, prob_nongaussian = ng_count / M,
                      in_nG_hat = seq_len(p) %in% nG_hat,
                      stringsAsFactors = FALSE)

  # Step 4: the distinct RPDAGs in decreasing order of posterior mass. Ties
  # go to the graph visited first.
  ord <- order(-counts_D, seq_len(U))
  E_hat <- .edge_keys(D_hat_pos, p)
  draws_of <- split(seq_len(M), factor(dec$idx, levels = seq_len(U)))
  rows <- lapply(ord, function(u) {
    # Edges added, removed and reversed relative to the point estimate.
    ed <- .compare_edges(E_hat, .edge_keys(dec$D_unique[[u]], p))
    rev_txt <- if (length(ed$reversed_hat) > 0L) {
      paste(.format_edge_keys(ed$reversed_hat, names), "=>",
            .format_edge_keys(ed$reversed_cand, names))
    } else {
      character(0)
    }
    # The most frequent non-Gaussian set among the draws of this graph.
    ng <- .parse_key(.mode_key(dec$nG_key[draws_of[[u]]]))
    list(class = .graph_class(length(ed$added), length(ed$removed),
                              length(ed$reversed_hat)),
         added = paste(.format_edge_keys(ed$added, names), collapse = " | "),
         removed = paste(.format_edge_keys(ed$removed, names), collapse = " | "),
         reversed = paste(rev_txt, collapse = " | "),
         nongaussian = paste(names[ng], collapse = ", "))
  })
  pick <- function(field) vapply(rows, `[[`, character(1), field)
  graphs <- data.frame(rank = seq_along(ord), mass = counts_D[ord] / M,
                       class = pick("class"), shd = shd_u[ord],
                       n_edges = as.integer(dec$edge_counts[ord]),
                       added = pick("added"), removed = pick("removed"),
                       reversed = pick("reversed"), nongaussian = pick("nongaussian"),
                       in_ball = shd_u[ord] <= radius,
                       stringsAsFactors = FALSE)

  list(D_hat = D_hat, nG_hat = names[nG_hat], nG_hat_idx = nG_hat,
       edges = edges, nodes = nodes, graphs = graphs, ball = ball,
       freq = freq, n_unique_rpdags = U)
}
