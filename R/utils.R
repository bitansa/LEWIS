# Small helpers used by several files.

# Return a if it is not NULL, and b otherwise.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Turn a key such as "1,4,7" into an integer vector. The empty key gives an
# empty vector.
.parse_key <- function(k) {
  # The empty key is the empty set. Other keys are split at the commas.
  if (is.na(k) || identical(k, "")) {
    integer(0)
  } else {
    as.integer(strsplit(k, ",", fixed = TRUE)[[1L]])
  }
}

# Turn 1-based edge positions into a p x p 0/1 integer matrix.
.positions_to_matrix <- function(pos, p, names = NULL) {
  # Set the entries at the given positions to 1 and add the names.
  A <- matrix(0L, p, p)
  A[pos] <- 1L
  if (!is.null(names)) dimnames(A) <- list(names, names)
  A
}

# The edges of a partially directed graph given by edge positions, as keys.
# A directed edge i -> j has the key "i->j". An undirected edge i - j with
# i < j has the key "i-j". The keys are ordered by i and then by j.
.edge_keys <- function(pos, p) {
  if (length(pos) == 0L) return(character(0))
  # Row (parent) and column (child) of each nonzero entry.
  i <- (pos - 1L) %% p + 1L
  j <- (pos - 1L) %/% p + 1L
  # An entry is undirected when the reverse entry is also present.
  rev_pos <- (i - 1L) * p + j
  undirected <- rev_pos %in% pos
  # Keep each undirected edge once, with i < j.
  keep <- !undirected | (i < j)
  i <- i[keep]
  j <- j[keep]
  undirected <- undirected[keep]
  key <- ifelse(undirected, paste0(i, "-", j), paste0(i, "->", j))
  key[order(i, j)]
}

# Write edge keys with node names: "a->b" for a directed edge and "a - b"
# for an undirected edge.
.format_edge_keys <- function(keys, names) {
  if (length(keys) == 0L) return(character(0))
  # Split each key into the indices of its two nodes.
  directed <- grepl("->", keys, fixed = TRUE)
  parts <- strsplit(keys, ifelse(directed, "->", "-"), fixed = TRUE)
  # Replace the indices by the node names.
  a <- names[as.integer(vapply(parts, `[`, character(1), 1L))]
  b <- names[as.integer(vapply(parts, `[`, character(1), 2L))]
  ifelse(directed, paste0(a, "->", b), paste0(a, " - ", b))
}

# Compare the edges of a graph with those of the point estimate. A reversal
# is a directed edge whose reverse is in the point estimate. It is written
# as the orientation in the point estimate, then " => ", then the
# orientation in the graph.
.compare_edges <- function(E_hat, E_cand) {
  # Edges only in the graph and edges only in the point estimate.
  added <- setdiff(E_cand, E_hat)
  removed <- setdiff(E_hat, E_cand)
  # Look for reversed directed edges among the added edges.
  dir_added <- added[grepl("->", added, fixed = TRUE)]
  if (length(dir_added) > 0L) {
    parts <- strsplit(dir_added, "->", fixed = TRUE)
    rev_e <- vapply(parts, function(v) paste0(v[2L], "->", v[1L]), character(1))
    hit <- rev_e %in% removed
  } else {
    rev_e <- character(0)
    hit <- logical(0)
  }
  reversed_hat <- rev_e[hit]
  reversed_cand <- dir_added[hit]
  # Remove the reversed edges from the added and removed lists.
  added <- setdiff(added, reversed_cand)
  removed <- setdiff(removed, reversed_hat)
  list(added = added, removed = removed,
       reversed_hat = reversed_hat, reversed_cand = reversed_cand)
}

# The class of a graph relative to the point estimate (Step 4 of the
# paper): identical, upper, lower or rearranged.
.graph_class <- function(n_added, n_removed, n_reversed) {
  if (n_added + n_removed + n_reversed == 0L) {
    "identical"
  } else if (n_removed == 0L && n_reversed == 0L) {
    "upper"
  } else if (n_added == 0L && n_reversed == 0L) {
    "lower"
  } else {
    "rearranged"
  }
}

# Derive the seed of chain k from the seed of the first chain.
.chain_seed <- function(seed, k) {
  # seed + k - 1, kept within the range of R's integers.
  as.integer((as.numeric(seed) + k - 1) %% .Machine$integer.max)
}
