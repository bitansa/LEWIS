# =============================================================================
# real_data_lewis.R
#
# LEWIS (complexity) on the patient-by-gene matrix written by
# prepare_tcga.R, for one network (RD_NETWORK) and one chain length
# (RD_N_ITER). Section 7.4 and Appendix I.
#
#   1. The reported chain (seed RD_SEED): the point estimate D_hat, the
#      estimated non-Gaussian genes, the credible ball, the posterior
#      probabilities of every pair of genes, and the distinct RPDAGs with
#      their posterior mass and class (Step 4 of the post-processing).
#   2. The Shapiro-Wilk p-values of Table 2.
#   3. The stability check: the whole procedure with the ten seeds
#      RD_STABILITY_SEEDS.
#   4. The convergence check: four additional chains with the seeds
#      RD_CONV_SEEDS, two from the TL estimate and two from the empty graph.
#
# Outputs, in <RD_OUT_DIR>/tables/ (<tag> = RD_TAG, for example
# PRAD_ar8_it10000k):
#   point_estimate_<tag>.csv   the edges of D_hat with their posterior
#                              probability, the agreement with the known
#                              regulation and the recurrence over the seeds
#   pair_posterior_<tag>.csv   every pair of genes: P(i->j), P(j->i),
#                              P(i - j), P(absent)
#   node_posterior_<tag>.csv   every gene: class, Shapiro-Wilk p-value,
#                              level, proportion of draws in which it is
#                              declared non-Gaussian
#   posterior_rpdag_<tag>.csv  every distinct RPDAG of the draws with its
#                              mass, SHD to D_hat, class and the edges that
#                              differ from D_hat
#   uq_<tag>/                  the 15 most probable RPDAGs, a summary of the
#                              chain and the gene names, read by
#                              plot_real_data_figures.R
#   stability_<tag>.csv        the point estimates of the ten other seeds
#   convergence_<tag>.csv      the four additional chains
#   summary_<tag>.csv          one row: size, Gaussian and non-Gaussian
#                              genes, agreement counts, effective sample size
#
# Run:  RD_NETWORK=PRAD_ar8 RD_N_ITER=10000000 Rscript real_data_lewis.R
#       (run prepare_tcga.R first with the same RD_NETWORK and RD_OUT_DIR)
# =============================================================================

suppressPackageStartupMessages(library(LEWIS))

# ---- Folder of this script and settings -------------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
source(file.path(.here, "real_data_config.R"))

# ---- Data -----------------------------------------------------------------------
expr_file <- file.path(RD_TAB_DIR, sprintf("expression_%s.csv", RD_NETWORK))
if (!file.exists(expr_file)) stop("Run prepare_tcga.R first. Missing file: ", expr_file)
dat <- read.csv(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(dat[, -1, drop = FALSE])
storage.mode(X) <- "double"
genes <- colnames(X)
map <- read.csv(file.path(RD_TAB_DIR, sprintf("gene_map_%s.csv", RD_NETWORK)),
                stringsAsFactors = FALSE)
stopifnot(identical(genes, map$symbol))
UQ_DIR <- file.path(RD_TAB_DIR, sprintf("uq_%s", RD_TAG))
dir.create(UQ_DIR, recursive = TRUE, showWarnings = FALSE)

# The size of the data and the schedule of each chain.
n <- nrow(X)
p <- ncol(X)
sched <- lewis_schedule(p, n_iter = RD_N_ITER, burn_in = RD_BURN_IN,
                        max_draws = RD_MAX_DRAWS)
cat(sprintf(paste0("%s: n = %d patients, p = %d genes. %d iterations, burn-in %d, ",
                   "thinning %d, %d draws.\n"),
            RD_TAG, n, p, sched$n_iter, sched$burn_in, sched$thin, sched$n_draws))

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
# A graph is a 0/1 matrix with the parent in the row. D[i, j] = 1 and
# D[j, i] = 0 is the edge i -> j. D[i, j] = D[j, i] = 1 is the undirected
# edge i - j.

# Run LEWIS (complexity) as in the paper. The data are centered already, so
# center = FALSE. The seed is set once by lewis() before the TL start and
# the chain. The check of the effective sample size inside lewis() is off
# (ess_min = 0), because the effective sample size is reported below.
fit_lewis <- function(seed, start = "TL") {
  lewis(X, d = LEWIS_D, alpha = LEWIS_ALPHA, start = start, n_iter = RD_N_ITER,
        burn_in = RD_BURN_IN, max_draws = RD_MAX_DRAWS,
        alpha_sw = LEWIS_ALPHA_SW, alpha_ball = RD_ALPHA_BALL, center = FALSE,
        ess_min = 0, seed = seed)
}

# The edges of a graph as text with node numbers, for example "1->2" and
# "3-4".
edges_of <- function(D) {
  out <- character(0)
  q <- nrow(D)
  for (i in seq_len(q)) for (j in seq_len(q)) if (i != j) {
    if (D[i, j] == 1L && D[j, i] == 0L) out <- c(out, sprintf("%d->%d", i, j))
    else if (D[i, j] == 1L && D[j, i] == 1L && i < j) out <- c(out, sprintf("%d-%d", i, j))
  }
  out
}

# Split the difference between two edge sets into added, removed and
# reversed edges. A reversal is stored "a<->b": the graph has a->b and D_hat
# has b->a.
compare_edges <- function(E_hat, E_cand) {
  added <- setdiff(E_cand, E_hat)
  removed <- setdiff(E_hat, E_cand)
  reversed <- character(0)
  for (e in added) if (grepl("->", e, fixed = TRUE)) {
    ij <- strsplit(e, "->", fixed = TRUE)[[1]]
    rev_e <- sprintf("%s->%s", ij[2L], ij[1L])
    if (rev_e %in% removed) {
      reversed <- c(reversed, sprintf("%s<->%s", ij[1L], ij[2L]))
      added <- setdiff(added, e)
      removed <- setdiff(removed, rev_e)
    }
  }
  list(added = added, removed = removed, reversed = reversed)
}

# The class of a graph relative to D_hat (Step 4 of the post-processing).
classify <- function(D_cand, D_hat) {
  ed <- compare_edges(edges_of(D_hat), edges_of(D_cand))
  na <- length(ed$added)
  nr <- length(ed$removed)
  nv <- length(ed$reversed)
  category <- if (na + nr + nv == 0L) "identical"
              else if (nr == 0L && nv == 0L) "upper"
              else if (na == 0L && nv == 0L) "lower"
              else "rearranged"
  data.frame(category = category, n_edges_added = na, n_edges_removed = nr,
             n_edges_reversed = nv,
             added_edges_str = paste(ed$added, collapse = "|"),
             removed_edges_str = paste(ed$removed, collapse = "|"),
             reversed_edges_str = paste(ed$reversed, collapse = "|"),
             stringsAsFactors = FALSE)
}

# A key of a graph: its entries as a string of 0 and 1 (column by column).
key01 <- function(D) paste(as.integer(D), collapse = "")

# Replace node numbers by gene names, for example "5->1" by "KLK2->AR".
relabel <- function(s, labels) {
  if (!nzchar(s)) return(s)
  m <- gregexpr("[0-9]+", s)
  regmatches(s, m) <- lapply(regmatches(s, m), function(v) labels[as.integer(v)])
  s
}

# The distinct RPDAGs of a fit as matrices, with their number of draws.
# lewis() lists each distinct RPDAG in fit$graphs by the edges that differ
# from D_hat: `added` and `removed` ("a->b" for a directed edge and "a - b"
# for an undirected edge) and `reversed` ("a->b => b->a", the orientation in
# D_hat first). The matrix is rebuilt from D_hat and these lists, and checked
# against the SHD and the number of edges that lewis() reports.
graphs_of <- function(fit) {
  D_hat <- unname(fit$D_hat)
  M <- fit$diagnostics$n_draws
  idx <- function(g) {
    k <- match(g, genes)
    if (is.na(k)) stop("Unknown gene in an edge list: ", g)
    k
  }
  # The two nodes of an edge and whether it is directed.
  ends <- function(e) {
    if (grepl("->", e, fixed = TRUE)) {
      ab <- strsplit(e, "->", fixed = TRUE)[[1]]
      list(a = idx(ab[1L]), b = idx(ab[2L]), directed = TRUE)
    } else {
      ab <- strsplit(e, " - ", fixed = TRUE)[[1]]
      list(a = idx(ab[1L]), b = idx(ab[2L]), directed = FALSE)
    }
  }
  split_list <- function(s) {
    if (is.na(s) || !nzchar(s)) character(0) else strsplit(s, " | ", fixed = TRUE)[[1]]
  }
  gr <- fit$graphs
  mats <- lapply(seq_len(nrow(gr)), function(r) {
    D <- D_hat
    # Remove the edges of D_hat that are not in the graph.
    for (e in split_list(gr$removed[r])) {
      z <- ends(e)
      D[z$a, z$b] <- 0L
      if (!z$directed) D[z$b, z$a] <- 0L
    }
    # Add the edges of the graph that are not in D_hat.
    for (e in split_list(gr$added[r])) {
      z <- ends(e)
      D[z$a, z$b] <- 1L
      if (!z$directed) D[z$b, z$a] <- 1L
    }
    # Reverse the reversed edges.
    for (e in split_list(gr$reversed[r])) {
      parts <- strsplit(e, " => ", fixed = TRUE)[[1]]
      old <- ends(parts[1L])
      new <- ends(parts[2L])
      D[old$a, old$b] <- 0L
      D[new$a, new$b] <- 1L
    }
    D
  })
  # Check every rebuilt graph against the SHD and the edge count of lewis().
  shd_chk <- vapply(mats, function(D) as.numeric(LEWIS::shd(D_hat, D)), numeric(1))
  edg_chk <- vapply(mats, function(D) length(edges_of(D)), numeric(1))
  if (!isTRUE(all(shd_chk == gr$shd)) || !isTRUE(all(edg_chk == gr$n_edges))) {
    stop("A graph rebuilt from fit$graphs does not match its SHD or its number of edges.")
  }
  # The number of draws of each graph.
  cnt <- as.integer(round(gr$mass * M))
  stopifnot(sum(cnt) == M)
  # The non-Gaussian set attached to each graph: the most frequent set among
  # its draws, as node numbers.
  ng <- lapply(gr$nongaussian, function(s) {
    if (is.na(s) || !nzchar(s)) integer(0) else idx_all(strsplit(s, ", ", fixed = TRUE)[[1]])
  })
  list(mats = mats, count = cnt, M = M, nongaussian = ng,
       keys = vapply(mats, key01, character(1)))
}
idx_all <- function(g) {
  k <- match(g, genes)
  if (anyNA(k)) stop("Unknown gene: ", paste(g[is.na(k)], collapse = ", "))
  k
}

# Posterior probabilities over the draws: the proportion of draws with a
# directed edge (dir) and with an adjacency (adj) at every pair.
edge_probs <- function(G) {
  dir <- matrix(0, p, p)
  adj <- matrix(0, p, p)
  for (u in seq_along(G$mats)) {
    D <- G$mats[[u]]
    dir <- dir + G$count[u] * ((D == 1L) & (t(D) == 0L))
    adj <- adj + G$count[u] * ((D == 1L) | (t(D) == 1L))
  }
  list(dir = dir / G$M, adj = adj / G$M)
}

# Effective sample size (Geyer's initial positive sequence estimator), the
# estimator lewis() uses for the log-score.
ess_geyer <- function(x) {
  N <- length(x)
  if (N < 3L || stats::var(x) == 0) return(N)
  r <- stats::acf(x, lag.max = min(N - 1L, 2000L), plot = FALSE)$acf[, 1, 1]
  tau <- -1
  m <- 0L
  while (2L * m + 2L <= length(r)) {
    G <- r[2L * m + 1L] + r[2L * m + 2L]
    if (G <= 0) break
    tau <- tau + 2 * G
    m <- m + 1L
  }
  N / max(tau, 1)
}

# Potential scale reduction factor (Gelman and Rubin) of a list of chains.
psrf <- function(chains) {
  N <- min(lengths(chains))
  x <- vapply(chains, function(v) as.numeric(v[seq_len(N)]), numeric(N))
  W <- mean(apply(x, 2, stats::var))
  B <- N * stats::var(colMeans(x))
  if (W == 0) return(if (B == 0) 1 else Inf)
  sqrt(((N - 1) / N * W + B / N) / W)
}

# -----------------------------------------------------------------------------
# 1. The reported chain
# -----------------------------------------------------------------------------
fit <- fit_lewis(RD_SEED)
D_hat <- unname(fit$D_hat)
nG_hat <- idx_all(fit$nG_hat)
G <- graphs_of(fit)
M <- G$M
cat(sprintf("Chain: acceptance %.3f, %d distinct DAGs, %d draws, %.0f s\n",
            fit$diagnostics$acceptance_rate, fit$diagnostics$n_unique_dags, M,
            fit$diagnostics$runtime))
cat(sprintf("D_hat: %d edges (%d directed, %d undirected). Non-Gaussian genes: %s. Radius %d.\n",
            length(edges_of(D_hat)), sum(D_hat == 1L & t(D_hat) == 0L),
            sum(D_hat == 1L & t(D_hat) == 1L) / 2L,
            paste(genes[nG_hat], collapse = ", "), fit$ball$radius))

# Effective sample size of the log-score and of the number of edges of the
# retained DAGs.
ess_score <- fit$diagnostics$ess_log_score
ess_edges <- ess_geyer(fit$trace$n_edges)
cat(sprintf("Effective sample size of the %d draws: %.0f (score), %.0f (number of edges)\n",
            M, ess_score, ess_edges))

# ---- Every pair of genes ------------------------------------------------------
# The proportion of draws in which the pair is i -> j, j -> i, undirected or
# absent. The posterior probability of an edge of D_hat is the proportion
# of draws in which the pair has its status in D_hat.
status_of <- function(D, i, j) {
  if (D[i, j] == 1L && D[j, i] == 1L) "undirected"
  else if (D[i, j] == 1L) "i->j"
  else if (D[j, i] == 1L) "j->i"
  else "absent"
}
levels_st <- c("i->j", "j->i", "undirected", "absent")
pair_rows <- list()
for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
  st <- vapply(G$mats, status_of, character(1), i = i, j = j)
  tb <- vapply(levels_st, function(s) sum(G$count[st == s]), numeric(1)) / M
  hat_st <- status_of(D_hat, i, j)
  pair_rows[[length(pair_rows) + 1L]] <- data.frame(
    i = i, j = j, gene_i = genes[i], gene_j = genes[j],
    p_i_to_j = tb[["i->j"]], p_j_to_i = tb[["j->i"]],
    p_undirected = tb[["undirected"]], p_absent = tb[["absent"]],
    p_adjacent = 1 - tb[["absent"]],
    status_in_D_hat = hat_st, p_status_in_D_hat = tb[[hat_st]],
    stringsAsFactors = FALSE)
}
pairs <- do.call(rbind, pair_rows)
write.csv(pairs, file.path(RD_TAB_DIR, sprintf("pair_posterior_%s.csv", RD_TAG)),
          row.names = FALSE)

# The edges of D_hat.
pe <- pairs[pairs$status_in_D_hat != "absent", ]
pe$edge <- ifelse(pe$status_in_D_hat == "i->j", paste(pe$gene_i, "->", pe$gene_j),
           ifelse(pe$status_in_D_hat == "j->i", paste(pe$gene_j, "->", pe$gene_i),
                  paste(pe$gene_i, "-", pe$gene_j)))

# ---- Agreement with the known regulation ----------------------------------------
# An edge from a regulator to a target agrees with the known regulation. An
# edge from a target to a regulator goes against it. An edge between genes
# of the same role ("within tier") or an undirected edge is not tested.
tier <- match(map$role, unique(map$role))
pe$known_order <- vapply(seq_len(nrow(pe)), function(r) {
  st <- pe$status_in_D_hat[r]
  if (st == "undirected") return("undirected")
  from <- if (st == "i->j") pe$i[r] else pe$j[r]
  to <- if (st == "i->j") pe$j[r] else pe$i[r]
  if (tier[from] < tier[to]) "agrees" else if (tier[from] > tier[to]) "against" else "within tier"
}, character(1))

# -----------------------------------------------------------------------------
# 2. The genes: classification and Shapiro-Wilk p-values (Table 2)
# -----------------------------------------------------------------------------
# The p-values are those of the residuals of each gene in the first retained
# DAG whose draw gives (D_hat, nG_hat). When D_hat has no undirected edge,
# this DAG is D_hat itself, because every directed edge of an RPDAG has the
# orientation of the DAG it comes from. Then nongaussian_nodes() on D_hat
# must return nG_hat, which is checked. When D_hat has undirected edges, a
# DAG consistent with D_hat is used. Its p-values can then differ from those
# of the first retained draw at the nodes of the undirected edges.
und <- D_hat == 1L & t(D_hat) == 1L
if (!any(und)) {
  dag0 <- D_hat
} else {
  message("D_hat has undirected edges. A DAG consistent with D_hat is used ",
          "for the Shapiro-Wilk p-values.")
  if (!requireNamespace("pcalg", quietly = TRUE)) {
    stop("The package pcalg is needed when D_hat has undirected edges.")
  }
  g0 <- methods::as(D_hat, "graphNEL")
  ext <- pcalg::pdag2dag(g0)
  if (!isTRUE(ext$success)) stop("D_hat has no consistent DAG extension.")
  dag0 <- matrix(as.integer(methods::as(ext$graph, "matrix") != 0), p, p)
}
dimnames(dag0) <- list(genes, genes)
ng0 <- nongaussian_nodes(X, dag0, alpha_sw = LEWIS_ALPHA_SW, center = FALSE)
if (!any(und) && !setequal(ng0, fit$nG_hat)) {
  stop("The non-Gaussian genes of D_hat differ from nG_hat.")
}
node_post <- data.frame(node = seq_len(p), gene = genes, role = map$role,
                        p_non_gaussian = fit$nodes$prob_nongaussian,
                        in_nG_hat = fit$nodes$in_nG_hat,
                        sw_p_value = unname(attr(ng0, "p_values")[genes]),
                        sw_level = attr(ng0, "level"),
                        stringsAsFactors = FALSE)
node_post$class <- ifelse(node_post$in_nG_hat, "non-Gaussian", "Gaussian")
write.csv(node_post, file.path(RD_TAB_DIR, sprintf("node_posterior_%s.csv", RD_TAG)),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# The distinct RPDAGs of the draws (Step 4)
# -----------------------------------------------------------------------------
# The graphs are ranked by their number of draws. Graphs with the same
# number of draws are ranked by their key.
keys_rep <- rep(G$keys, G$count)
tab <- sort(table(keys_rep), decreasing = TRUE)
u_of <- match(names(tab), G$keys)
post <- do.call(rbind, lapply(seq_along(tab), function(r) {
  u <- u_of[r]
  D <- G$mats[[u]]
  nG_r <- G$nongaussian[[u]]
  data.frame(rank = r, D_key = names(tab)[r], mass = as.numeric(tab[r]) / M,
             is_D_hat = identical(names(tab)[r], key01(D_hat)),
             shd_to_Dhat = as.numeric(LEWIS::shd(D_hat, D)), classify(D, D_hat),
             nG_added_str = paste(setdiff(nG_r, nG_hat), collapse = ","),
             nG_removed_str = paste(setdiff(nG_hat, nG_r), collapse = ","),
             n_edges = length(edges_of(D)), stringsAsFactors = FALSE)
}))
post$cum_mass <- cumsum(post$mass)
post$inside_ball <- post$shd_to_Dhat <= fit$ball$radius
post$edits <- vapply(seq_len(nrow(post)), function(r) {
  b <- character(0)
  if (post$n_edges_added[r]) b <- c(b, paste0("+", relabel(post$added_edges_str[r], genes)))
  if (post$n_edges_removed[r]) b <- c(b, paste0("-", relabel(post$removed_edges_str[r], genes)))
  if (post$n_edges_reversed[r]) b <- c(b, paste0("rev ", relabel(post$reversed_edges_str[r], genes)))
  if (!length(b)) "(identical)" else paste(b, collapse = "; ")
}, character(1))
stopifnot(any(post$is_D_hat))
write.csv(post, file.path(RD_TAB_DIR, sprintf("posterior_rpdag_%s.csv", RD_TAG)),
          row.names = FALSE)

# The tables of the posterior figure: the RD_TOP_K most probable graphs, a
# summary of the chain and the gene names.
write.csv(head(post, RD_TOP_K), file.path(UQ_DIR, "posterior_uncertainty_table.csv"),
          row.names = FALSE)
write.csv(data.frame(tag = RD_TAG, n = n, p = p, seed = RD_SEED, alpha_ball = RD_ALPHA_BALL,
                     eps_alpha = fit$ball$radius, M = M,
                     n_distinct_dags = fit$diagnostics$n_unique_dags,
                     n_distinct_rpdags = nrow(post),
                     acceptance_rate = fit$diagnostics$acceptance_rate,
                     mass_D_hat = post$mass[post$is_D_hat][1],
                     mass_outside_ball = fit$ball$mass_outside,
                     n_edges_D_hat = length(edges_of(D_hat)),
                     nG_hat = paste(genes[nG_hat], collapse = ";"),
                     n_iter = sched$n_iter, burn_in = sched$burn_in, thin = sched$thin,
                     runtime_s = fit$diagnostics$runtime,
                     ess_log_score = ess_score, ess_edges = ess_edges),
          file.path(UQ_DIR, "decision_rule_summary.csv"), row.names = FALSE)
write.csv(data.frame(node = seq_len(p), symbol = genes),
          file.path(UQ_DIR, "node_labels.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Stability across seeds
# -----------------------------------------------------------------------------
# The whole procedure, start at the TL estimate included, with ten other
# seeds. For each edge of D_hat, the recurrence is the proportion of these
# runs whose point estimate has the pair in the same status.
stab <- lapply(RD_STABILITY_SEEDS, function(s) {
  f <- fit_lewis(s)
  D <- unname(f$D_hat)
  list(seed = s, D = D, shd = as.numeric(LEWIS::shd(D_hat, D)), n_edges = length(edges_of(D)),
       nG = paste(f$nG_hat, collapse = ";"))
})
same_status <- function(D, i, j) (D[i, j] == D_hat[i, j]) && (D[j, i] == D_hat[j, i])
pe$recurrence <- vapply(seq_len(nrow(pe)), function(r) {
  mean(vapply(stab, function(z) same_status(z$D, pe$i[r], pe$j[r]), logical(1)))
}, numeric(1))
write.csv(pe[, c("edge", "i", "j", "gene_i", "gene_j", "status_in_D_hat", "known_order",
                 "p_status_in_D_hat", "recurrence", "p_adjacent", "p_i_to_j", "p_j_to_i",
                 "p_undirected")],
          file.path(RD_TAB_DIR, sprintf("point_estimate_%s.csv", RD_TAG)), row.names = FALSE)
write.csv(data.frame(seed = vapply(stab, `[[`, numeric(1), "seed"),
                     shd_to_reported = vapply(stab, `[[`, numeric(1), "shd"),
                     n_edges = vapply(stab, `[[`, numeric(1), "n_edges"),
                     nG_hat = vapply(stab, `[[`, character(1), "nG")),
          file.path(RD_TAB_DIR, sprintf("stability_%s.csv", RD_TAG)), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Convergence
# -----------------------------------------------------------------------------
# Four additional chains of the same length: the first half from the TL
# estimate and the second half from the empty graph. The reported chain is
# not one of them. The checks are
#   (i)   the PSRF of the score and of the number of edges of the draws.
#   (ii)  the highest score of the retained draws of each chain.
#   (iii) the largest difference between the chains in the posterior
#         probability of an adjacency and of an orientation.
#   (iv)  the posterior mass of the reported D_hat in each chain.
n_conv <- length(RD_CONV_SEEDS)
from_tl <- seq_len(n_conv) <= ceiling(n_conv / 2)
# Run each chain and keep its trace, its edge probabilities and the mass of
# the reported D_hat.
conv <- lapply(seq_len(n_conv), function(k) {
  f <- fit_lewis(RD_CONV_SEEDS[k], start = if (from_tl[k]) "TL" else "empty")
  Gk <- graphs_of(f)
  pr <- edge_probs(Gk)
  list(seed = RD_CONV_SEEDS[k], start = if (from_tl[k]) "TL" else "empty",
       acceptance = f$diagnostics$acceptance_rate,
       score = f$trace$log_score, n_edges_dag = f$trace$n_edges,
       mean_edges_rpdag = sum(Gk$count * vapply(Gk$mats, function(D) length(edges_of(D)),
                                                numeric(1))) / Gk$M,
       adj = pr$adj, dir = pr$dir,
       mass_D_hat = sum(Gk$count[Gk$keys == key01(D_hat)]) / Gk$M)
})
# PSRF of the score and of the number of edges. The RPDAG of a DAG has the
# same skeleton as the DAG, so a draw and its RPDAG have the same number of
# edges. The PSRF of the number of edges is therefore the same for the DAG
# draws and for the projected draws, and it is written in both columns.
rhat_score <- psrf(lapply(conv, `[[`, "score"))
rhat_dag <- psrf(lapply(conv, `[[`, "n_edges_dag"))
# The largest difference between the chains over a set of entries.
spread <- function(field, sel) {
  v <- vapply(conv, function(z) z[[field]][sel], numeric(sum(sel)))
  max(apply(v, 1, function(r) diff(range(r))))
}
max_d_adj <- spread("adj", upper.tri(D_hat))
max_d_dir <- spread("dir", row(D_hat) != col(D_hat))
conv_tab <- data.frame(chain = seq_len(n_conv), seed = RD_CONV_SEEDS,
                       start = vapply(conv, `[[`, character(1), "start"),
                       acceptance = vapply(conv, `[[`, numeric(1), "acceptance"),
                       mean_edges_dag = vapply(conv, function(z) mean(z$n_edges_dag), numeric(1)),
                       mean_edges_rpdag = vapply(conv, `[[`, numeric(1), "mean_edges_rpdag"),
                       mass_D_hat = vapply(conv, `[[`, numeric(1), "mass_D_hat"),
                       mean_log_score = vapply(conv, function(z) mean(z$score), numeric(1)),
                       max_log_score = vapply(conv, function(z) max(z$score), numeric(1)),
                       psrf_log_score = rhat_score, psrf_edges_dag = rhat_dag,
                       psrf_edges_rpdag = rhat_dag,
                       max_diff_p_adjacent = max_d_adj, max_diff_p_directed = max_d_dir)
write.csv(conv_tab, file.path(RD_TAB_DIR, sprintf("convergence_%s.csv", RD_TAG)),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
ko <- table(factor(pe$known_order, levels = c("agrees", "against", "within tier", "undirected")))
write.csv(data.frame(tag = RD_TAG, n = n, p = p, n_edges = length(edges_of(D_hat)),
                     n_directed = sum(D_hat == 1L & t(D_hat) == 0L),
                     n_non_gaussian = length(nG_hat),
                     non_gaussian = paste(genes[nG_hat], collapse = ";"),
                     n_gaussian = p - length(nG_hat),
                     gaussian = paste(setdiff(genes, genes[nG_hat]), collapse = ";"),
                     agrees = as.integer(ko["agrees"]), against = as.integer(ko["against"]),
                     within_tier = as.integer(ko["within tier"]),
                     undirected = as.integer(ko["undirected"]),
                     mass_D_hat = post$mass[post$is_D_hat][1],
                     eps_alpha = fit$ball$radius,
                     acceptance = fit$diagnostics$acceptance_rate, n_retained = M,
                     ess_log_score = ess_score, ess_edges = ess_edges,
                     n_iter = sched$n_iter, burn_in = sched$burn_in, thin = sched$thin,
                     agreeing_edges = paste(pe$edge[pe$known_order == "agrees"], collapse = "; "),
                     against_edges = paste(pe$edge[pe$known_order == "against"], collapse = "; ")),
          file.path(RD_TAB_DIR, sprintf("summary_%s.csv", RD_TAG)), row.names = FALSE)

# Print the main results.
cat("\n==== ", RD_TAG, " ====\n", sep = "")
cat(sprintf("n = %d patients, p = %d genes, %d draws, %d distinct RPDAGs\n",
            n, p, M, nrow(post)))
cat(sprintf("Point estimate: %d edges, posterior mass %.1f%%. Radius %d. %.1f%% of the mass outside the ball.\n",
            length(edges_of(D_hat)), 100 * post$mass[post$is_D_hat][1], fit$ball$radius,
            100 * fit$ball$mass_outside))
cat(sprintf("Gaussian genes: %s\n", paste(setdiff(genes, genes[nG_hat]), collapse = ", ")))
cat(sprintf("\nShapiro-Wilk p-values (level %.5f):\n", LEWIS_ALPHA_SW / p))
print(node_post[, c("gene", "class", "sw_p_value", "p_non_gaussian")], row.names = FALSE, digits = 3)
cat("\nEdges of the point estimate:\n")
print(pe[order(-pe$p_status_in_D_hat), c("edge", "known_order", "p_status_in_D_hat", "recurrence")],
      row.names = FALSE, digits = 3)
cat("\nStability: SHD of each other seed to D_hat:",
    paste(vapply(stab, `[[`, numeric(1), "shd"), collapse = " "), "\n")
cat(sprintf("Convergence: PSRF %.4f (score), %.4f (edges). Largest difference %.3f (adjacency), %.3f (orientation).\n",
            rhat_score, rhat_dag, max_d_adj, max_d_dir))
cat("Highest score of each chain minus the overall highest:",
    paste(sprintf("%.2f", conv_tab$max_log_score - max(conv_tab$max_log_score)), collapse = " "), "\n")
