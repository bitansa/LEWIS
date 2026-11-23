# =============================================================================
# check_orientation_conventions.R
#
# Checks that every method wrapper returns its graph with the convention
# A[k, j] = 1 for the edge k -> j (the row is the parent).
#
#     Rscript reproducibility/common/check_orientation_conventions.R
#
# WHY. The packages of the methods use different conventions for their
# matrices. A wrapper that reads a matrix the wrong way still returns a
# plausible graph, but with every edge reversed. This cannot be seen in the
# tables. The check uses a graph whose edges are all identifiable.
#
# TEST GRAPH (p = 5): 1 -> 3 <- 2, 3 -> 4, 4 -> 5.
#   1 -> 3 <- 2 is a v-structure. Meek's rule 1 then orients 3 -> 4 and
#   4 -> 5. So the CPDAG of this DAG is the DAG itself, and every method
#   should return the four edges with the right orientation at n = 5000.
#   The errors are Student's t3, so that the LiNGAM methods also apply.
#
# Verdict for each method:
#   OK          the four edges are present with the right orientation
#   TRANSPOSED  the edges are reversed: the wrapper reads the matrix of the
#               package the wrong way. Fix it before using any result.
#   WEAK        edges are missing or undirected. This is not a problem of
#               orientation.
#
# The script stops with an error if a method is TRANSPOSED.
# =============================================================================

# ---- Folder of this script and shared code ----------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- .here
# Use the Python on the PATH when RETICULATE_PYTHON is not set.
if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
  pybin <- Sys.which("python3")
  if (nzchar(pybin)) Sys.setenv(RETICULATE_PYTHON = pybin)
}
suppressPackageStartupMessages({
  source(file.path(COMMON_DIR, "scenario_instances.R"))
  source(file.path(COMMON_DIR, "method_wrappers.R"))
})

# ---- The test graph and data ------------------------------------------------
p <- 5L
n <- 5000L
TRUE_DAG <- matrix(0L, p, p)
TRUE_DAG[1, 3] <- 1L
TRUE_DAG[2, 3] <- 1L
TRUE_DAG[3, 4] <- 1L
TRUE_DAG[4, 5] <- 1L

# Data from the test graph with t3 errors scaled to unit variance.
set.seed(20260901)
Z <- matrix(rt(n * p, df = 3) / sqrt(3), n, p)
X <- matrix(0, n, p)
X[, 1] <- Z[, 1]
X[, 2] <- Z[, 2]
X[, 3] <- 0.8 * X[, 1] - 0.7 * X[, 2] + Z[, 3]
X[, 4] <- 0.9 * X[, 3] + Z[, 4]
X[, 5] <- 0.6 * X[, 4] + Z[, 5]

# ---- Verdict for one graph --------------------------------------------------
verdict <- function(A) {
  A <- (A != 0) * 1L
  correct <- sum(A == 1L & TRUE_DAG == 1L & t(A) == 0L)
  reversed <- sum(A == 1L & t(TRUE_DAG) == 1L & t(A) == 0L)
  if (correct == 4L && reversed == 0L) "OK"
  else if (reversed >= 3L && correct == 0L) "TRANSPOSED"
  else "WEAK"
}

# Compare two graphs by their 0/1 pattern only.
same_graph <- function(A, B) {
  identical(matrix(as.integer(A != 0), nrow(A)), matrix(as.integer(B != 0), nrow(B)))
}

# ---- The reference graphs ---------------------------------------------------
cat("\nTrue DAG   :", edge_list(TRUE_DAG), "\n")
cat("True CPDAG :", edge_list(to_cpdag(TRUE_DAG)), "\n")
stopifnot(same_graph(to_cpdag(TRUE_DAG), TRUE_DAG))
stopifnot(same_graph(rpdag(TRUE_DAG), TRUE_DAG))
stopifnot(same_graph(rpdag(TRUE_DAG, seq_len(p)), TRUE_DAG))
cat("to_cpdag() and LEWIS::rpdag() agree with the true graph.\n")

# The node order must be kept at p = 100. With every node non-Gaussian the
# RPDAG of a DAG is the DAG itself, and the CPDAG of pcalg must agree with
# the RPDAG of LEWIS with no non-Gaussian node.
local({
  inst <- build_instance("ER", 100L)
  stopifnot(same_graph(rpdag(inst$dag, seq_len(100L)), inst$dag))
  stopifnot(same_graph(to_cpdag(inst$dag), rpdag(inst$dag)))
  cat("Node order kept at p = 100. to_cpdag() equals rpdag() with no",
      "non-Gaussian node.\n\n")
})

# ---- Every method -----------------------------------------------------------
cat(sprintf("%-18s %-40s %s\n", "Method", "Edges", "Verdict"))
results <- character(0)
for (m in c("pc", "ges", "lingam", "direct_lingam", "rcd", "tl_zhao", "notears")) {
  require_method(m)
  out <- run_method(m, X)
  if (isTRUE(out$skipped)) stop("Method ", m, " failed: ", out$error)
  A <- est_partial_dag(m, out$dag)
  results[m] <- verdict(A)
  cat(sprintf("%-18s %-40s %s\n", m, edge_list(A), results[m]))
}

# LEWIS on a short chain. Only the orientation is checked here.
fit <- lewis(X, n_iter = 2000, center = FALSE, ess_min = 0, seed = 1)
results["lewis_complexity"] <- verdict(fit$D_hat)
cat(sprintf("%-18s %-40s %s\n", "lewis_complexity", edge_list(fit$D_hat),
            results["lewis_complexity"]))

# The start graph of LEWIS reads the matrix of tl_lingam() and must agree with
# the TL wrapper.
tl <- run_method("tl_zhao", X)$dag
if (same_graph(fit$diagnostics$start_dag, tl)) {
  cat("The start graph of lewis() agrees with the TL wrapper.\n\n")
} else {
  cat("WARNING: the start graph of lewis() differs from the TL wrapper.\n",
      "Check that both read the matrix B of tl_lingam() with the parent in the row.\n\n")
}

# Stop if a wrapper reads its matrix the wrong way.
if (any(results == "TRANSPOSED")) {
  stop("Some wrappers return reversed edges: ",
       paste(names(results)[results == "TRANSPOSED"], collapse = ", "))
}
cat("No wrapper returns reversed edges.\n")
