# =============================================================================
# method_wrappers.R
#
# One wrapper for each method of the comparative study (Section 7.2) and of
# Figure 1. Each wrapper takes the data matrix X and returns
#
#     list(dag = <p x p adjacency matrix>, runtime = <seconds>,
#          skipped = <TRUE if the method failed>, error = <message or "">)
#
# with A[k, j] = 1 for the edge k -> j and A[i, j] = A[j, i] = 1 for an
# undirected edge.
#
# FAILURES. A wrapper catches an error of the method, prints the message at
# once with message(), and returns the empty graph with skipped = TRUE and
# the message in `error`. The calling script counts these rows and exits
# with a non-zero status if there is any. A missing package is not a
# failure of the method: require_method() stops before any dataset is run.
#
# METHODS AND SOFTWARE (Appendix G.1 and Appendix J of the paper).
#     lewis_complexity  LEWIS::lewis() with d = 0.03, alpha = 0.75   RPDAG
#     lewis_uniform     LEWIS::lewis() with d = 0                    RPDAG
#     pc                pcalg::pc, Gaussian tests at level 0.05      CPDAG
#     ges               pcalg::ges, Gaussian BIC score               CPDAG
#     notears           Python notears.linear.notears_linear         DAG
#     lingam            Python lingam.ICALiNGAM                      DAG
#     direct_lingam     Python lingam.DirectLiNGAM                   DAG
#     rcd               Python lingam.RCD                            DAG
#     tl_zhao           LEWIS::tl_lingam (TL code of TransGraph 1.1.0) DAG
#     pairwise_lingam   MATLAB pwling. It is run by the export, MATLAB
#                       and import scripts of 02_comparison and
#                       01_motivation, not by run_method().
#
# LEWIS is run with center = FALSE, as in the paper. The chain consumes the
# random numbers of the R session from the state left by draw_dataset(), as
# in the runs of the paper, so no seed is passed to lewis(). The check of
# the effective sample size is switched off (ess_min = 0). It does not
# change the results.
#
# EVALUATION (Section 7.2). est_partial_dag() maps the output of each method
# to the graph that is compared with the true RPDAG:
#     LEWIS ............ the RPDAG D_hat, as returned
#     PC, GES .......... the CPDAG, as returned
#     NOTEARS .......... the CPDAG of the returned DAG (least-squares loss
#                        identifies only the Markov equivalence class)
#     ICA-LiNGAM, DirectLiNGAM, RCD, TL, Pairwise-LiNGAM
#                        the returned DAG (these methods assume that all
#                        errors are non-Gaussian, and then the RPDAG of a
#                        DAG is the DAG itself)
#
# The calling script must define COMMON_DIR before it sources this file.
# Dependencies: LEWIS, pcalg, reticulate with the Python
# packages lingam and notears.
# =============================================================================

stopifnot(exists("COMMON_DIR"))
source(file.path(COMMON_DIR, "settings.R"))
source(file.path(COMMON_DIR, "graph_tools.R"))

# pcalg also has a function shd(). LEWIS is put ahead of pcalg on the search
# path, so that shd() is LEWIS::shd(). If LEWIS is attached already, it is
# detached and attached again. Scripts still call LEWIS::shd() explicitly.
suppressPackageStartupMessages({
  library(methods)
  library(pcalg)
  if ("package:LEWIS" %in% search()) detach("package:LEWIS")
  library(LEWIS)
})

# -----------------------------------------------------------------------------
# Python modules, imported once per R session.
# -----------------------------------------------------------------------------
.py_cache <- new.env(parent = emptyenv())

# The Python package lingam (ICA-LiNGAM, DirectLiNGAM, RCD).
.import_lingam <- function() {
  if (is.null(.py_cache$lingam)) {
    .py_cache$lingam <- reticulate::import("lingam", delay_load = FALSE)
  }
  .py_cache$lingam
}

# The function notears_linear of the authors' Python package notears 3.0
# (module notears.linear).
.import_notears <- function() {
  if (is.null(.py_cache$notears)) {
    mod <- reticulate::import("notears.linear", delay_load = FALSE)
    .py_cache$notears <- mod$notears_linear
  }
  .py_cache$notears
}

# -----------------------------------------------------------------------------
# require_method: stop at once if the software of a method is missing.
# -----------------------------------------------------------------------------
require_method <- function(name) {
  # Check that the needed package can be loaded, and stop if it cannot.
  need <- function(ok, what) {
    if (!isTRUE(ok)) stop("Method '", name, "' needs ", what, ".", call. = FALSE)
  }
  # Try to import a Python module. Return "" or the error message.
  import_error <- function(import) {
    tryCatch({
      import()
      ""
    }, error = function(e) conditionMessage(e))
  }
  switch(name,
    lewis_complexity = , lewis_uniform =
      need(requireNamespace("LEWIS", quietly = TRUE), "the R package LEWIS"),
    pc = , ges =
      need(requireNamespace("pcalg", quietly = TRUE), "the R package pcalg"),
    tl_zhao =
      need(requireNamespace("LEWIS", quietly = TRUE), "the R package LEWIS"),
    lingam = , direct_lingam = , rcd = {
      msg <- import_error(.import_lingam)
      need(!nzchar(msg), paste0("the Python package lingam (through reticulate). ",
                                "The import failed with the message: ", msg))
    },
    notears = {
      msg <- import_error(.import_notears)
      need(!nzchar(msg), paste0("the Python package notears 3.0 (module notears.linear). ",
                                "The import failed with the message: ", msg))
    },
    stop("Unknown method: ", name, call. = FALSE))
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Helpers.
# -----------------------------------------------------------------------------

# Run `expr`, time it, and turn an error into a flagged empty graph. The
# message is printed at once, so that it appears in the log of the job.
.run_timed <- function(name, p, expr) {
  t0 <- Sys.time()
  res <- tryCatch(list(dag = expr, error = ""),
                  error = function(e) {
                    msg <- conditionMessage(e)
                    message(sprintf("Method %s failed: %s", name, msg))
                    list(dag = matrix(0L, p, p), error = msg)
                  })
  list(dag = res$dag,
       runtime = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       skipped = nzchar(res$error), error = res$error)
}

# A coefficient matrix of the Python package lingam has B[i, j] != 0 for the
# edge j -> i (the row is the child). Convert it to A[k, j] = 1 for k -> j.
.b_to_adjacency <- function(Bmat, p) {
  A <- matrix(0L, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) {
    if (i != j && abs(Bmat[i, j]) > 1e-8) A[j, i] <- 1L
  }
  A
}

# LEWIS::tl_lingam and notears_linear return a matrix with B[i, j] != 0
# for the edge i -> j (the row is the parent). Convert it to a 0/1 matrix.
.b_to_adjacency_rowparent <- function(Bmat, p) {
  A <- matrix(0L, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) {
    if (i != j && abs(Bmat[i, j]) > 1e-8) A[i, j] <- 1L
  }
  A
}

# -----------------------------------------------------------------------------
# Wrappers.
# -----------------------------------------------------------------------------

# LEWIS. The point estimate D_hat is an RPDAG.
run_lewis <- function(X, prior = c("complexity", "uniform")) {
  prior <- match.arg(prior)
  name <- paste0("lewis_", prior)
  .run_timed(name, ncol(X), {
    fit <- if (prior == "complexity") {
      lewis(X, d = LEWIS_D, alpha = LEWIS_ALPHA, center = FALSE, ess_min = 0)
    } else {
      lewis(X, d = 0, center = FALSE, ess_min = 0)
    }
    unname(fit$D_hat)
  })
}

# PC (Spirtes et al.) with Gaussian partial-correlation tests at level 0.05.
# The CPDAG is read from the graphNEL object, whose matrix form has
# M[a, b] = 1 for the edge a -> b. The coercion as(fit, "amat") is not
# used, because pcalg's amat type has the parent in the column.
run_pc <- function(X) {
  p <- ncol(X)
  n <- nrow(X)
  .run_timed("pc", p, {
    suff <- list(C = cor(X), n = n)
    fit <- pc(suffStat = suff, indepTest = gaussCItest, alpha = 0.05,
              labels = as.character(seq_len(p)), verbose = FALSE)
    M <- as(fit@graph, "matrix")
    M <- matrix(as.integer(M != 0), nrow = p, ncol = p)
    diag(M) <- 0L
    M
  })
}

# GES (Chickering 2002) with the Gaussian BIC score. The CPDAG is read from
# the parent lists of the essential graph. An undirected edge appears in the
# parent lists of both of its nodes.
run_ges <- function(X) {
  p <- ncol(X)
  .run_timed("ges", p, {
    score <- new("GaussL0penObsScore", data = X)
    fit <- ges(score)
    in_edges <- fit$essgraph$.in.edges
    stopifnot(is.list(in_edges), length(in_edges) == p)
    A <- matrix(0L, p, p)
    for (j in seq_len(p)) {
      pa <- as.integer(in_edges[[j]])
      if (length(pa) > 0L) A[pa, j] <- 1L
    }
    diag(A) <- 0L
    A
  })
}

# ICA-LiNGAM (Shimizu et al. 2006) through the Python package lingam. It
# prunes the coefficients by the adaptive lasso. FastICA is random. Its seed
# is drawn from the random numbers of the R session, which the calling
# script has set for each dataset, so the fit is reproducible.
run_lingam <- function(X) {
  p <- ncol(X)
  py <- .import_lingam()
  .run_timed("lingam", p, {
    seed <- sample.int(.Machine$integer.max, 1L)
    model <- py$ICALiNGAM(random_state = as.integer(seed))
    model$fit(X)
    Bmat <- as.matrix(model$adjacency_matrix_)
    Bmat[!is.finite(Bmat)] <- 0
    .b_to_adjacency(Bmat, p)
  })
}

# DirectLiNGAM (Shimizu et al. 2011) through the Python package lingam.
run_direct_lingam <- function(X) {
  p <- ncol(X)
  py <- .import_lingam()
  .run_timed("direct_lingam", p, {
    model <- py$DirectLiNGAM()
    model$fit(X)
    Bmat <- as.matrix(model$adjacency_matrix_)
    Bmat[!is.finite(Bmat)] <- 0
    .b_to_adjacency(Bmat, p)
  })
}

# RCD (Maeda and Shimizu 2020) through the Python package lingam. RCD marks
# the pairs it finds confounded with NaN. These pairs are treated as
# non-adjacent (Appendix G.1).
run_rcd <- function(X) {
  p <- ncol(X)
  py <- .import_lingam()
  .run_timed("rcd", p, {
    model <- py$RCD()
    model$fit(X)
    Bmat <- as.matrix(model$adjacency_matrix_)
    Bmat[!is.finite(Bmat)] <- 0
    .b_to_adjacency(Bmat, p)
  })
}

# TL (Zhao et al. 2022) through LEWIS::tl_lingam, a copy of the TL code of
# the R package TransGraph 1.1.0. Its matrix B has the parent in the row.
run_tl_zhao <- function(X) {
  p <- ncol(X)
  .run_timed("tl_zhao", p, {
    fit <- LEWIS::tl_lingam(as.matrix(X))
    Bmat <- as.matrix(fit$B)
    Bmat[!is.finite(Bmat)] <- 0
    .b_to_adjacency_rowparent(Bmat, p)
  })
}

# NOTEARS (Zheng et al. 2018), linear version, with least-squares loss,
# l1 penalty 0.1 and edge threshold 0.3. The matrix W has the parent in the
# row.
run_notears <- function(X) {
  p <- ncol(X)
  fn <- .import_notears()
  .run_timed("notears", p, {
    W <- as.matrix(fn(X, lambda1 = 0.1, loss_type = "l2", w_threshold = 0.3))
    .b_to_adjacency_rowparent(W, p)
  })
}

# -----------------------------------------------------------------------------
# Dispatch.
# -----------------------------------------------------------------------------
run_method <- function(name, X) {
  # Every study has n > p. Stop if a dataset does not.
  if (nrow(X) <= ncol(X)) {
    stop(sprintf("The dataset has n = %d and p = %d. All studies need n > p.",
                 nrow(X), ncol(X)))
  }
  switch(name,
    lewis_complexity = run_lewis(X, "complexity"),
    lewis_uniform    = run_lewis(X, "uniform"),
    pc               = run_pc(X),
    ges              = run_ges(X),
    lingam           = run_lingam(X),
    direct_lingam    = run_direct_lingam(X),
    rcd              = run_rcd(X),
    notears          = run_notears(X),
    tl_zhao          = run_tl_zhao(X),
    stop("Unknown method: ", name))
}

# -----------------------------------------------------------------------------
# est_partial_dag: the graph of a method that is compared with the truth.
# -----------------------------------------------------------------------------
.LEWIS_METHODS <- c("lewis_complexity", "lewis_uniform")
.CPDAG_METHODS <- c("pc", "ges")
.DAG_TO_CPDAG_METHODS <- c("notears")
.DAG_METHODS <- c("lingam", "direct_lingam", "rcd", "tl_zhao",
                  "pairwise_lingam")

# Map the output of a method to the graph that is compared with the true RPDAG.
est_partial_dag <- function(method_name, est_dag) {
  est_dag <- matrix(as.integer(est_dag != 0), nrow(est_dag), ncol(est_dag))
  if (method_name %in% .LEWIS_METHODS) return(est_dag)
  if (method_name %in% .CPDAG_METHODS) return(est_dag)
  if (method_name %in% .DAG_TO_CPDAG_METHODS) return(to_cpdag(est_dag))
  if (method_name %in% .DAG_METHODS) return(est_dag)
  stop("No evaluation rule for method '", method_name, "'.")
}
