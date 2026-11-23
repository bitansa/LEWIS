# Input validation. Every check stops with a short message that says what
# is wrong and what to do.

# Check the data and return a list with the numeric matrix, its size and
# the variable names. With center = TRUE the column means are subtracted.
.check_data <- function(X, center = TRUE) {
  # Convert a data frame with numeric columns to a matrix.
  if (is.data.frame(X)) {
    is_num <- vapply(X, is.numeric, logical(1))
    if (!all(is_num)) {
      stop("All columns of X must be numeric. Remove or recode the other columns.",
           call. = FALSE)
    }
    X <- as.matrix(X)
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix or a data frame with numeric columns.",
         call. = FALSE)
  }
  n <- nrow(X)
  p <- ncol(X)

  # Check the dimensions.
  if (p < 2L) {
    stop("X must have at least 2 columns (variables).", call. = FALSE)
  }
  if (n < 3L) {
    stop("X must have at least 3 rows (observations).", call. = FALSE)
  }
  if (n <= p) {
    stop(sprintf("X has n = %d rows and p = %d columns. LEWIS needs n > p. ", n, p),
         "Use more observations or fewer variables.", call. = FALSE)
  }

  # Check for missing and infinite values.
  if (anyNA(X) || any(!is.finite(X))) {
    stop("X contains missing or infinite values. Remove or impute them first.",
         call. = FALSE)
  }

  # Variable names. Missing names become V1, ..., Vp.
  nm <- colnames(X)
  if (is.null(nm)) nm <- paste0("V", seq_len(p))
  empty <- is.na(nm) | nm == ""
  nm[empty] <- paste0("V", seq_len(p))[empty]
  if (anyDuplicated(nm)) {
    stop("The column names of X must be unique.", call. = FALSE)
  }
  storage.mode(X) <- "double"
  dimnames(X) <- list(NULL, nm)

  # Constant columns carry no information about the graph.
  is_const <- apply(X, 2L, function(v) all(v == v[1L]))
  if (any(is_const)) {
    stop("These columns of X are constant: ", paste(nm[is_const], collapse = ", "),
         ". Remove them.", call. = FALSE)
  }

  # Duplicated columns make the regressions degenerate.
  dup <- duplicated(t(X))
  if (any(dup)) {
    stop("These columns of X duplicate earlier columns: ",
         paste(nm[dup], collapse = ", "), ". Remove them.", call. = FALSE)
  }

  # Center the columns if requested.
  if (isTRUE(center)) {
    X <- X - matrix(colMeans(X), n, p, byrow = TRUE)
  }

  # Linearly dependent columns make the score infinite.
  if (qr(X, tol = 1e-7)$rank < p) {
    stop("The columns of X are linearly dependent. ",
         "Remove columns that are linear combinations of others.", call. = FALSE)
  }

  list(X = X, n = n, p = p, names = nm)
}

# Check a single logical value.
.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
  }
  invisible(x)
}

# Check a single finite number.
.check_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("%s must be a single finite number.", name), call. = FALSE)
  }
  invisible(x)
}

# Check a single whole number of at least `min`.
.check_count <- function(x, name, min = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x) ||
      x < min || x > .Machine$integer.max) {
    stop(sprintf("%s must be a whole number of at least %d.", name, as.integer(min)),
         call. = FALSE)
  }
  invisible(as.integer(x))
}

# Check a level strictly between 0 and 1.
.check_level <- function(x, name) {
  .check_number(x, name)
  if (x <= 0 || x >= 1) {
    stop(sprintf("%s must lie strictly between 0 and 1.", name), call. = FALSE)
  }
  invisible(x)
}

# Check the constants d and alpha of the complexity prior. alpha is only
# checked when d > 0, because d = 0 gives the uniform prior.
.check_prior <- function(d, alpha) {
  .check_number(d, "d")
  if (d < 0) {
    stop("d must be zero or positive. Use d = 0 for the uniform prior.", call. = FALSE)
  }
  if (d > 0) {
    .check_number(alpha, "alpha")
    if (alpha <= 0.5 || alpha >= 1) {
      stop("alpha must lie strictly between 0.5 and 1. ",
           "The consistency theorem for the complexity prior in the paper requires this.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Check an adjacency matrix of a DAG with p nodes and return it as an
# integer matrix without names.
.check_dag <- function(dag, p, names = NULL, name = "dag") {
  if (is.data.frame(dag)) dag <- as.matrix(dag)
  if (!is.matrix(dag) || !(is.numeric(dag) || is.logical(dag))) {
    stop(sprintf("%s must be a numeric %d x %d adjacency matrix.", name, p, p),
         call. = FALSE)
  }
  if (nrow(dag) != p || ncol(dag) != p) {
    stop(sprintf("%s must be a %d x %d matrix, one row and one column per variable.",
                 name, p, p), call. = FALSE)
  }
  if (anyNA(dag) || !all(dag == 0 | dag == 1)) {
    stop(sprintf("%s must contain only 0 and 1.", name), call. = FALSE)
  }
  if (any(diag(dag) != 0)) {
    stop(sprintf("%s must have a zero diagonal.", name), call. = FALSE)
  }
  # Names, if present, must match the variable names.
  if (!is.null(names) && !is.null(colnames(dag)) && !identical(colnames(dag), names)) {
    stop(sprintf("The column names of %s must match the column names of X.", name),
         call. = FALSE)
  }
  A <- matrix(as.integer(dag), p, p)
  if (!lewis_is_acyclic_cpp(A)) {
    stop(sprintf("%s has a directed cycle. It must be a DAG.", name), call. = FALSE)
  }
  A
}

# Check a partially directed graph (0/1 matrix with a zero diagonal).
.check_pdag <- function(D, name) {
  if (is.data.frame(D)) D <- as.matrix(D)
  if (!is.matrix(D) || !(is.numeric(D) || is.logical(D)) || nrow(D) != ncol(D)) {
    stop(sprintf("%s must be a square adjacency matrix.", name), call. = FALSE)
  }
  if (nrow(D) < 2L) {
    stop(sprintf("%s must have at least 2 nodes.", name), call. = FALSE)
  }
  if (anyNA(D) || !all(D == 0 | D == 1)) {
    stop(sprintf("%s must contain only 0 and 1.", name), call. = FALSE)
  }
  if (any(diag(D) != 0)) {
    stop(sprintf("%s must have a zero diagonal.", name), call. = FALSE)
  }
  matrix(as.integer(D), nrow(D), ncol(D))
}

# Check the start argument of lewis() and return a list with the type and,
# for a user matrix, the matrix.
.check_start <- function(start, p, names) {
  if (is.character(start)) {
    if (length(start) != 1L || !(start %in% c("TL", "empty"))) {
      stop('start must be "TL", "empty" or a p x p adjacency matrix of a DAG.',
           call. = FALSE)
    }
    return(list(type = start, dag = NULL))
  }
  list(type = "user", dag = .check_dag(start, p, names, name = "start"))
}

# Check the seed: NULL or a single whole number.
.check_seed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed) || abs(seed) > .Machine$integer.max) {
    stop("seed must be NULL or a single whole number.", call. = FALSE)
  }
  invisible(as.integer(seed))
}

# Turn a set of non-Gaussian nodes into sorted 1-based indices. The set may
# be given by indices, by names or by a logical vector of length p.
.check_nodes <- function(nG, p, names) {
  if (is.null(nG) || length(nG) == 0L) return(integer(0))
  if (is.logical(nG)) {
    if (length(nG) != p || anyNA(nG)) {
      stop(sprintf("A logical nG must have length %d and no missing values.", p),
           call. = FALSE)
    }
    return(which(nG))
  }
  if (is.character(nG)) {
    idx <- match(nG, names)
    if (anyNA(idx)) {
      stop("Unknown node names in nG: ", paste(nG[is.na(idx)], collapse = ", "), ".",
           call. = FALSE)
    }
    return(sort(unique(idx)))
  }
  if (is.numeric(nG)) {
    if (anyNA(nG) || any(nG != round(nG)) || any(nG < 1) || any(nG > p)) {
      stop(sprintf("Numeric nG must contain node indices between 1 and %d.", p),
           call. = FALSE)
    }
    return(sort(unique(as.integer(nG))))
  }
  stop("nG must be node indices, node names or a logical vector.", call. = FALSE)
}
