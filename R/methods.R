#' Print a LEWIS fit
#'
#' Prints a short description of a fit of [lewis()]: the data size, the
#' prior, the schedule of the chain and the effective sample size of the
#' score, the edges and non-Gaussian nodes of the point estimate, the
#' radius of the credible ball and any warnings.
#'
#' @param x An object of class `"lewis"`.
#' @param ... Not used.
#'
#' @return `x`, invisibly.
#'
#' @seealso [lewis()], [summary.lewis()]
#'
#' @examples
#' # Simulated data from the chain x1 -> x2 -> x3.
#' set.seed(1)
#' n <- 300
#' x1 <- rexp(n) - rexp(n)
#' x2 <- 0.8 * x1 + rnorm(n)
#' x3 <- 0.8 * x2 + rexp(n) - rexp(n)
#' fit <- lewis(cbind(x1, x2, x3), seed = 1)
#' print(fit)
#'
#' @export
print.lewis <- function(x, ...) {
  s <- x$settings
  dg <- x$diagnostics
  # Data and prior.
  cat("LEWIS fit\n")
  cat(sprintf("  Data: n = %d observations, p = %d variables%s\n", s$n, s$p,
              if (isTRUE(s$center)) " (centered)" else ""))
  if (s$d == 0) {
    cat("  Prior: uniform (d = 0)\n")
  } else {
    cat(sprintf("  Prior: complexity, d = %g, alpha = %g\n", s$d, s$alpha))
  }
  # Chain.
  cat(sprintf("  Chain: %d iterations, burn-in %d, thinning %d, %d draws, start %s\n",
              dg$n_iter, dg$burn_in, dg$thin, dg$n_draws, dg$start))
  cat(sprintf("  Effective sample size of the score: %.0f\n",
              as.numeric(dg$ess_log_score)))
  # Point estimate.
  n_dir <- sum(x$edges$type == "directed")
  n_und <- sum(x$edges$type == "undirected")
  cat(sprintf("  Point estimate: %d edges (%d directed, %d undirected)\n",
              n_dir + n_und, n_dir, n_und))
  cat(sprintf("  Non-Gaussian nodes: %s\n",
              if (length(x$nG_hat) == 0L) "none" else paste(x$nG_hat, collapse = ", ")))
  cat(sprintf("  Credible ball (%g%%): radius %d, %.1f%% of the mass outside\n",
              100 * (1 - x$ball$alpha_ball), x$ball$radius, 100 * x$ball$mass_outside))
  # Warnings.
  if (length(dg$warnings) > 0L) {
    cat("  Warnings:\n")
    for (w in dg$warnings) cat("   -", w, "\n")
  }
  invisible(x)
}

#' Summarize a LEWIS fit
#'
#' Collects the main tables of a fit of [lewis()]: the edges of the point
#' estimate, the nodes, the most probable graphs and the diagnostics.
#'
#' @param object An object of class `"lewis"`.
#' @param top Positive whole number. The number of graphs to show. The
#'   default is 10.
#' @param ... Not used.
#'
#' @return `summary()` returns an object of class `"summary.lewis"`, a list
#'   with elements `edges`, `nodes`, `graphs` (the first `top` rows of
#'   `object$graphs`), `ball`, `diagnostics`, `settings` and `class_mass`
#'   (the posterior mass of each class of graphs). Its `print()` method
#'   prints these tables and returns `x` invisibly.
#'
#' @seealso [lewis()], [print.lewis()]
#'
#' @examples
#' # Simulated data from the chain x1 -> x2 -> x3.
#' set.seed(1)
#' n <- 300
#' x1 <- rexp(n) - rexp(n)
#' x2 <- 0.8 * x1 + rnorm(n)
#' x3 <- 0.8 * x2 + rexp(n) - rexp(n)
#' fit <- lewis(cbind(x1, x2, x3), seed = 1)
#' summary(fit)
#'
#' @export
summary.lewis <- function(object, top = 10, ...) {
  top <- .check_count(top, "top", min = 1)
  # Posterior mass of each class of graphs.
  classes <- c("identical", "upper", "lower", "rearranged")
  class_mass <- vapply(classes, function(k) sum(object$graphs$mass[object$graphs$class == k]),
                       numeric(1))
  structure(list(edges = object$edges, nodes = object$nodes,
                 graphs = utils::head(object$graphs, top), ball = object$ball,
                 diagnostics = object$diagnostics, settings = object$settings,
                 class_mass = class_mass),
            class = "summary.lewis")
}

#' @rdname summary.lewis
#' @param x An object of class `"summary.lewis"`.
#' @param digits Positive whole number. The number of digits for
#'   probabilities. The default is 3.
#' @export
print.summary.lewis <- function(x, digits = 3, ...) {
  digits <- .check_count(digits, "digits", min = 1)
  dg <- x$diagnostics
  s <- x$settings
  # Header.
  cat(sprintf("LEWIS fit: n = %d, p = %d, %s\n\n", s$n, s$p,
              if (s$d == 0) "uniform prior" else
                sprintf("complexity prior (d = %g, alpha = %g)", s$d, s$alpha)))

  # Edges of the point estimate.
  cat("Edges of the point estimate:\n")
  if (nrow(x$edges) == 0L) {
    cat("  (none)\n")
  } else {
    e <- x$edges
    e$prob <- round(e$prob, digits)
    e$prob_adjacent <- round(e$prob_adjacent, digits)
    print(e, row.names = FALSE)
  }

  # Nodes.
  cat("\nNodes:\n")
  nd <- x$nodes
  nd$prob_nongaussian <- round(nd$prob_nongaussian, digits)
  print(nd, row.names = FALSE)

  # Most probable graphs.
  cat(sprintf("\nMost probable graphs (%d shown):\n", nrow(x$graphs)))
  g <- x$graphs
  g$mass <- round(g$mass, digits)
  print(g[, c("rank", "mass", "class", "shd", "n_edges", "added", "removed",
              "reversed", "nongaussian")], row.names = FALSE)
  cat("\nPosterior mass by class:\n")
  print(round(x$class_mass, digits))

  # Credible ball and diagnostics.
  cat(sprintf("\nCredible ball (%g%%): radius %d, %.1f%% of the mass outside\n",
              100 * (1 - x$ball$alpha_ball), x$ball$radius, 100 * x$ball$mass_outside))
  cat("\nDiagnostics:\n")
  cat(sprintf("  iterations %d, burn-in %d, thinning %d, draws %d, start %s\n",
              dg$n_iter, dg$burn_in, dg$thin, dg$n_draws, dg$start))
  cat(sprintf("  acceptance rate %.3f, distinct DAGs %d, distinct RPDAGs %d\n",
              dg$acceptance_rate, dg$n_unique_dags, dg$n_unique_rpdags))
  cat(sprintf("  effective sample size of the score %.0f\n",
              as.numeric(dg$ess_log_score)))
  if (!is.null(dg$psrf_log_score)) {
    cat(sprintf("  PSRF: score %.3f, number of edges %.3f\n",
                dg$psrf_log_score, dg$psrf_n_edges))
    cat(sprintf("  largest difference between chains: adjacency %.3f, orientation %.3f\n",
                dg$max_diff_adjacency, dg$max_diff_orientation))
  }
  cat(sprintf("  runtime %.1f seconds\n", dg$runtime))
  if (length(dg$warnings) > 0L) {
    cat("\nWarnings:\n")
    for (w in dg$warnings) cat("  -", w, "\n")
  }
  invisible(x)
}

#' Plot a LEWIS fit
#'
#' Draws one of four plots of a fit of [lewis()] with base graphics.
#'
#' @details
#' `type = "graph"` draws the point estimate with the nodes on a circle.
#' Directed edges have arrows and undirected edges have none. An edge is
#' solid if its posterior probability is at least 0.9 and dashed
#' otherwise. Nodes declared non-Gaussian are filled black.
#'
#' `type = "ball"` draws the posterior mass of the projected draws by their
#' SHD from the point estimate, split by class. The dashed line marks the
#' radius of the credible ball.
#'
#' `type = "trace"` draws the score \eqn{s_n} of the retained draws against
#' the iteration. It should show no trend after the burn-in.
#'
#' `type = "edges"` draws the posterior probability of each edge of the
#' point estimate as a bar, and its adjacency probability as a point.
#'
#' @param x An object of class `"lewis"`.
#' @param y Not used. The generic [base::plot()] has this argument. Giving
#'   `y` is an error. Choose the plot with `type`.
#' @param type Character. One of `"graph"` (the default), `"ball"`,
#'   `"trace"` or `"edges"`.
#' @param ... Further arguments passed to [graphics::barplot()] for
#'   `type = "ball"` and `type = "edges"`, and to [base::plot()] for
#'   `type = "trace"`. They are not used for `type = "graph"`.
#'
#' @return `x`, invisibly.
#'
#' @seealso [lewis()]
#'
#' @examples
#' # Simulated data from the chain x1 -> x2 -> x3.
#' set.seed(1)
#' n <- 300
#' x1 <- rexp(n) - rexp(n)
#' x2 <- 0.8 * x1 + rnorm(n)
#' x3 <- 0.8 * x2 + rexp(n) - rexp(n)
#' fit <- lewis(cbind(x1, x2, x3), seed = 1)
#' plot(fit)
#' plot(fit, type = "ball")
#' plot(fit, type = "trace")
#' plot(fit, type = "edges")
#'
#' @export
plot.lewis <- function(x, y, type = c("graph", "ball", "trace", "edges"), ...) {
  # The second argument of the generic is not used.
  if (!missing(y)) {
    stop('The argument y is not used. Choose the plot with type, ',
         'for example plot(fit, type = "ball").', call. = FALSE)
  }
  # Draw the chosen plot.
  type <- match.arg(type)
  switch(type,
         graph = .plot_graph(x),
         ball = .plot_ball(x, ...),
         trace = .plot_trace(x, ...),
         edges = .plot_edges(x, ...))
  invisible(x)
}

#' Draw the point estimate on a circle
#' @noRd
.plot_graph <- function(x) {
  D <- x$D_hat
  nm <- colnames(D)
  p <- length(nm)
  # Node positions on a circle, the first node at the top.
  ang <- pi / 2 - 2 * pi * (seq_len(p) - 1) / p
  px <- cos(ang)
  py <- sin(ang)
  r <- min(0.15, 0.9 * sin(pi / p))
  # Empty plot region.
  op <- graphics::par(mar = c(1, 1, 2, 1))
  on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25), asp = 1)
  graphics::title(main = "Point estimate")
  # Edges. Each segment stops at the boundary of the node circles.
  e <- x$edges
  if (nrow(e) > 0L) {
    a <- match(e$from, nm)
    b <- match(e$to, nm)
    dx <- px[b] - px[a]
    dy <- py[b] - py[a]
    len <- sqrt(dx^2 + dy^2)
    x0 <- px[a] + r * dx / len
    y0 <- py[a] + r * dy / len
    x1 <- px[b] - r * dx / len
    y1 <- py[b] - r * dy / len
    lty <- ifelse(e$prob >= 0.9, 1, 2)
    for (k in seq_len(nrow(e))) {
      if (e$type[k] == "directed") {
        graphics::arrows(x0[k], y0[k], x1[k], y1[k], length = 0.1, lty = lty[k], lwd = 1.5)
      } else {
        graphics::segments(x0[k], y0[k], x1[k], y1[k], lty = lty[k], lwd = 1.5)
      }
    }
  }
  # Nodes. Non-Gaussian nodes are black with white labels.
  ng <- nm %in% x$nG_hat
  graphics::symbols(px, py, circles = rep(r, p), inches = FALSE, add = TRUE,
                    bg = ifelse(ng, "black", "white"), fg = "black")
  graphics::text(px, py, labels = nm, col = ifelse(ng, "white", "black"),
                 cex = if (p > 10) 0.7 else 0.9)
  invisible(NULL)
}

#' Draw the posterior mass by SHD from the point estimate
#' @noRd
.plot_ball <- function(x, ...) {
  g <- x$graphs
  classes <- c("identical", "upper", "lower", "rearranged")
  cols <- c(identical = "black", upper = "#1b9e77", lower = "#d95f02",
            rearranged = "#7570b3")
  # Mass by SHD and class.
  shd_vals <- 0:max(g$shd, x$ball$radius)
  tab <- vapply(shd_vals, function(s) {
    vapply(classes, function(k) sum(g$mass[g$shd == s & g$class == k]), numeric(1))
  }, numeric(length(classes)))
  tab <- matrix(tab, nrow = length(classes), dimnames = list(classes, shd_vals))
  # Stacked bars and the radius of the ball.
  mids <- graphics::barplot(tab, col = cols[classes], border = NA,
                            xlab = "SHD from the point estimate",
                            ylab = "Posterior mass", ...)
  step <- if (length(mids) > 1L) mids[2L] - mids[1L] else 1
  graphics::abline(v = mids[x$ball$radius + 1L] + step / 2, lty = 2)
  graphics::legend("topright", legend = classes, fill = cols[classes], border = NA,
                   bty = "n")
  invisible(NULL)
}

#' Draw the trace of the score
#' @noRd
.plot_trace <- function(x, ...) {
  # The score of each retained draw against its iteration.
  tr <- x$trace
  graphics::plot(tr$iteration, tr$log_score, type = "l", xlab = "Iteration",
                 ylab = "Score", ...)
  invisible(NULL)
}

#' Draw the posterior probabilities of the edges of the point estimate
#' @noRd
.plot_edges <- function(x, ...) {
  # An empty plot with a title when there is no edge.
  e <- x$edges
  if (nrow(e) == 0L) {
    graphics::plot.new()
    graphics::title(main = "The point estimate has no edges")
    return(invisible(NULL))
  }
  # Edge labels.
  lab <- ifelse(e$type == "directed", paste0(e$from, " -> ", e$to),
                paste0(e$from, " - ", e$to))
  # Horizontal bars for the edge probabilities.
  op <- graphics::par(mar = c(4, max(6, 0.6 * max(nchar(lab))), 1, 1))
  on.exit(graphics::par(op))
  mids <- graphics::barplot(e$prob, names.arg = lab, horiz = TRUE, las = 1,
                            xlim = c(0, 1), col = "grey70", border = NA,
                            xlab = "Posterior probability", ...)
  # Points for the adjacency probabilities.
  graphics::points(e$prob_adjacent, mids, pch = 19)
  graphics::abline(v = 0.9, lty = 2)
  graphics::legend("bottomright", legend = c("edge", "adjacency"),
                   fill = c("grey70", NA), border = NA, pch = c(NA, 19),
                   bty = "n")
  invisible(NULL)
}
