# =============================================================================
# comparison_cells.R
#
# The 54 settings of the comparative study (Section 7.2 and Appendix G.1).
# Each setting is a graph type and values of (p, n, rho, s). All settings
# use the fixed instances of common/scenario_instances.R.
#
#   core          rho = 0.5 and s = 1, with the 13 pairs (p, n)
#                     p = 10  : n in {50, 100, 500, 1000}
#                     p = 20  : n in {100, 250, 500, 1000}
#                     p = 50  : n in {250, 500, 1000}
#                     p = 100 : n in {500, 1000}
#                 and both graph types (26 settings).
#   effect_scale  s in {0.25, 0.5, 2} at (p, n) in {(20, 500), (50, 500)},
#                 rho = 0.5, both graph types (12 settings). The setting
#                 with s = 1 is a core setting.
#   mixture       rho in {0, 0.25, 0.75, 1} at (p, n) in {(20, 500),
#                 (50, 500)}, s = 1, both graph types (16 settings). The
#                 setting with rho = 0.5 is a core setting.
#
# Each setting has N_REPS = 100 datasets, which differ only in the errors.
#
# This file is sourced by the scripts of the study. It does not run
# anything.
# =============================================================================

# Number of datasets per setting.
N_REPS <- 100L

# The 54 settings as a data frame with one row per setting.
comparison_cells <- function() {
  # The core settings.
  core_pn <- rbind(
    expand.grid(p = 10L,  n = c(50L, 100L, 500L, 1000L)),
    expand.grid(p = 20L,  n = c(100L, 250L, 500L, 1000L)),
    expand.grid(p = 50L,  n = c(250L, 500L, 1000L)),
    expand.grid(p = 100L, n = c(500L, 1000L))
  )
  core <- merge(expand.grid(graph_type = c("ER", "SF"), stringsAsFactors = FALSE),
                core_pn)
  core$rho <- 0.5
  core$effect_scale <- 1.0
  core$experiment <- "core"

  # The settings that vary the effect-size scale s.
  es <- expand.grid(graph_type = c("ER", "SF"), p = c(20L, 50L),
                    effect_scale = c(0.25, 0.5, 2.0), stringsAsFactors = FALSE)
  es$n <- 500L
  es$rho <- 0.5
  es$experiment <- "effect_scale"

  # The settings that vary the fraction rho of non-Gaussian errors.
  mx <- expand.grid(graph_type = c("ER", "SF"), p = c(20L, 50L),
                    rho = c(0.0, 0.25, 0.75, 1.0), stringsAsFactors = FALSE)
  mx$n <- 500L
  mx$effect_scale <- 1.0
  mx$experiment <- "mixture"

  cols <- c("experiment", "graph_type", "p", "n", "rho", "effect_scale")
  cells <- rbind(core[, cols], es[, cols], mx[, cols])

  # Every setting has n >= 5 p.
  stopifnot(nrow(cells) == 54L, all(cells$n >= 5L * cells$p))

  # Order the settings from the cheapest to the most expensive. The order
  # also defines the setting index used in the file names of the
  # Pairwise-LiNGAM datasets, so it must not change.
  cells[order(cells$p, cells$n, cells$experiment, cells$graph_type), , drop = FALSE]
}

# Key of one setting, formatted with sprintf("%g") so that it is the same
# whether it is built from the grid or read back from a CSV file.
cell_key <- function(experiment, graph_type, p, n, rho, effect_scale) {
  sprintf("%s|%s|p%d|n%d|rho%g|s%g", experiment, graph_type,
          as.integer(p), as.integer(n), rho, effect_scale)
}
