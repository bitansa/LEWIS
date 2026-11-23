# =============================================================================
# paths.R
#
# Output locations shared by all scripts.
#
# The figures and LaTeX tables of the paper are written to
#
#     <figure folder>/<study subfolder>/
#
# The figure folder is the value of the environment variable LEWIS_FIG_DIR.
# If it is not set, the folder is reproducibility/figures/. For example, to
# write the figures straight into the folder of the paper, run
#
#     export LEWIS_FIG_DIR=/path/to/paper/figures
#
# before the plot scripts. The study subfolders motivation/, comparison/,
# comparison/tables/, sensitivity/, uq/ and real_data/TCGA/ are the ones
# that the paper reads. The subfolder laplace_vs_exact/ holds the body of
# Table 1, which the paper types in.
#
# The calling script must define COMMON_DIR, the path of this folder,
# before it sources this file.
# =============================================================================

stopifnot(exists("COMMON_DIR"))

# The root of the reproducibility folder.
REPRO_DIR <- normalizePath(file.path(COMMON_DIR, ".."), mustWork = TRUE)

# Return the figure folder of one study, and create it if needed.
lewis_fig_dir <- function(...) {
  # Read the environment variable, and fall back to reproducibility/figures.
  root <- Sys.getenv("LEWIS_FIG_DIR", "")
  if (!nzchar(root)) root <- file.path(REPRO_DIR, "figures")
  # Build and create the subfolder.
  d <- normalizePath(file.path(root, ...), mustWork = FALSE)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}
