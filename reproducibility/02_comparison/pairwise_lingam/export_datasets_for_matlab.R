# =============================================================================
# export_datasets_for_matlab.R   (Pairwise-LiNGAM, stage 1 of 3)
#
# Pairwise-LiNGAM (Hyvarinen and Smith 2013) runs in MATLAB with the
# authors' function pwling. This stage writes every dataset of the
# comparative study to a CSV file that MATLAB reads:
#
#     inputs/<experiment>_c<setting index>_r<dataset>_X.csv   (n x p, no header)
#
# and the list of datasets to manifest.csv, with the columns
#
#     experiment, cell_index, rep_id, graph_type, p, n, rho, effect_scale,
#     X_file, LR_file
#
# LR_file is the name of the file that MATLAB writes. The setting index is
# the row of comparison_cells(). The datasets are those that every other
# method sees: draw_dataset() with the regime tag "mixture".
#
# The three stages:
#     1. this script (R)                                    writes inputs/
#     2. ../../common/pairwise_lingam/run_pwling_batch.m     writes outputs/
#     3. import_pairwise_results.R (R)                      writes
#                                     ../results/results_pairwise_lingam.csv
#
# The commands of the three stages are in ../README.md.
#
# The files of inputs/ take about 1.9 GB. They are written again
# identically by this script, so they are not archived. Stage 3 reads them,
# so keep them until stage 3 has run.
#
# Run from any folder, for example from the root of the repository:
#     Rscript reproducibility/02_comparison/pairwise_lingam/export_datasets_for_matlab.R
# On a cluster, see slurm/export_pairwise_comparison.slurm.
# =============================================================================

# ---- Folder of this script and shared code ----------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- file.path(.here, "..", "..", "common")
source(file.path(COMMON_DIR, "scenario_instances.R"))
source(file.path(.here, "..", "comparison_cells.R"))

# ---- Write the datasets -------------------------------------------------------
# Create inputs/ and write one file per setting and dataset.
in_dir <- file.path(.here, "inputs")
dir.create(in_dir, showWarnings = FALSE, recursive = TRUE)

# The settings, and one row of the manifest for each setting and dataset.
cells <- comparison_cells()
man <- vector("list", nrow(cells) * N_REPS)
k <- 0L

# Write the datasets of each setting.
for (ci in seq_len(nrow(cells))) {
  cell <- cells[ci, ]
  inst <- build_instance(cell$graph_type, cell$p)
  view <- instance_view(inst, rho = cell$rho, effect_scale = cell$effect_scale)
  for (r in seq_len(N_REPS)) {
    # Write the dataset unless the file exists already.
    stem <- sprintf("%s_c%03d_r%03d", cell$experiment, ci, r)
    x_file <- paste0(stem, "_X.csv")
    x_path <- file.path(in_dir, x_file)
    if (!file.exists(x_path)) {
      X <- draw_dataset(view, n = cell$n, rep = r, regime_tag = "mixture")
      write.table(X, x_path, sep = ",", row.names = FALSE, col.names = FALSE)
    }
    # One row of the manifest.
    k <- k + 1L
    man[[k]] <- data.frame(experiment = cell$experiment, cell_index = ci,
                           rep_id = r, graph_type = cell$graph_type,
                           p = cell$p, n = cell$n, rho = cell$rho,
                           effect_scale = cell$effect_scale,
                           X_file = x_file, LR_file = paste0(stem, "_LR.csv"),
                           stringsAsFactors = FALSE)
  }
  cat(sprintf("Exported setting %d of %d (%s %s p = %d n = %d)\n", ci, nrow(cells),
              cell$experiment, cell$graph_type, cell$p, cell$n))
}

# ---- Write the manifest --------------------------------------------------------
# Check that every dataset is listed once and that its file exists.
manifest <- do.call(rbind, man[seq_len(k)])
stopifnot(nrow(manifest) == nrow(cells) * N_REPS,
          !anyDuplicated(manifest$LR_file),
          all(file.exists(file.path(in_dir, manifest$X_file))))
write.csv(manifest, file.path(.here, "manifest.csv"), row.names = FALSE)
cat(sprintf("Wrote manifest.csv with %d datasets.\n", nrow(manifest)))
