# =============================================================================
# export_motivation_for_matlab.R   (Figure 1, Pairwise-LiNGAM stage 1 of 3)
#
# Pairwise-LiNGAM (Hyvarinen and Smith 2013) runs in MATLAB with the
# authors' code pwling. This stage writes every dataset of Figure 1 to a
# CSV file, so that MATLAB runs on the same datasets as the other methods:
#
#     inputs/<regime>_n<n>_r<rep>_X.csv     (n x 12, no header)
#
# and lists them in manifest.csv, one row per (regime, n, dataset), with
# the columns
#
#     experiment, cell_index, rep_id, graph_type, p, n, rho, effect_scale,
#     X_file, LR_file
#
# Here `experiment` is the regime name, `cell_index` is the position of the
# regime in motivation_config.R, `graph_type` is "paper12" and `rho` is
# k / 12. `LR_file` is the name of the output file of MATLAB. A dataset is
# identified by LR_file (or by experiment, n and rep_id). cell_index alone
# does not identify it, because it is shared by the three values of n.
#
# The datasets are generated from their seeds, so they are the same on
# every run and every machine with the same R version. Run this script with
# the R installation that ran motivation_experiment.R.
#
# The three stages:
#     1. this script (R, writes inputs/ and manifest.csv)
#     2. common/pairwise_lingam/run_pwling_batch.m (MATLAB, run from this
#        folder, writes outputs/)
#     3. import_motivation_pairwise.R (R, evaluates the estimates and adds
#        them to ../results/motivation_results.csv)
# The commands of the three stages are in ../README.md.
#
# Run from any folder, for example from the root of the repository:
#     Rscript reproducibility/01_motivation/pairwise_lingam/export_motivation_for_matlab.R
# On a cluster, see slurm/motivation_pairwise_export.slurm.
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
source(file.path(COMMON_DIR, "paper12_instance.R"))
source(file.path(.here, "..", "motivation_config.R"))

# ---- Write the datasets -------------------------------------------------------
# Create inputs/ and write one file per regime, n and dataset.
in_dir <- file.path(.here, "inputs")
dir.create(in_dir, showWarnings = FALSE, recursive = TRUE)

# One row of the manifest for each regime, n and dataset.
man <- vector("list", length(REGIMES) * length(N_VALUES) * N_REPS)
k <- 0L
for (n in N_VALUES) {
  for (ri in seq_along(REGIMES)) {
    rg <- REGIMES[[ri]]
    for (r in seq_len(N_REPS)) {
      # File names of this dataset.
      stem <- sprintf("%s_n%d_r%03d", rg$regime, n, r)
      x_file <- paste0(stem, "_X.csv")
      x_path <- file.path(in_dir, x_file)
      # Write the dataset unless it exists already.
      if (!file.exists(x_path)) {
        X <- draw_motivation_dataset(rg, n = n, rep = r)
        write.table(X, x_path, sep = ",", row.names = FALSE, col.names = FALSE)
      }
      # Add its row to the manifest.
      k <- k + 1L
      man[[k]] <- data.frame(experiment = rg$regime, cell_index = ri,
                             rep_id = r, graph_type = "paper12",
                             p = P12, n = n,
                             rho = rg$k_ng / P12, effect_scale = 1.0,
                             X_file = x_file,
                             LR_file = paste0(stem, "_LR.csv"),
                             stringsAsFactors = FALSE)
    }
    cat(sprintf("Exported n = %d, regime %d of %d (%s)\n", n, ri,
                length(REGIMES), rg$regime))
  }
}

# ---- Write the manifest --------------------------------------------------------
# Check that every dataset is listed once and that its file exists.
manifest <- do.call(rbind, man[seq_len(k)])
stopifnot(nrow(manifest) == length(REGIMES) * length(N_VALUES) * N_REPS,
          !anyDuplicated(manifest$LR_file),
          all(file.exists(file.path(in_dir, manifest$X_file))))
write.csv(manifest, file.path(.here, "manifest.csv"), row.names = FALSE)
cat(sprintf("Wrote manifest.csv with %d datasets.\n", nrow(manifest)))
