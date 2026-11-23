# 01_motivation: Figure 1

This folder reproduces Figure 1 of the paper (Section 1). The design is in
Appendix G.1, paragraph "Design of Figure 1".

Figure 1 shows the mean SHD over 100 datasets between the estimate of each
method and the true RPDAG, on the 12-node DAG of Figure 7(a), as the errors
become closer to Gaussian. In both rows, the first k nodes of the order
1, 2, 3, 5, 10, 4, 6, 7, 8, 9, 11, 12 have Student's t5 errors, for
k = 12, 11, ..., 0. In the top row the other nodes have Gaussian errors. In
the bottom row they have t15 errors. The columns are n = 100, 250 and 500.

## Regimes

The study has 25 regimes, defined in `motivation_config.R`. With the three
sample sizes, this gives 75 settings of 100 datasets each, 7500 datasets in
total.

| Regime | Errors |
| --- | --- |
| `all_t5` | t5 on all 12 nodes (k = 12, both rows) |
| `count_t5g_k01` to `count_t5g_k11` | t5 on the first k nodes, Gaussian on the others (top row) |
| `gaussian` | Gaussian on all nodes (k = 0, top row) |
| `count_t5t15_k01` to `count_t5t15_k11` | t5 on the first k nodes, t15 on the others (bottom row) |
| `t15` | t15 on all nodes (k = 0, bottom row) |

The seed of a dataset depends only on the regime name, n and the dataset
number. The regime names are part of the seeds and must not be renamed.

## Methods

LEWIS (complexity), LEWIS (uniform), PC, GES, NOTEARS, ICA-LiNGAM,
DirectLiNGAM, RCD and TL run in R (`motivation_experiment.R`).
Pairwise-LiNGAM runs in MATLAB on the same datasets (`pairwise_lingam/`).
The settings of the methods and the evaluation are those of Section 7.2
(see `../02_comparison/README.md`). Of the two versions of LEWIS, Figure 1
shows only LEWIS (complexity), with the legend entry "LEWIS". The rows of
LEWIS (uniform) are kept in the results. The SID is not computed in this
study.

The nine methods in R run one after the other on each dataset, in the
order of `motivation_experiment.R`, and share one random stream that
starts at the seed of the dataset. The methods are not seeded one by one,
because a separate seed for each method would change the stored results.

## Software

R, the package LEWIS, the R packages and the Python packages (`lingam`
1.12.2 and `notears` 3.0, called through `reticulate`) are listed in
`../README.md`. Set `RETICULATE_PYTHON` to the Python that has `lingam`
and `notears`. Pairwise-LiNGAM needs MATLAB and the authors' code `pwling`.
Where to get it and how to set `PWLING_PATH` is explained in the section
"Software" of `../README.md`. The three stages of Pairwise-LiNGAM are
described in the section "Pairwise-LiNGAM" of
`../02_comparison/README.md`.

## Files

| File | Purpose |
| --- | --- |
| `motivation_config.R` | The 25 regimes, the rows of Figure 1 and the data generation |
| `motivation_experiment.R` | Runs the nine methods in R on all datasets. Writes `results/motivation_results.csv` |
| `pairwise_lingam/export_motivation_for_matlab.R` | Pairwise-LiNGAM stage 1: writes the datasets to `pairwise_lingam/inputs/` and the list `pairwise_lingam/manifest.csv` |
| `../common/pairwise_lingam/run_pwling_batch.m` | Pairwise-LiNGAM stage 2 (MATLAB): writes the likelihood-ratio matrices to `pairwise_lingam/outputs/` |
| `pairwise_lingam/import_motivation_pairwise.R` | Pairwise-LiNGAM stage 3: adds the Pairwise-LiNGAM rows to `results/motivation_results.csv` |
| `plot_motivation.R` | Draws Figure 1 and writes the two summary tables |
| `results/motivation_cell_summary.csv` | Mean, standard deviation, median and quartiles of the SHD for each regime, n and method (750 rows) |
| `results/motivation_gaussianity_table.csv` | Median Kolmogorov-Smirnov p-value and mean excess kurtosis of the errors for each regime and n (75 rows) |

Every R script finds its own folder, so it can be started from any
folder. The commands below are given from the root of the repository (the
folder that contains `reproducibility/`).

## Run order

Steps 1, 2 and 4 are cluster jobs. Each SLURM file only calls the R
script of the step, so the same step can be run without SLURM by the
plain `Rscript` command shown. `MC_CORES` sets the number of datasets that
run at the same time (default: all cores but one).

1. The nine methods in R (about 80 CPU hours, 1 to 2 hours on 100 cores):

   ```sh
   mkdir -p reproducibility/slurm/logs
   sbatch reproducibility/slurm/motivation.slurm
   # or, without SLURM:
   export RETICULATE_PYTHON=/path/to/python3
   MC_CORES=100 Rscript reproducibility/01_motivation/motivation_experiment.R
   ```

   The script can be restarted. It skips every dataset whose nine rows are
   present without failure. A dataset with a failed method is run again
   from the start, so that the random stream is the same.

2. Pairwise-LiNGAM, stage 1 (one core, a few minutes, about 460 MB of CSV
   files):

   ```sh
   sbatch reproducibility/slurm/motivation_pairwise_export.slurm
   # or, without SLURM:
   Rscript reproducibility/01_motivation/pairwise_lingam/export_motivation_for_matlab.R
   ```

3. Pairwise-LiNGAM, stage 2, in MATLAB (a few minutes). If MATLAB runs on
   another machine, copy `pairwise_lingam/inputs/` and
   `pairwise_lingam/manifest.csv` there first, and copy
   `pairwise_lingam/outputs/` back afterwards. Then, in a terminal:

   ```sh
   export PWLING_PATH=/path/to/pwling
   cd reproducibility/01_motivation/pairwise_lingam
   matlab -batch "addpath('../../common/pairwise_lingam'), run_pwling_batch"
   ```

4. Pairwise-LiNGAM, stage 3 (one core, under an hour). It reads
   `inputs/` to build the skeleton, so keep `inputs/` until this step is
   done. It needs the rows of step 1.

   ```sh
   sbatch reproducibility/slurm/motivation_pairwise_import.slurm
   # or, without SLURM:
   Rscript reproducibility/01_motivation/pairwise_lingam/import_motivation_pairwise.R
   ```

5. Figure 1 and the two tables, on any machine (a minute):

   ```sh
   Rscript reproducibility/01_motivation/plot_motivation.R
   ```

   It writes `<figure folder>/motivation/motivation_intro_2x3.pdf` and
   `.png` (Figure 1), `results/motivation_cell_summary.csv` and
   `results/motivation_gaussianity_table.csv`, and prints the mean SHD of
   every panel. The figure folder is `reproducibility/figures/` unless the
   environment variable `LEWIS_FIG_DIR` is set.

Every script stops, or exits with a non-zero status, if a dataset or a
method fails. It prints the error messages.

## Raw results from Zenodo

To draw Figure 1 without running the study, get the archive
`lewis_raw_01_motivation.tar.gz` (see `../ZENODO.md`) and unpack it in
`reproducibility/`. It contains:

| File | Content |
| --- | --- |
| `results/motivation_results.csv` | One row per regime, n, method and dataset (75,000 rows) |
| `pairwise_lingam/outputs/*_LR.csv` and `pairwise_lingam/outputs/runtimes.csv` | The likelihood-ratio matrices of pwling for the 7500 datasets |

Only `results/motivation_results.csv` is needed for step 5. The
likelihood-ratio matrices are needed only to run step 4 again without
MATLAB. In that case, run step 2 first (it writes `inputs/` and
`manifest.csv` again, identically), and remove the rows with
`method = pairwise_lingam` from `results/motivation_results.csv`, because
step 4 skips the rows that are present.

## Hardware and run time

The study ran on a cluster node with 100 cores and 128 GB of memory, one
dataset per core, in R 4.3.2. Step 1 took about 80 CPU hours in total,
almost all of it for the two versions of LEWIS. These times were measured
with an implementation of LEWIS in R only. The package LEWIS uses compiled
code and is faster. Steps 2 to 5 take minutes.
