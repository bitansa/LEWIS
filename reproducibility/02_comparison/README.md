# 02_comparison: the comparative study

This folder reproduces the comparative study of the paper:

* Section 7.2: Figure 9 (mean SHD against p), Figure 10 (mean SHD against
  rho) and every number quoted in the section.
* Appendix G.1 (design and evaluation) and Appendix G.3: Figure 15 (mean
  SHD against n), Figure 16 (mean SHD against s), the SID results and the
  computing times of the existing methods.
* Appendix E.1, paragraph "Computing time": the run times of LEWIS.
* Section S1 of the online appendix: Tables S1 to S14 (mean and standard
  deviation of the SHD and of the two SID bounds in every setting).

## Design

There are 54 settings (`comparison_cells.R`): 26 core settings with
rho = 0.5 and s = 1, 12 settings that vary s and 16 settings that vary rho.
Each setting has 100 datasets. For each graph type and each p, the graph,
its coefficients and its error scales are drawn once and kept fixed
(`../common/scenario_instances.R`). The datasets of a setting differ only
in the errors. Every method sees the same datasets, because each dataset
is generated from a seed that depends only on the setting and the dataset
number.

Each estimate is compared with the true RPDAG by the SHD (`LEWIS::shd`),
and with the true DAG by the lower and upper SID bounds (package `SID`).
The evaluation rules are in `../common/method_wrappers.R`
(`est_partial_dag()`) and `../common/scoring.R`.

## Methods

| Method | Software | Settings (Appendix G.1) |
| --- | --- | --- |
| LEWIS (complexity) | `LEWIS::lewis()` | `d = 0.03`, `alpha = 0.75`, `center = FALSE`, start at the TL estimate, iterations from `lewis_schedule(p)` |
| LEWIS (uniform) | `LEWIS::lewis()` | `d = 0`, otherwise as above |
| PC | R package `pcalg` | Gaussian partial-correlation tests at level 0.05 |
| GES | R package `pcalg` | Gaussian BIC score |
| NOTEARS | Python package `notears` 3.0, `notears.linear.notears_linear` | least-squares loss, l1 penalty 0.1, edge threshold 0.3 |
| ICA-LiNGAM | Python package `lingam` 1.12.2, `ICALiNGAM` | default settings (adaptive lasso pruning) |
| DirectLiNGAM | Python package `lingam` 1.12.2, `DirectLiNGAM` | default settings (adaptive lasso pruning) |
| RCD | Python package `lingam` 1.12.2, `RCD` | default settings. Pairs flagged as confounded are non-adjacent |
| TL | `LEWIS::tl_lingam()`, a copy of `TLLiNGAM()` of the R package TransGraph 1.1.0 | default settings |
| Pairwise-LiNGAM | MATLAB code `pwling` of the authors | skeleton by partial-correlation tests at level 0.05, each edge oriented by the sign of the likelihood ratio |

The wrappers are in `../common/method_wrappers.R`. The R packages, the
Python packages and their versions are listed once, in `../README.md`. The
Python methods are called from R through `reticulate`. Set the environment
variable `RETICULATE_PYTHON` to the Python that has `lingam` and `notears`.
Pairwise-LiNGAM is described in the last section of this README.

## Files

| File | Purpose |
| --- | --- |
| `comparison_cells.R` | The 54 settings and the number of datasets (100) |
| `run_one_method.R` | Runs one method on all settings. Writes `results/results_<method>.csv` |
| `pairwise_lingam/export_datasets_for_matlab.R` | Pairwise-LiNGAM stage 1: writes the datasets to `pairwise_lingam/inputs/` and the list `pairwise_lingam/manifest.csv` |
| `../common/pairwise_lingam/run_pwling_batch.m` | Pairwise-LiNGAM stage 2 (MATLAB): writes `pairwise_lingam/outputs/` |
| `pairwise_lingam/import_pairwise_results.R` | Pairwise-LiNGAM stage 3: writes `results/results_pairwise_lingam.csv` |
| `aggregate_and_plot.R` | Checks the results and writes the summary, the figures and the tables |
| `measure_lewis_runtime.R` | Run times of LEWIS for Appendix E.1. Writes `results/lewis_runtime_complexity.csv` |
| `variance_by_scale.R` | Variances and correlations of the fixed graphs for each s (the correlations up to 0.999 on SF graphs at s = 2 in Section 7.2) |

Every R script finds its own folder, so it can be started from any
folder. The commands below are given from the root of the repository (the
folder that contains `reproducibility/`).

## Run order

The computations are cluster jobs. Each SLURM file in
`reproducibility/slurm/` only sets the resources, loads R (and Python) and
calls one R script. The same step can therefore be run without SLURM by
the plain `Rscript` command shown. `MC_CORES` sets the number of datasets
that run at the same time (default: all cores but one). Before the first
SLURM job, run `mkdir -p reproducibility/slurm/logs`, and edit the
partition, modules and resources of each SLURM file for your cluster.

1. One job per method. They are independent and can run at the same time.

   ```sh
   sbatch reproducibility/slurm/comparison_lewis_complexity.slurm
   sbatch reproducibility/slurm/comparison_lewis_uniform.slurm
   sbatch reproducibility/slurm/comparison_pc.slurm
   sbatch reproducibility/slurm/comparison_ges.slurm
   sbatch reproducibility/slurm/comparison_lingam.slurm
   sbatch reproducibility/slurm/comparison_direct_lingam.slurm
   sbatch reproducibility/slurm/comparison_rcd.slurm
   sbatch reproducibility/slurm/comparison_notears.slurm
   sbatch reproducibility/slurm/comparison_tl_zhao.slurm
   ```

   Without SLURM, run for each method

   ```sh
   export RETICULATE_PYTHON=/path/to/python3   # for lingam, direct_lingam, rcd, notears
   MC_CORES=100 Rscript reproducibility/02_comparison/run_one_method.R <method>
   ```

   with `<method>` one of `lewis_complexity`, `lewis_uniform`, `pc`, `ges`,
   `lingam`, `direct_lingam`, `rcd`, `notears`, `tl_zhao`.

   Each run writes `02_comparison/results/results_<method>.csv` (5400
   rows). A run that is started again continues where it stopped, and runs
   again the datasets on which the method failed. A run exits with status
   1 if a setting has fewer than 100 datasets or if the method failed on a
   dataset (rows with `skipped = TRUE`). The error messages are printed
   (with SLURM, in `reproducibility/slurm/logs/`).

2. Pairwise-LiNGAM, stage 1. It writes the 5400 datasets (about 1.9 GB) to
   `pairwise_lingam/inputs/` and the list `pairwise_lingam/manifest.csv`.

   ```sh
   sbatch reproducibility/slurm/export_pairwise_comparison.slurm
   # or, without SLURM:
   Rscript reproducibility/02_comparison/pairwise_lingam/export_datasets_for_matlab.R
   ```

3. Pairwise-LiNGAM, stage 2, in MATLAB. If MATLAB runs on another machine,
   copy `pairwise_lingam/inputs/` and `pairwise_lingam/manifest.csv` there
   first, and copy `pairwise_lingam/outputs/` back afterwards. Then, in a
   terminal:

   ```sh
   export PWLING_PATH=/path/to/pwling
   cd reproducibility/02_comparison/pairwise_lingam
   matlab -batch "addpath('../../common/pairwise_lingam'), run_pwling_batch"
   ```

   The script uses all cores but one if the Parallel Computing Toolbox is
   installed, and runs serially otherwise. It can be restarted and skips
   the datasets that are done. It stops with an error if a dataset failed.

4. Pairwise-LiNGAM, stage 3. It needs `inputs/` and `manifest.csv`,
   because the skeleton is computed from the exported datasets.

   ```sh
   sbatch reproducibility/slurm/import_pairwise_comparison.slurm
   # or, without SLURM:
   MC_CORES=100 Rscript reproducibility/02_comparison/pairwise_lingam/import_pairwise_results.R
   ```

5. Run times of LEWIS (Appendix E.1). The script also checks that the SHD
   of every run equals the SHD in `results/results_lewis_complexity.csv`
   and prints "SHD agrees with the comparative study on 1200 of 1200
   datasets."

   ```sh
   sbatch reproducibility/slurm/lewis_runtime.slurm
   # or, without SLURM:
   MC_CORES=100 Rscript reproducibility/02_comparison/measure_lewis_runtime.R
   ```

6. Summary, figures and tables, on any machine (a few minutes):

   ```sh
   Rscript reproducibility/02_comparison/aggregate_and_plot.R
   ```

   It stops unless all ten result files are complete and prints "Checks
   passed: 10 methods, 54 settings, 100 datasets each, no missing value."
   It writes `results/comparison_all.csv`,
   `results/comparison_cell_summary.csv`, and, in
   `<figure folder>/comparison/`:

   * `comparison_shd_vs_p` (Figure 9), `comparison_shd_vs_rho` (Figure
     10), `comparison_shd_vs_n` (Figure 15) and `comparison_shd_vs_scale`
     (Figure 16), each as PDF and PNG.
   * in `tables/`: `shd_core_ER.tex`, `shd_core_SF.tex`,
     `shd_mixture_ER.tex`, `shd_mixture_SF.tex`, `shd_effect_ER.tex`,
     `shd_effect_SF.tex`, `sid_lower_core_ER.tex`, `sid_upper_core_ER.tex`,
     `sid_lower_core_SF.tex`, `sid_upper_core_SF.tex`,
     `sid_mixture_ER.tex`, `sid_mixture_SF.tex`, `sid_effect_ER.tex` and
     `sid_effect_SF.tex` (Tables S1 to S14, in this order).

   The figure folder is `reproducibility/figures/` unless the environment
   variable `LEWIS_FIG_DIR` is set.

7. The numbers quoted in the text are printed by
   `Rscript reproducibility/paper_numbers.R`. The correlations of Section
   7.2 are printed by `Rscript reproducibility/02_comparison/variance_by_scale.R`
   (no data, a few seconds).

Optional check while the jobs run:
`Rscript reproducibility/check_progress.R comparison`.

## Results in this folder and on Zenodo

`results/` holds the small files that the paper uses:

* `comparison_cell_summary.csv`: one row per setting and method (540
  rows) with the mean, standard deviation, median and quartiles of the SHD
  (`shd_mean`, ...), the mean, standard deviation and median of the lower
  and upper SID bounds (`sid_lo_mean`, `sid_up_mean`, ...) and the mean
  run time.
* `lewis_runtime_complexity.csv`: one row per dataset of the run-time
  measurement (1200 rows). `ms_per_iter` is the time per iteration net of
  the TL initialization, `1000 * (runtime_s - warm_start_s) / n_iter`. The
  column `engine` names the software that ran LEWIS. In the stored file it
  reads `rcpp`, the compiled code that the package LEWIS contains. A new
  run writes the name and version of the package, for example
  `LEWIS 0.1.0`.

The results of each dataset are in the archive
`lewis_raw_02_comparison.tar.gz` on Zenodo (see `../ZENODO.md`). Unpack it
in `reproducibility/`. It contains:

* `results/results_<method>.csv` for the ten methods (5400 rows each). The
  summary, the figures, the tables and `paper_numbers.R` need them.
* `pairwise_lingam/manifest.csv` and `pairwise_lingam/outputs/` (the
  likelihood-ratio matrices of MATLAB and `runtimes.csv`). They are needed
  only to run stage 3 again without MATLAB. Stage 3 also needs `inputs/`,
  which stage 1 writes again, identically.

`results/comparison_all.csv` is not archived. `aggregate_and_plot.R`
writes it.

## Hardware and run time

The jobs ran on a cluster node with 100 cores and 128 GB of memory, one
dataset per core, in R 4.3.2. The total CPU time of each method over the
5400 datasets, without the SID, was about 400 hours for LEWIS (complexity),
825 hours for LEWIS (uniform), 300 hours for RCD, 120 hours for NOTEARS,
80 hours for TL, 40 hours for DirectLiNGAM, and less than 3 hours for each
of PC, GES, ICA-LiNGAM and Pairwise-LiNGAM. The SID adds hours at p = 50
and p = 100. The times of LEWIS in `results_lewis_*.csv` were measured
with an implementation of LEWIS in R only. The package LEWIS uses compiled
code and is several times faster.

The run-time measurement of step 5 (1200 calls of `lewis()` and 1200
separate TL calls) takes about 70 CPU hours. Stage 1 takes about an hour
on one core. Stage 2 takes minutes. Stage 3 takes a few hours on 100
cores, because of the SID.

## Pairwise-LiNGAM

Pairwise-LiNGAM uses the MATLAB function `pwling` of

> A. Hyvarinen and S. M. Smith (2013). Pairwise likelihood ratios for
> estimation of non-Gaussian structural equation models. Journal of
> Machine Learning Research, 14, 111-152.

This code is not part of the repository. The section "Software" of
`../README.md` explains where to download it and how to set the
environment variable `PWLING_PATH`.

The study runs in three stages, because MATLAB may run on another
machine than R:

1. R writes each dataset to a CSV file (`n` rows, `p` columns, no header)
   in `pairwise_lingam/inputs/`, and lists the datasets in
   `pairwise_lingam/manifest.csv`.
2. MATLAB (`../common/pairwise_lingam/run_pwling_batch.m`) runs
   `LR = pwling(X', 1)` on each dataset and writes the `p` x `p` matrix
   `LR` to `pairwise_lingam/outputs/` (CSV, no header), with the run times
   in `outputs/runtimes.csv`.
3. R builds the skeleton from the dataset by Fisher z tests of the partial
   correlations at level 0.05, and orients each edge {i, j} as i -> j if
   `LR[i, j] > 0` and j -> i if `LR[i, j] < 0`. The graph is used as
   returned, so it may contain a directed cycle. An exact zero leaves the
   edge undirected.

The same three stages, with the same MATLAB script, produce the
Pairwise-LiNGAM results of Figure 1 (see `../01_motivation/README.md`).
