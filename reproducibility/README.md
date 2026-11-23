# Reproducibility of the paper

This folder reproduces the simulations, the real-data analysis, the figures
and the tables of the paper

> A. Chaudhuri, A. Bhattacharya, Y. Ni and B. Sarkar. Consistent DAG
> selection for Bayesian causal discovery under general error
> distributions.

with the method LEWIS of the R package `LEWIS` in this repository. Every
study calls `LEWIS::lewis()` with `center = FALSE`, as in the paper. The
folder is not part of the R package build.

## Where everything is

- **This folder (GitHub).** Every file of this folder is in the
  repository: all the code, and the small summary files from which the
  numbers, tables and figures of the paper are drawn. Each of these files
  is smaller than 5 MB.
- **Zenodo.** The results of each dataset of the simulations are too large
  for the repository. They are deposited on Zenodo. `ZENODO.md` lists the
  files and explains how to unpack them so that the scripts find them.
  `make_zenodo_archive.sh` is the script that made the archive.
- **TCGA data.** The real-data analysis uses public data of The Cancer
  Genome Atlas (TCGA). They are not in this repository and not in the
  Zenodo archive. `05_real_data/download_tcga.R` downloads them and checks
  their md5 sums. The preprocessing is done by `05_real_data/prepare_tcga.R`,
  which `05_real_data/run_reported_analyses.R` calls.

## What each folder reproduces

Each folder has its own README with the design, the inputs, the run order,
the commands, the outputs and the run time.

| Folder | Paper |
| --- | --- |
| `01_motivation/` | Figure 1 (Section 1). Its design is in Appendix G.1, paragraph "Design of Figure 1". |
| `02_comparison/` | Section 7.2 (comparative study) with Figures 9 and 10, Appendix G.1 and G.3 with Figures 15 and 16, the computing times of Appendix E.1, and Tables S1 to S14 of the online appendix. |
| `03_calibration/` | The remark "Choice of alpha and d" in Section 5, Appendix G.2 with Figure 14, and Figures S1 to S4 of the online appendix. |
| `04_uq/` | Section 7.3 (uncertainty quantification) with Figure 11(b), and Appendix H with Figure 17. |
| `05_real_data/` | Section 7.4 (real-data analysis of TCGA gene expression data), Appendix I with Table 2 and Figures 18 and 19, and Table S15 of the online appendix. Figure 12 is drawn in LaTeX from the tables of this folder. |
| `06_laplace_vs_exact/` | Appendix F (accuracy of the Laplace approximation) and its Table 1. |
| `common/` | Code shared by the studies: data generation, method wrappers, evaluation, figure style, the two chains of Appendix F, the MATLAB stage of Pairwise-LiNGAM, and the record of the software versions. |
| `slurm/` | One SLURM file for each cluster job. |

Files at the top of this folder:

| File | Purpose |
| --- | --- |
| `smoke_test.R` | Checks in a few minutes that the package and the shared code are installed and reproduce stored results of the paper. Run it first. |
| `paper_numbers.R` | Prints the numbers quoted in the text, computed from the saved results. |
| `check_progress.R` | Shows how far the simulation jobs have come. |
| `ZENODO.md` | The content of the Zenodo archive and how to unpack it. |
| `make_zenodo_archive.sh` | Builds the Zenodo archive from a copy of this folder in which the studies were run. |

## Software

All the software needed is listed here.

**R and the package.** Install R (4.1 or later) and a C++ compiler (the
Xcode command line tools on macOS, gcc on Linux). First install the R
packages that LEWIS needs. Then install LEWIS from the root of the
repository:

```sh
Rscript -e 'install.packages(c("Rcpp", "MASS", "huge", "dcov", "remotes"))'
R CMD INSTALL .
```

The package can also be installed from GitHub with
`remotes::install_github("bitansa/LEWIS")`.

**R packages.** The studies also need these R packages:

```r
install.packages(c("igraph", "SID", "reticulate", "ggplot2",
                   "patchwork", "dplyr", "tidyr", "scales", "data.table",
                   "R.utils", "quantreg"))
install.packages("BiocManager")
BiocManager::install(c("graph", "RBGL"))
install.packages("pcalg")
```

If a package is no longer on CRAN, install the version of the paper from
the CRAN archive with the package `remotes`, for example
`remotes::install_version("SID", "1.1")`.

**Python packages.** ICA-LiNGAM, DirectLiNGAM and RCD come from the Python
package `lingam`. NOTEARS comes from the authors' Python package `notears`
(version 3.0, module `notears.linear`). Both are called from R through
`reticulate`. Install them in one Python environment:

```sh
pip install lingam==1.12.2 numpy==1.26.4 scipy==1.13.0 scikit-learn==1.4.1.post1
pip install git+https://github.com/xunzheng/notears.git
```

Then set the environment variable `RETICULATE_PYTHON` to the Python of
this environment, for example
`export RETICULATE_PYTHON="$(command -v python3)"`. The wrappers are in
`common/method_wrappers.R`. The settings of each method are in
`02_comparison/README.md`.

**MATLAB and pwling.** Pairwise-LiNGAM runs in MATLAB with the function
`pwling` of Hyvarinen and Smith (2013, Journal of Machine Learning
Research, 14, 111-152). It is not part of this repository. Download it from
the authors' code page, https://www.cs.helsinki.fi/u/ahyvarin/code/pwcausal/,
as the file https://www.cs.helsinki.fi/u/ahyvarin/code/pwcausal/pwling.zip.
Unpack it into a folder of your choice, and set the environment variable
`PWLING_PATH` to the folder that contains `pwling.m` before you start
MATLAB. The three stages of Pairwise-LiNGAM (R, MATLAB, R) are described
in `02_comparison/README.md`.

**Operating system.** The scripts run on macOS and Linux. They run the
datasets in parallel with `parallel::mclapply`, which does not run in
parallel on Windows.

### Versions used for the paper

The simulations were run in R 4.3.2 on a computing cluster. The real-data
analysis was run in R 4.6.1 on a laptop. The R packages were `Rcpp` 1.1.2,
`quantreg` 6.1, `igraph` 2.3.3, `pcalg` 2.7.12, `TransGraph` 1.1.0, `SID`
1.1 and `reticulate` 1.46.0. The Python methods ran under Python 3.11.5
with `lingam` 1.12.2, `notears` 3.0, `numpy` 1.26.4, `scipy` 1.13.0 and
`scikit-learn` 1.4.1.post1. Pairwise-LiNGAM ran under MATLAB R2026a.

TL was run with the function `TLLiNGAM()` of `TransGraph` 1.1.0. This
package is archived on CRAN. The function `LEWIS::tl_lingam()` is a copy
of its TL code, and the scripts of this folder use it.

`common/record_versions.R` writes the versions of a machine to a file.
Run it in the environment of the jobs, so that `reticulate` finds the
same Python:

```sh
export RETICULATE_PYTHON="$(command -v python3)"
Rscript reproducibility/common/record_versions.R
```

This writes `common/versions.txt`.

## Before any run

```sh
Rscript reproducibility/smoke_test.R
Rscript reproducibility/common/check_orientation_conventions.R
```

The smoke test ends with "All checks passed." The second script needs the
Python packages. It checks that every method wrapper returns its graph
with the parent in the row, and it ends with "No wrapper returns reversed
edges."

## Three ways to use this folder

1. **Check the numbers of the paper.** Unpack the Zenodo archive in
   `reproducibility/` (see `ZENODO.md`) and run
   `Rscript reproducibility/paper_numbers.R`. Without the archive, the
   script still prints every number that comes from the summary files of
   this folder, and names the blocks that need the archive.
2. **Draw the figures and tables again without the simulations.** Unpack
   the Zenodo archive and run the scripts of the column "Figures and
   tables" of the table below. They take minutes on a laptop.
3. **Run everything again.** Run the steps of the table below in order.
   The simulations need a computing cluster or a large machine (see
   "Hardware and run times").

## Run order

All commands start from the root of the repository (the folder that
contains `reproducibility/`). Every R script finds its own folder, so the
commands also work from other folders.

**With SLURM.** Each file ending in `.slurm` in the column "Computation"
is in `slurm/`. Create the folder of the logs once with
`mkdir -p reproducibility/slurm/logs`. Edit the partition, modules and
resources of each SLURM file for your cluster. Then submit it, for
example with `sbatch reproducibility/slurm/motivation.slurm`. The logs are
written to `reproducibility/slurm/logs/`.

**Without SLURM.** Each SLURM file only sets the resources, loads R (and
Python) and runs one R script. The header of each SLURM file gives the
same command as a plain `Rscript` command. For example:

```sh
MC_CORES=8 Rscript reproducibility/02_comparison/run_one_method.R lewis_complexity
```

`MC_CORES` is the number of datasets that run at the same time. Without
it, all cores but one are used. The jobs can be stopped and started again.
They continue where they stopped. Each job exits with a non-zero status if
a dataset fails or a result is missing.

The steps of one study run from left to right. In the column "Output", a
path that starts with `motivation/`, `comparison/`, `sensitivity/`, `uq/`,
`real_data/` or `laplace_vs_exact/` is in the figure folder (see "Where
the outputs go").

| Study | Computation (SLURM file or command) | Figures and tables | Output |
| --- | --- | --- | --- |
| Figure 1 | `motivation.slurm`, then `motivation_pairwise_export.slurm`, the MATLAB stage (see `01_motivation/README.md`), and `motivation_pairwise_import.slurm` | `Rscript reproducibility/01_motivation/plot_motivation.R` | `motivation/motivation_intro_2x3.pdf` |
| Comparative study | `comparison_<method>.slurm` for the nine methods that run in R and Python, then `export_pairwise_comparison.slurm`, the MATLAB stage (see `02_comparison/README.md`), and `import_pairwise_comparison.slurm`. Then `lewis_runtime.slurm` for Appendix E.1. | `Rscript reproducibility/02_comparison/aggregate_and_plot.R` | `comparison/comparison_shd_vs_{p,rho,n,scale}.pdf`, `comparison/tables/*.tex`, `02_comparison/results/comparison_cell_summary.csv` |
| Calibration | `calibration_d_sweep.slurm` and `calibration_alpha_sweep.slurm` | `Rscript reproducibility/03_calibration/plot_calibration.R` | `sensitivity/sensitivity_regret.pdf`, `sensitivity/sensitivity_shd_vs_{d,d_pure,alpha,alpha_pure}.pdf`, the tables of `03_calibration/results/` |
| Uncertainty quantification | `Rscript reproducibility/04_uq/choose_seed_for_demo.R`, then `Rscript reproducibility/04_uq/uq_decision_rule_demo.R` (a few minutes on one core, or `uq_choose_seed.slurm` and `uq_demo.slurm`) | `Rscript reproducibility/04_uq/replot_uq_figures.R` | `uq/hasse_partition_compact.pdf`, `uq/posterior_marginal_annotated.pdf`, the tables of `04_uq/results/` |
| Real data | On a laptop: `cd reproducibility/05_real_data`, then `Rscript download_tcga.R` and `Rscript run_reported_analyses.R` | run by `run_reported_analyses.R`, or `Rscript plot_real_data_figures.R` | `real_data/TCGA/real_data_uq.pdf`, `real_data/TCGA/posterior_marginal_annotated_*.pdf`, the tables of `05_real_data/results/tables/` |
| Appendix F | `laplace_vs_exact.slurm` | `Rscript reproducibility/06_laplace_vs_exact/make_table_lve.R` | `06_laplace_vs_exact/results/exact_vs_laplace_summary.csv`, `exact_vs_laplace_runtime_ratio.csv`, `laplace_vs_exact/table_lve_body.tex` (the body of Table 1) |
| All | | `Rscript reproducibility/paper_numbers.R` | the numbers quoted in the text, printed |

The correlations "up to 0.999" of Section 7.2 are printed by
`Rscript reproducibility/02_comparison/variance_by_scale.R`. It needs no
data and takes seconds.

### All steps without SLURM

These are all the steps, in order, as plain commands from the root of the
repository. Set `MC_CORES`, `RETICULATE_PYTHON` and `PWLING_PATH` first
(see "Software"). The simulations take thousands of CPU hours (see
"Hardware and run times").

```sh
# Figure 1 (01_motivation)
Rscript reproducibility/01_motivation/motivation_experiment.R
Rscript reproducibility/01_motivation/pairwise_lingam/export_motivation_for_matlab.R
(cd reproducibility/01_motivation/pairwise_lingam && \
  matlab -batch "addpath('../../common/pairwise_lingam'), run_pwling_batch")
Rscript reproducibility/01_motivation/pairwise_lingam/import_motivation_pairwise.R
Rscript reproducibility/01_motivation/plot_motivation.R

# Comparative study (02_comparison)
for m in lewis_complexity lewis_uniform pc ges lingam direct_lingam rcd notears tl_zhao; do
  Rscript reproducibility/02_comparison/run_one_method.R "$m"
done
Rscript reproducibility/02_comparison/pairwise_lingam/export_datasets_for_matlab.R
(cd reproducibility/02_comparison/pairwise_lingam && \
  matlab -batch "addpath('../../common/pairwise_lingam'), run_pwling_batch")
Rscript reproducibility/02_comparison/pairwise_lingam/import_pairwise_results.R
Rscript reproducibility/02_comparison/measure_lewis_runtime.R
Rscript reproducibility/02_comparison/aggregate_and_plot.R
Rscript reproducibility/02_comparison/variance_by_scale.R

# Calibration (03_calibration)
Rscript reproducibility/03_calibration/calibration_sweep.R d
Rscript reproducibility/03_calibration/calibration_sweep.R alpha
Rscript reproducibility/03_calibration/plot_calibration.R

# Uncertainty quantification (04_uq)
Rscript reproducibility/04_uq/choose_seed_for_demo.R
Rscript reproducibility/04_uq/uq_decision_rule_demo.R
Rscript reproducibility/04_uq/replot_uq_figures.R

# Real data (05_real_data), run in its folder
(cd reproducibility/05_real_data && Rscript download_tcga.R && \
  Rscript run_reported_analyses.R)

# Appendix F (06_laplace_vs_exact)
Rscript reproducibility/06_laplace_vs_exact/laplace_vs_exact.R
Rscript reproducibility/06_laplace_vs_exact/make_table_lve.R

# The numbers quoted in the text
Rscript reproducibility/paper_numbers.R
```

## Where the outputs go

Results are written to `<study>/results/`. The figures (a vector PDF and a
300-dpi PNG) and the LaTeX tables are written to

    <figure folder>/<subfolder>/

where the figure folder is the value of the environment variable
`LEWIS_FIG_DIR`, or `reproducibility/figures/` if the variable is not set.
The folder is created when needed. The subfolders (`motivation/`,
`comparison/`, `comparison/tables/`, `sensitivity/`, `uq/`,
`real_data/TCGA/`) are those that the LaTeX source of the paper reads. The
subfolder `laplace_vs_exact/` holds the body of Table 1, which the paper
types in. To write the figures straight into the figure folder of the
paper, set for example

```sh
export LEWIS_FIG_DIR=/path/to/paper/figures
```

## Hardware and run times

The simulations ran on a cluster node with 100 cores and 128 GB of memory,
one dataset per core. The times of LEWIS in the stored per-dataset results
of Figure 1, the comparative study and the calibration study were measured
with a plain R implementation of the same algorithm. The package uses
compiled code and is faster. The CPU hours below count the time of all
datasets.

| Study | Run time |
| --- | --- |
| Figure 1 | about 80 CPU hours (1 to 2 hours on 100 cores). The MATLAB stage takes minutes. |
| Comparative study | about 400 CPU hours for LEWIS (complexity), 825 for LEWIS (uniform), 300 for RCD, 120 for NOTEARS, 80 for TL, 40 for DirectLiNGAM, and less than 3 for each other method. These times do not include the SID, which adds hours at p = 50 and p = 100. The run-time job of Appendix E.1 (1200 runs of the package) takes about 70 CPU hours. |
| Calibration | about 8,500 CPU hours for each of the two sweeps (about 4 days on 100 cores). |
| Uncertainty quantification | a few minutes on one core. |
| Real data | on a laptop (Apple M2 Pro, 16 GB of memory), one chain of 10^7 iterations takes about 13 minutes (prostate) or 10 minutes (liver). Each analysis runs 15 chains, about 3.3 and 2.5 hours. The optional screen of the 18 gene sets takes several days. |
| Appendix F | about 55 CPU hours (about 1 hour on 100 cores). |
| Figures, tables and `paper_numbers.R` | minutes on a laptop. |

## Checking the numbers of the paper

Before any rerun, save the output of `paper_numbers.R` with the Zenodo
archive in place:

```sh
Rscript reproducibility/paper_numbers.R > paper_numbers_before.txt
```

After the rerun, run it again into `paper_numbers_after.txt` and compare
the two files with `diff`. The datasets and the chains of LEWIS are
seeded, so the SHD and SID values are expected to be the same. The
computing times depend on the machine.
