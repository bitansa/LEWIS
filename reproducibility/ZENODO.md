# Results of each dataset (Zenodo archive)

The folders `results/` of this repository hold only the small summary files
from which the numbers, tables and figures of the paper are drawn. The
results of each dataset of the simulations are too large for the
repository. They are deposited on Zenodo:

> <https://doi.org/10.5281/zenodo.24514618>

The archive has one file per study. Each file holds paths relative to
`reproducibility/`. To use it, download the files of the deposit. First
check them with the md5 sums in `MD5SUMS.txt` of the deposit. Run the
check in the folder of the download (`md5sum -c MD5SUMS.txt` on Linux, or
`md5 -r lewis_raw_*.tar.gz` on macOS and compare with `MD5SUMS.txt`). Then
unpack the files in `reproducibility/`:

```sh
cd reproducibility
for f in /path/to/download/lewis_raw_*.tar.gz; do tar -xzf "$f"; done
```

Every file then lands where the scripts read it. `FILES_MD5.txt` of the
deposit gives the md5 sum of every unpacked file.

With the archive in place, every figure and table of the paper can be drawn
again on a laptop, without running the simulations.

## What is not in the archive

- **The TCGA data.** They are public and are not archived. The script
  `05_real_data/download_tcga.R` downloads them from their public sources
  and checks their md5 sums. The preprocessing is done by
  `05_real_data/prepare_tcga.R`. The tables of the real-data analysis are
  small and are in `05_real_data/results/tables/` of this repository.
- **The inputs of Pairwise-LiNGAM** (`pairwise_lingam/inputs/`, about
  0.5 GB for Figure 1 and 1.9 GB for the comparative study). The export
  scripts write them again, exactly, from the seeds.
- **The list of exported datasets of Figure 1**
  (`01_motivation/pairwise_lingam/manifest.csv`). The export script writes
  it again.
- **`02_comparison/results/comparison_all.csv`.** `aggregate_and_plot.R`
  writes it from the files of each method.

## Contents

Sizes are those of the unpacked files.

| Archive | Files | Size | Read by |
| --- | --- | --- | --- |
| `lewis_raw_01_motivation.tar.gz` (4.9 MB) | `01_motivation/results/motivation_results.csv`: one row per regime, sample size, dataset and method, for the 25 regimes of Figure 1 and the 10 methods (75,000 rows) | 8.9 MB | `01_motivation/plot_motivation.R`, `check_progress.R` |
| | `01_motivation/pairwise_lingam/outputs/`: 7,500 likelihood-ratio matrices of Pairwise-LiNGAM (`<regime>_n<n>_r<dataset>_LR.csv`) and `runtimes.csv` | 10.7 MB | `01_motivation/pairwise_lingam/import_motivation_pairwise.R` |
| `lewis_raw_02_comparison.tar.gz` (35.5 MB) | `02_comparison/results/results_<method>.csv`: one file per method, one row per setting and dataset (5,400 rows each, 10 methods) | 4.2 MB | `02_comparison/aggregate_and_plot.R`, `02_comparison/measure_lewis_runtime.R` (check of the SHD), `paper_numbers.R`, `check_progress.R` |
| | `02_comparison/pairwise_lingam/manifest.csv`: the list of the exported datasets | 0.5 MB | `02_comparison/pairwise_lingam/import_pairwise_results.R` |
| | `02_comparison/pairwise_lingam/outputs/`: 5,400 likelihood-ratio matrices of Pairwise-LiNGAM and `runtimes.csv` | 105.6 MB | `02_comparison/pairwise_lingam/import_pairwise_results.R` |
| `lewis_raw_03_calibration.tar.gz` (2.3 MB) | `03_calibration/results/d_sweep_results.csv` and `alpha_sweep_results.csv`: one row per setting, value and dataset (84,000 rows each) | 10.8 MB | `03_calibration/plot_calibration.R`, `check_progress.R` |
| `lewis_raw_04_uq.tar.gz` (0.4 MB) | `04_uq/results/tables/samples_per_dag.csv`: one row per retained draw of the chain on the selected dataset | 1.9 MB | `paper_numbers.R` (mass outside the credible ball) |
| | `04_uq/results/seed_sweep/seed_01/` to `seed_10/tables/`: the tables of each of the 10 datasets of the seed rule | 18.6 MB | none (for reference) |
| `lewis_raw_06_laplace_vs_exact.tar.gz` (0.2 MB) | `06_laplace_vs_exact/results/exact_vs_laplace.csv`: one row per dataset, prior and chain (5,200 rows) | 1.7 MB | `06_laplace_vs_exact/make_table_lve.R`, `check_progress.R` |

## Notes

- The column `shd_rpdag` is the SHD from the true RPDAG. In the files of
  `04_uq`, the column `distinct_rpdags` and the file name
  `top_rpdag_posterior.csv` use the same term.
- The tables of `04_uq` in the archive are those of the paper. Some of
  their columns differ from those that the scripts of `04_uq` write (see
  `04_uq/README.md`).

## How the archive was made

The archive was made by `make_zenodo_archive.sh` from a copy of this folder
in which the studies were run:

```sh
bash reproducibility/make_zenodo_archive.sh <reproducibility folder> <output folder>
```

Users of the repository do not need to run it.
