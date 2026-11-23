# 05_real_data: the TCGA analysis

This folder reproduces the real-data analysis of the paper:

- Section 7.4. Figure 12 shows the two point estimates. It is drawn in
  LaTeX (TikZ) from `results/tables/point_estimate_*.csv` and
  `results/tables/node_posterior_*.csv`.
- Appendix I: the data and the preprocessing, and the choice of the gene
  sets (I.1), the chain length and the convergence checks (I.2), Table 2
  (I.3), and Figures 18 and 19 (I.4).
- Table S15 in Section S3 of the online appendix: the posterior
  probabilities of the edges of the two point estimates.

The two reported analyses are `PRAD_ar8` (prostate cancer,
androgen-receptor program, n = 462) and `LIHC_chol8` (liver cancer, SREBP-2
cholesterol program, n = 348). Both have p = 8 genes.

LEWIS (complexity) is run as in the simulations. It uses the complexity
prior with alpha = 0.75 and d = 0.03, the Shapiro-Wilk level 0.05 with the
Bonferroni correction, and the start at the TL estimate. Every call is
`lewis(..., center = FALSE)`, because `prepare_tcga.R` has already centered
and scaled the genes. Each chain has 10^7 iterations with a burn-in of 20
percent. Every 1000th state after the burn-in is kept, which gives 8000
draws. The reported chain uses the seed 1.

## Data

The TCGA data are public and are not in this repository or the Zenodo
archive. `download_tcga.R` downloads them, and `prepare_tcga.R` (called by
`run_reported_analyses.R`) does the preprocessing.

`download_tcga.R` saves the files in `data/raw/` and checks each file
against the md5 sum in `data/raw_files.csv`. These are the md5 sums of the
files used for the paper. The file `data/raw_files.csv` also gives the size
and the download address of each file. The addresses are on three public
servers: the UCSC Xena GDC hub (`https://gdc.xenahubs.net/download/...`),
the UCSC Xena Pan-Cancer Atlas hub
(`https://pancanatlas.xenahubs.net/download/...`) and the GDC data portal
(`https://api.gdc.cancer.gov/data/...`). No login is needed.

| File | Content | Source | Needed for |
| --- | --- | --- | --- |
| `TCGA-<cohort>.star_counts.tsv.gz` | gene counts, log2(count + 1), STAR pipeline, GENCODE v36 | UCSC Xena GDC hub, GDC data release 40 | PRAD and LIHC: the reported analyses. BRCA, KIRC, LUAD and STAD: the screen only |
| `TCGA-<cohort>.clinical.tsv.gz` | sample type, preservation, age, sex | UCSC Xena GDC hub, GDC data release 41 | as above |
| `gencode.v36.annotation.gtf.gene.probemap` | gene identifier to gene symbol | UCSC Xena GDC hub | all analyses |
| `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | tumor purity by ABSOLUTE | Pan-Cancer Atlas supplemental data (GDC) | all analyses |
| `merged_sample_quality_annotations.tsv` | sample quality flags (`Do_not_use`) | Pan-Cancer Atlas supplemental data (GDC) | all analyses |
| `TCGASubtype.20170308.tsv.gz` | molecular subtypes | UCSC Xena Pan-Cancer Atlas hub | all analyses (the share of samples with a known subtype) |

The files of the two reported analyses take about 115 MB. All files,
including those that only the screen needs, take about 450 MB.

The UCSC Xena hubs may replace a file by a newer data release. The md5
check then fails and the script stops. `TCGA_ALLOW_MD5_MISMATCH=1` keeps
such a file, but the results may then differ from the paper.

## Preprocessing (`prepare_tcga.R`)

1. Samples: primary tumors (sample type 01), one sample per patient. The
   samples preserved in formalin (FFPE) and the samples flagged
   `Do_not_use` for RNA sequencing are removed.
2. Quality control: a sample is removed if its library size is below one
   quarter of the median, or more than 3 scaled median absolute deviations
   above the median.
3. Normalization: the counts are scaled to the median library size and
   transformed as log(1 + x).
4. Genes: the genes of the set, mapped from their symbols to the GENCODE
   v36 identifiers.
5. Covariates: the samples without a value of a covariate that is used are
   removed. Each gene is replaced by its residual from a linear regression
   on age, sex, tumor purity and tissue source site. A sex with fewer than
   5 percent of the samples is removed. This happens only in the breast
   cancer cohort of the screen. Sex is used only when both sexes occur.
   Of the two reported analyses, this is the liver analysis. The sites with
   fewer than 10 samples are pooled. The molecular subtype would be used if
   it were known for at least 90 percent of the samples. It is known for 67
   percent of the prostate samples and 49 percent of the liver samples, so
   it is not used.
6. Each gene is centered and scaled to unit standard deviation.

Every count is recorded in `results/tables/qc_summary_<network>.csv`.

## The gene sets and the screen (`screen_networks.R`)

The 18 gene sets are defined in `real_data_config.R`. Six sets have 14
genes and twelve sets have 8 genes. Each set has the main regulators of a
known program and some of their best-documented targets. The sets were
chosen in two steps.

1. Twelve sets were fixed on biological grounds before the screen. These
   are six programs with 14 genes (`BRCA_er`, `PRAD_ar`, `KIRC_hif`,
   `LUAD_nrf2`, `LIHC_hep` and `STAD_ifn`) and a version of each with 8
   genes (`BRCA_er8`, `PRAD_ar8`, `KIRC_hif8`, `LUAD_nrf28`, `LIHC_hep8`
   and `STAD_ifn8`).
2. Six sets with 8 genes were added during the screen, so that several
   genes of a set would be likely to have Gaussian errors. These are
   `BRCA_er8m`, `PRAD_chol8`, `PRAD_ifn8`, `LIHC_chol8`, `LIHC_ppara8` and
   `LIHC_myc8`. Some of them were chosen after a look at the normality of
   the data. This look was a normality test of the covariate-adjusted
   expression, or the genes declared Gaussian in an earlier run of the
   screen (for `BRCA_er8m`, the run of `BRCA_er`). The reported liver set
   `LIHC_chol8` is one of these six sets.

The screen analyzes each set with `real_data_lewis.R` at increasing chain
lengths: the default of `lewis_schedule()` (10^4 iterations for p <= 20),
then 10^5, 10^6 and 10^7 iterations. The sets with 14 genes stop at 10^6
iterations. A set stops at the first length at which it passes four
checks:

1. The potential scale reduction factor of the score of the four
   additional chains is below 1.05.
2. The highest scores of the four additional chains agree. They differ by
   at most 0.5 in the units of the score (`RD_CONV_MAX_GAP`).
3. The posterior probability of every adjacency and every orientation
   differs between the four chains by less than 0.1.
4. All ten repetitions with other seeds return the same point estimate.

The summary of the screen is `results/tables/screen_summary.csv`. It
records the following runs.

- The six sets with 14 genes were run with up to 10^6 iterations. None of
  them passed.
- Seven sets with 8 genes were run with up to 10^7 iterations. Two of them
  passed at 10^7 iterations: `PRAD_ar8` and `LIHC_chol8`.
- `LIHC_hep8` passed at 10^6 iterations.
- The other four sets with 8 genes (`PRAD_chol8`, `PRAD_ifn8`,
  `LIHC_ppara8` and `LIHC_myc8`) were stopped after 10^6 iterations, once
  three sets had passed. Their runs with 10^7 iterations were not
  completed.

So three sets passed: `PRAD_ar8`, `LIHC_chol8` and `LIHC_hep8`. The paper
reports `PRAD_ar8` and `LIHC_chol8`. A full rerun of the screen also runs
the four stopped sets with 10^7 iterations and adds their rows to the
summary.

## Software and run order

The scripts need R, the package LEWIS of this repository, and the R
packages `ggplot2` and `scales` (figures), `pcalg` (only if a point
estimate has undirected edges), and `data.table` with `R.utils`
(optional, for faster reading of the large files). They use the files
`../common/settings.R`, `../common/paths.R` and `../common/plot_theme.R`.
The run scripts start the other scripts with `Rscript` and pass settings
through environment variables, so they need macOS or Linux. No SLURM file
is needed. The analysis runs on a laptop.

Run the scripts with `Rscript`. Each script finds its own folder, so it can
be started from any folder. The commands below are run in the folder
`reproducibility/05_real_data`. From the root of the repository, go there
first:

```sh
cd reproducibility/05_real_data
Rscript download_tcga.R              # the files of the reported analyses (about 115 MB)
Rscript run_reported_analyses.R      # prepare, analyze and plot the two reported analyses
```

`run_reported_analyses.R` runs, for `PRAD_ar8` and then `LIHC_chol8`,

1. `prepare_tcga.R`: the patient-by-gene matrix.
2. `real_data_lewis.R`: the reported chain, the ten stability runs and the
   four convergence chains, each with 10^7 iterations.

and then `plot_real_data_figures.R`. The output of each step is written to
a log file in `results/logs/`. The new tables replace the tables stored in
`results/tables/`, so `git diff results/tables` shows any difference. Run
times, effective sample sizes and the last digits of the scores can differ.

The two analyses can also run in parallel in two terminals:

```sh
REPORTED_ONLY=PRAD_ar8 Rscript run_reported_analyses.R
REPORTED_ONLY=LIHC_chol8 Rscript run_reported_analyses.R
```

Each run redraws the figures from the tables present in `results/tables/`.
The figures of the run that finishes last use the new tables of both
analyses.

Each step can also be run alone:

```sh
RD_NETWORK=PRAD_ar8 Rscript prepare_tcga.R
RD_NETWORK=PRAD_ar8 RD_N_ITER=10000000 Rscript real_data_lewis.R
Rscript plot_real_data_figures.R     # the figures from the stored tables, in seconds
```

On a Mac, `caffeinate -i Rscript run_reported_analyses.R` keeps the
computer awake during the run.

The screen is optional. It needs all files and takes several days on a
laptop:

```sh
Rscript download_tcga.R all                                 # all files (about 450 MB)
SCREEN_RUN=all Rscript screen_networks.R                    # the screen of the 18 gene sets
SCREEN_RUN=PRAD_ar8,LIHC_chol8 Rscript screen_networks.R    # some gene sets only
Rscript screen_networks.R                                   # rebuild the summary from screen/
```

The analyses of the screen are written to `screen/` (or to the folder in
`RD_SCREEN_DIR`). The summary in `results/tables/screen_summary.csv` is
rebuilt from the gene sets found in `screen/`. If `screen/` holds no
analysis, the stored summary is kept. If `screen/` holds only some gene
sets, the new summary has rows for these sets only. It then replaces the
stored summary of the 18 sets. The command
`git checkout results/tables/screen_summary.csv` restores the stored
summary.

## Outputs

Tables, in `results/tables/`. Here `<network>` is `PRAD_ar8` or
`LIHC_chol8`, and `<tag>` is `PRAD_ar8_it10000k` or `LIHC_chol8_it10000k`.

| File | Content | Used for |
| --- | --- | --- |
| `expression_<network>.csv` | the patient-by-gene matrix | input of LEWIS |
| `gene_map_<network>.csv` | gene symbols and identifiers | Appendix I.1 |
| `qc_summary_<network>.csv` | the number of samples at each step | Appendix I.1 |
| `gene_normality_<network>.csv` | normality of each gene across patients (descriptive) | not reported, and not used by the analysis |
| `summary_<tag>.csv` | size, Gaussian genes, agreement counts, effective sample sizes | Section 7.4, Appendix I.2 |
| `point_estimate_<tag>.csv` | the edges of D_hat with their posterior probability and their order relative to the known regulation | Figure 12, Table S15 |
| `pair_posterior_<tag>.csv` | the posterior probabilities of the four states (i -> j, j -> i, undirected, absent) of every pair of genes | not reported (it extends Table S15 to all pairs) |
| `node_posterior_<tag>.csv` | Shapiro-Wilk p-values and the proportion of draws in which each gene is declared non-Gaussian | Figure 12, Table 2 |
| `posterior_rpdag_<tag>.csv` | every distinct RPDAG with its mass, SHD to D_hat and class | Section 7.4, Appendix I.4 |
| `uq_<tag>/` | the 15 most probable RPDAGs, the summary of the chain (radius, mass outside the ball, run time, effective sample size) and the gene names | Figure 19, Appendix I.2 and I.4 |
| `stability_<tag>.csv` | the point estimates of the ten other seeds | Appendix I.2 |
| `convergence_<tag>.csv` | the four additional chains | Appendix I.2 |
| `real_data_uq_by_shd.csv` | the masses of Figure 18 | Figure 18 |
| `screen_summary.csv` | the screen of the 18 gene sets | Appendix I.1 and I.2 |

The column `known_order` of `point_estimate_<tag>.csv` is "agrees",
"against" or "within tier". The paper calls the last case "same level":
an edge between two regulators or between two targets. The same counts are
in the columns `agrees`, `against` and `within_tier` of `summary_<tag>.csv`.

Figures, in `<LEWIS_FIG_DIR>/real_data/TCGA/` (by default
`reproducibility/figures/real_data/TCGA/`), each as a PDF and a PNG:

| File | Figure |
| --- | --- |
| `real_data_uq.pdf` | Figure 18 |
| `posterior_marginal_annotated_PRAD_ar8_it10000k.pdf` | Figure 19(a) |
| `posterior_marginal_annotated_LIHC_chol8_it10000k.pdf` | Figure 19(b) |

`paper_numbers.R` (in the parent folder) prints the numbers of Section 7.4
and Appendix I from these tables.

## Run time and hardware

The analysis was run in R 4.6.1 on a laptop (Apple M2 Pro, 16 GB of
memory). One chain of 10^7 iterations takes about 13 minutes for the
prostate analysis and 10 minutes for the liver analysis.
`real_data_lewis.R` runs 15 such chains (the reported chain, ten stability
runs and four convergence chains). So `run_reported_analyses.R` takes about
3.3 hours for the prostate analysis and 2.5 hours for the liver analysis.
The download takes a few minutes. The figures take seconds.
