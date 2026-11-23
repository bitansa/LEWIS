# =============================================================================
# real_data_config.R
#
# Every setting of the real-data analysis (Section 7.4 and Appendix I).
# Every script of this folder except download_tcga.R sources this file
# after it has defined .here, the folder of the script. The only settings
# read from the environment are
#
#     RD_NETWORK   the gene set to prepare or analyze (default PRAD_ar8)
#     RD_N_ITER    the number of iterations of each chain (default: the
#                  schedule of lewis_schedule(), 10^4 for p <= 20)
#     RD_OUT_DIR   the output folder (default results/ in this folder). The
#                  screen sets it to screen/.
#
# DATA. Bulk RNA sequencing of primary tumors from The Cancer Genome Atlas
# (TCGA), as processed by the Genomic Data Commons (GDC data release 40,
# STAR pipeline, GENCODE v36) and distributed by the UCSC Xena platform.
# The observations are patients and the variables are genes. The files are
# listed in data/raw_files.csv and downloaded by download_tcga.R into
# data/raw/.
#
# The data are public. They are not in the repository. download_tcga.R
# downloads them.
#
# DESIGN. One cohort and one gene set (a "network") at a time. Each gene
# set is a documented transcriptional program of the tumor type: its main
# regulators and some of their best-documented targets. Twelve gene sets
# were fixed on biological grounds before the screen. Six more were added
# during the screen (see "Gene sets" below).
# =============================================================================

stopifnot(exists(".here"))

# ---- Shared settings: LEWIS_ALPHA, LEWIS_D and LEWIS_ALPHA_SW ---------------
source(file.path(.here, "..", "common", "settings.R"))

# ---- Input files ------------------------------------------------------------
RD_RAW_DIR       <- file.path(.here, "data", "raw")
RD_PROBEMAP_FILE <- "gencode.v36.annotation.gtf.gene.probemap"
RD_PURITY_FILE   <- "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
RD_SUBTYPE_FILE  <- "TCGASubtype.20170308.tsv.gz"
RD_QUALITY_FILE  <- "merged_sample_quality_annotations.tsv"
rd_counts_file   <- function(cohort) sprintf("TCGA-%s.star_counts.tsv.gz", cohort)
rd_clinical_file <- function(cohort) sprintf("TCGA-%s.clinical.tsv.gz", cohort)

# ---- Gene sets ----------------------------------------------------------------
# Each network is one cohort with one gene set. The genes are listed in the
# order in which the nodes are numbered in every output. Each gene has a
# role, and the roles are in regulatory order: regulators before targets.
# An edge of the point estimate from a regulator to a target agrees with
# the known regulation. An edge from a target to a regulator goes against
# it. An edge between two genes of the same role is not tested. The symbols
# are those of GENCODE v36 (for example ACP3, formerly ACPP).
#
# There are 18 gene sets, six with 14 genes and twelve with 8 genes.
#
# Twelve sets were fixed on biological grounds before the screen. These
# are six programs with 14 genes,
#   BRCA_er    breast cancer, estrogen-receptor program
#   PRAD_ar    prostate cancer, androgen-receptor program
#   KIRC_hif   kidney clear-cell carcinoma, hypoxia program (HIF)
#   LUAD_nrf2  lung adenocarcinoma, oxidative-stress program (NRF2). SRXN1
#              is not used, because it is at zero in 27 percent of the
#              samples.
#   LIHC_hep   liver cancer, hepatocyte program
#   STAD_ifn   stomach cancer, interferon-gamma response
# and a version of each with 8 genes, the main regulators of the program
# and their best-documented direct targets:
#   BRCA_er8, PRAD_ar8, KIRC_hif8, LUAD_nrf28, LIHC_hep8, STAD_ifn8.
#
# Six sets with 8 genes were added during the screen, so that several
# genes of a set would be likely to have Gaussian errors:
#   BRCA_er8m (the estrogen-receptor regulators with five other targets),
#   PRAD_chol8 and LIHC_chol8 (the SREBP-2 cholesterol program),
#   PRAD_ifn8 (the interferon-gamma response in prostate cancer),
#   LIHC_ppara8 (the PPAR-alpha fatty-acid oxidation program),
#   LIHC_myc8 (the MYC program).
# Some of them were chosen after a look at the normality of the data: a
# normality test of the covariate-adjusted expression, or the genes
# declared Gaussian in an earlier run of the screen (for BRCA_er8m, the run
# of BRCA_er). LIHC_chol8 is one of the two reported sets.
rd_set <- function(regulators, targets) {
  data.frame(symbol = c(regulators, targets),
             role = c(rep("regulator", length(regulators)),
                      rep("target", length(targets))),
             stringsAsFactors = FALSE)
}
RD_NETWORKS <- list(
  BRCA_er = list(cohort = "BRCA",
    label = "Breast cancer, estrogen-receptor program",
    genes = rd_set(c("ESR1", "FOXA1", "GATA3"),
                   c("PGR", "TFF1", "TFF3", "XBP1", "GREB1", "MYB", "CCND1",
                     "AGR2", "CA12", "PDZK1", "SLC39A6"))),
  PRAD_ar = list(cohort = "PRAD",
    label = "Prostate cancer, androgen-receptor program",
    genes = rd_set(c("AR", "FOXA1", "HOXB13"),
                   c("NKX3-1", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "SLC45A3",
                     "STEAP2", "PMEPA1", "ACP3", "FOLH1", "ZBTB16"))),
  KIRC_hif = list(cohort = "KIRC",
    label = "Kidney clear-cell carcinoma, hypoxia program",
    genes = rd_set(c("HIF1A", "EPAS1", "ARNT"),
                   c("VEGFA", "CA9", "SLC2A1", "EGLN3", "NDUFA4L2", "PDK1",
                     "ANGPTL4", "LOX", "ADM", "BHLHE40", "PLIN2"))),
  LUAD_nrf2 = list(cohort = "LUAD",
    label = "Lung adenocarcinoma, NRF2 program",
    genes = rd_set(c("NFE2L2", "KEAP1", "MAFG"),
                   c("NQO1", "GCLM", "GCLC", "TXNRD1", "AKR1C1", "AKR1C2",
                     "AKR1C3", "G6PD", "HMOX1", "ME1", "CYP4F11"))),
  LIHC_hep = list(cohort = "LIHC",
    label = "Liver cancer, hepatocyte program",
    genes = rd_set(c("HNF4A", "HNF1A", "FOXA2", "CEBPA", "ONECUT1"),
                   c("ALB", "APOA1", "APOB", "TTR", "TF", "SERPINA1", "FGB",
                     "CYP3A4", "APOC3"))),
  STAD_ifn = list(cohort = "STAD",
    label = "Stomach cancer, interferon-gamma response",
    genes = rd_set(c("STAT1", "IRF1", "IRF9", "STAT2", "CIITA"),
                   c("CXCL9", "CXCL10", "GBP1", "GBP4", "IDO1", "TAP1", "PSMB9",
                     "HLA-DRA", "IFI27"))),
  BRCA_er8 = list(cohort = "BRCA",
    label = "Breast cancer, estrogen-receptor program (8 genes)",
    genes = rd_set(c("ESR1", "FOXA1", "GATA3"),
                   c("PGR", "GREB1", "TFF1", "MYB", "XBP1"))),
  LUAD_nrf28 = list(cohort = "LUAD",
    label = "Lung adenocarcinoma, NRF2 program (8 genes)",
    genes = rd_set(c("NFE2L2", "KEAP1"),
                   c("NQO1", "GCLM", "GCLC", "TXNRD1", "HMOX1", "AKR1C1"))),
  STAD_ifn8 = list(cohort = "STAD",
    label = "Stomach cancer, interferon-gamma response (8 genes)",
    genes = rd_set(c("STAT1", "IRF1"),
                   c("CXCL9", "CXCL10", "GBP1", "IDO1", "TAP1", "PSMB9"))),
  PRAD_ar8 = list(cohort = "PRAD",
    label = "Prostate cancer, androgen-receptor program (8 genes)",
    genes = rd_set(c("AR", "FOXA1", "HOXB13"),
                   c("KLK3", "KLK2", "TMPRSS2", "NKX3-1", "FKBP5"))),
  KIRC_hif8 = list(cohort = "KIRC",
    label = "Kidney clear-cell carcinoma, hypoxia program (8 genes)",
    genes = rd_set(c("EPAS1", "HIF1A", "ARNT"),
                   c("CA9", "VEGFA", "SLC2A1", "EGLN3", "NDUFA4L2"))),
  LIHC_hep8 = list(cohort = "LIHC",
    label = "Liver cancer, hepatocyte program (8 genes)",
    genes = rd_set(c("HNF4A", "HNF1A", "FOXA2"),
                   c("ALB", "APOA1", "TTR", "TF", "SERPINA1"))),
  BRCA_er8m = list(cohort = "BRCA",
    label = "Breast cancer, estrogen-receptor program (8 genes, second set)",
    genes = rd_set(c("ESR1", "FOXA1", "GATA3"),
                   c("PGR", "GREB1", "TFF3", "PDZK1", "SLC39A6"))),
  PRAD_chol8 = list(cohort = "PRAD",
    label = "Prostate cancer, SREBP-2 cholesterol program (8 genes)",
    genes = rd_set("SREBF2",
                   c("HMGCR", "HMGCS1", "LDLR", "SQLE", "INSIG1", "FDFT1", "DHCR7"))),
  PRAD_ifn8 = list(cohort = "PRAD",
    label = "Prostate cancer, interferon-gamma response (8 genes)",
    genes = rd_set(c("STAT1", "IRF1"),
                   c("CXCL9", "CXCL10", "GBP1", "IDO1", "TAP1", "PSMB9"))),
  LIHC_chol8 = list(cohort = "LIHC",
    label = "Liver cancer, SREBP-2 cholesterol program (8 genes)",
    genes = rd_set("SREBF2",
                   c("HMGCR", "HMGCS1", "LDLR", "SQLE", "INSIG1", "FDFT1", "PCSK9"))),
  LIHC_ppara8 = list(cohort = "LIHC",
    label = "Liver cancer, PPAR-alpha fatty-acid oxidation program (8 genes)",
    genes = rd_set("PPARA",
                   c("CPT1A", "ACOX1", "HMGCS2", "ACADM", "PDK4", "CYP4A11", "ANGPTL4"))),
  LIHC_myc8 = list(cohort = "LIHC",
    label = "Liver cancer, MYC program (8 genes)",
    genes = rd_set("MYC",
                   c("NPM1", "NCL", "ODC1", "LDHA", "CDK4", "NOP56", "SRM")))
)

# ---- The network of this run ----------------------------------------------------
RD_NETWORK <- Sys.getenv("RD_NETWORK", "PRAD_ar8")
if (!RD_NETWORK %in% names(RD_NETWORKS)) {
  stop("Unknown network RD_NETWORK = ", RD_NETWORK, ". Choose one of: ",
       paste(names(RD_NETWORKS), collapse = ", "))
}
RD_COHORT <- RD_NETWORKS[[RD_NETWORK]]$cohort
RD_GENES <- RD_NETWORKS[[RD_NETWORK]]$genes

# ---- Samples -------------------------------------------------------------------
# Primary tumors only (TCGA sample type 01). One sample per patient (the
# first in barcode order). Samples preserved in formalin (FFPE) and samples
# flagged "Do_not_use" for RNA sequencing in the Pan-Cancer Atlas quality
# annotations are removed. A sex with fewer than RD_MIN_SEX_FRAC of the
# samples is removed (in practice the male breast cancers). Otherwise sex
# is a covariate when both sexes occur.
RD_SAMPLE_TYPE <- "01"
RD_MIN_SEX_FRAC <- 0.05

# ---- Sample quality control and normalization ----------------------------------
# The library size of a sample is its total count over all genes. A sample
# is kept if its library size is at least RD_LIB_MIN_FRAC of the median and
# at most RD_LIB_MAD_K scaled median absolute deviations above the median.
# The counts are then scaled to the median library size and transformed as
# log(1 + x). A gene with more than RD_MAX_ZERO_FRAC of the samples at zero
# is dropped and recorded (no gene of the reported sets is dropped).
RD_LIB_MIN_FRAC <- 0.25
RD_LIB_MAD_K <- 3
RD_MAX_ZERO_FRAC <- 0.10

# ---- Covariates ------------------------------------------------------------------
# Each gene is replaced by its residual from a linear regression on
#   age      age at diagnosis (years)
#   sex      when both sexes remain
#   purity   the tumor purity estimated by ABSOLUTE
#   site     the tissue source site (characters 6 and 7 of the barcode).
#            Sites with fewer than RD_MIN_SITE samples are pooled.
#   subtype  the Pan-Cancer Atlas molecular subtype, when it is known for at
#            least RD_MIN_SUBTYPE_FRAC of the samples. It is known for 67
#            percent of the prostate samples and 49 percent of the liver
#            samples, so it is not used in the reported analyses.
# The samples without a value of a covariate that is used are removed.
RD_MIN_SITE <- 10L
RD_MIN_SUBTYPE_FRAC <- 0.90

# ---- LEWIS -----------------------------------------------------------------------
# As in the simulations: the complexity prior with alpha = 0.75 and
# d = 0.03, the Shapiro-Wilk level 0.05 with the Bonferroni correction, the
# start at the TL estimate, a burn-in of 20 percent and at most 8000 draws.
# With 10^7 iterations, lewis_schedule() then keeps every 1000th state after
# the burn-in, which gives M = 8000 draws (Appendix I). The data are not
# centered again by lewis() (center = FALSE), because prepare_tcga.R has
# centered and scaled every gene.
# The number of iterations is read from RD_N_ITER. If it is not set, the
# schedule of lewis_schedule() is used.
RD_N_ITER <- local({
  v <- Sys.getenv("RD_N_ITER", "")
  if (nzchar(v)) as.integer(v) else NULL
})
RD_BURN_IN <- 0.2
RD_MAX_DRAWS <- 8000L
# One minus the level of the credible ball, and the number of graphs shown
# in Figure 19.
RD_ALPHA_BALL <- 0.05
RD_TOP_K <- 15L

# ---- Seeds -----------------------------------------------------------------------
# The reported chain uses RD_SEED. The stability check repeats the whole
# procedure with the ten seeds RD_STABILITY_SEEDS. The convergence check runs
# four additional chains with the seeds RD_CONV_SEEDS. The first two start
# at the TL estimate and the last two at the empty graph.
RD_SEED <- 1L
RD_STABILITY_SEEDS <- 2:11
RD_CONV_SEEDS <- 101:104

# ---- Convergence checks (Appendix I) -------------------------------------------
# A gene set passes when the four conditions hold:
#   1. the potential scale reduction factor (PSRF) of the score of the four
#      additional chains is below RD_CONV_MAX_PSRF.
#   2. the highest scores of the chains agree: the highest scores reached
#      by the four additional chains differ by at most RD_CONV_MAX_GAP = 0.5,
#      in the units of the score s_n (a log posterior).
#   3. the posterior probability of every adjacency and every orientation
#      differs between the chains by less than RD_CONV_MAX_PDIFF.
#   4. all ten repetitions with the seeds RD_STABILITY_SEEDS return the same
#      point estimate D_hat as the reported chain.
RD_CONV_MAX_PSRF <- 1.05
RD_CONV_MAX_GAP <- 0.5
RD_CONV_MAX_PDIFF <- 0.1

# ---- The screen ------------------------------------------------------------------
# screen_networks.R analyzes a gene set at the chain lengths
# RD_SCREEN_LENGTHS in turn (NA is the schedule of lewis_schedule(), then
# 10^5, 10^6 and 10^7 iterations) and stops at the first length at which it
# passes. The sets with 14 genes are run with up to 10^6 iterations, the
# sets with 8 genes with up to 10^7 iterations.
RD_SCREEN_LENGTHS <- c(NA, 100000L, 1000000L, 10000000L)
rd_screen_max_length <- function(network) {
  if (nrow(RD_NETWORKS[[network]]$genes) >= 14L) 1000000L else 10000000L
}

# ---- The reported analyses ---------------------------------------------------------
# Three gene sets pass the checks: PRAD_ar8 (at 10^7 iterations), LIHC_hep8
# (10^6) and LIHC_chol8 (10^7). Two are reported, both with 10^7 iterations:
# one with a concentrated posterior (prostate) and one with a spread-out
# posterior (liver).
RD_REPORTED <- data.frame(network = c("PRAD_ar8", "LIHC_chol8"),
                          n_iter = c(10000000L, 10000000L),
                          stringsAsFactors = FALSE)

# ---- Output names and folders ----------------------------------------------------
# The prepared data do not depend on the chain length. Their files carry the
# name of the network. The outputs of the analysis carry RD_TAG, which adds
# the chain length in thousands of iterations (for example PRAD_ar8_it10000k).
rd_length_label <- function(n_iter) {
  if (is.null(n_iter) || is.na(n_iter)) "sched" else paste0(n_iter %/% 1000L, "k")
}
rd_tag <- function(network, n_iter) sprintf("%s_it%s", network, rd_length_label(n_iter))
RD_TAG <- rd_tag(RD_NETWORK, RD_N_ITER)

# The output folder of this run and its folder of tables.
RD_OUT_DIR <- Sys.getenv("RD_OUT_DIR", file.path(.here, "results"))
RD_TAB_DIR <- file.path(RD_OUT_DIR, "tables")
dir.create(RD_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
