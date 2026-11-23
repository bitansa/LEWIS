# =============================================================================
# prepare_tcga.R
#
# Builds the patient-by-gene matrix that LEWIS is fitted to, for one network
# (RD_NETWORK) of real_data_config.R. The steps are those of Appendix I
# ("Data and preprocessing"):
#
#   1. samples: primary tumors, not preserved in formalin, not flagged
#      "Do_not_use", one sample per patient.
#   2. quality control: samples with an outlying library size are removed.
#   3. normalization: the counts are scaled to the median library size and
#      transformed as log(1 + x).
#   4. genes: the genes of the set, mapped from their symbols to the
#      GENCODE v36 identifiers of the data.
#   5. covariate adjustment: each gene is replaced by its residual on age,
#      sex, tumor purity, tissue source site and, if it is known for enough
#      samples, the molecular subtype.
#   6. centering and scaling: each gene to mean 0 and standard deviation 1.
# Outlying values are not trimmed, because the non-Gaussianity of the errors
# is what orients the edges.
#
# Outputs, in <RD_OUT_DIR>/tables/:
#   expression_<network>.csv      patients x genes, the matrix LEWIS is
#                                 fitted to
#   gene_map_<network>.csv        symbol, role, identifier in the data,
#                                 number of identifiers, share of samples
#                                 at zero
#   qc_summary_<network>.csv      the number of samples at each step and the
#                                 choices made
#   gene_normality_<network>.csv  Shapiro-Wilk p-value, skewness and excess
#                                 kurtosis of each gene across patients.
#                                 This is descriptive only. LEWIS classifies
#                                 the genes from residuals (Step 1 of the
#                                 post-processing).
#
# Run:  Rscript prepare_tcga.R                    (the network PRAD_ar8)
#       RD_NETWORK=LIHC_chol8 Rscript prepare_tcga.R
# =============================================================================

# ---- Folder of this script and settings -------------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
source(file.path(.here, "real_data_config.R"))

# ---- Helpers ------------------------------------------------------------------
# The quality-control counts, collected into qc_summary_<network>.csv.
qc <- list()
note <- function(key, value) {
  qc[[key]] <<- value
  invisible(value)
}

# The path of a raw file. Stop if it is missing.
raw <- function(f) {
  path <- file.path(RD_RAW_DIR, f)
  if (!file.exists(path)) {
    stop("The input file ", path, " is missing. Run download_tcga.R first.")
  }
  path
}

# Read a tab-separated file. data.table::fread is used when it is installed,
# because it is much faster on the large expression files. The quality
# annotations have quoted text fields, some over several lines, so that file
# is read with quoting. The other files have no quotes.
read_tsv <- function(path, quote = "") {
  if (requireNamespace("data.table", quietly = TRUE) &&
      (!grepl("\\.gz$", path) || requireNamespace("R.utils", quietly = TRUE))) {
    return(as.data.frame(data.table::fread(path, sep = "\t", quote = quote,
                                           data.table = FALSE, check.names = FALSE,
                                           na.strings = c("", "NA"))))
  }
  utils::read.delim(if (grepl("\\.gz$", path)) gzfile(path) else path,
                    check.names = FALSE, stringsAsFactors = FALSE, quote = quote,
                    na.strings = c("", "NA"))
}

# Print and record the network of this run.
cat(sprintf("Network %s (%d genes): %s\n", RD_NETWORK, nrow(RD_GENES),
            RD_NETWORKS[[RD_NETWORK]]$label))
note("network", RD_NETWORK)
note("cohort", RD_COHORT)
note("expression", "counts")

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
# The expression file has one row per gene and one column per sample. Its
# values are log2(count + 1).
t0 <- Sys.time()
ex <- read_tsv(raw(rd_counts_file(RD_COHORT)))
ids <- ex[[1L]]
ex <- as.matrix(ex[, -1L, drop = FALSE])
rownames(ex) <- ids
cat(sprintf("Counts: %d genes x %d samples, read in %.0f s\n", nrow(ex), ncol(ex),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
note("n_genes_in_data", nrow(ex))
note("n_samples_in_data", ncol(ex))

# The clinical table, the quality annotations, the purity, the subtypes and
# the map from identifiers to gene symbols.
clin <- read_tsv(raw(rd_clinical_file(RD_COHORT)))
clin <- clin[!duplicated(clin$sample), , drop = FALSE]
rownames(clin) <- clin$sample
qual <- read_tsv(raw(RD_QUALITY_FILE), quote = "\"")
if (nrow(qual) < 79000L) {
  stop("The quality annotations were read incompletely (", nrow(qual), " rows).")
}
pur <- read_tsv(raw(RD_PURITY_FILE))
subt <- read_tsv(raw(RD_SUBTYPE_FILE))
pm <- read_tsv(raw(RD_PROBEMAP_FILE))

# -----------------------------------------------------------------------------
# 1. Samples
# -----------------------------------------------------------------------------
# A sample barcode such as TCGA-A1-A0SB-01A has the patient in characters 1
# to 12, the sample type in characters 14 and 15 (01 = primary tumor) and
# the vial in character 16.
s_all <- colnames(ex)
s1 <- s_all[substr(s_all, 14L, 15L) == RD_SAMPLE_TYPE]
note("n_primary_tumor", length(s1))

# Every primary tumor must be in the clinical table.
not_in_clin <- setdiff(s1, rownames(clin))
if (length(not_in_clin) > 0L) {
  stop(length(not_in_clin), " samples are not in the clinical table, for example ",
       paste(head(not_in_clin, 3L), collapse = ", "), ".")
}

# Remove the samples preserved in formalin (FFPE).
pres <- clin[s1, "preservation_method.samples"]
s2 <- s1[is.na(pres) | pres != "FFPE"]
note("removed_ffpe", length(s1) - length(s2))

# Remove the samples flagged "Do_not_use" for RNA sequencing.
bad <- qual$aliquot_barcode[grepl("RNASeq", qual$platform) &
                            qual$Do_not_use %in% c(TRUE, "True", "TRUE")]
s3 <- s2[!s2 %in% substr(bad, 1L, 16L)]
note("removed_do_not_use", length(s2) - length(s3))

# Keep the first sample of each patient in barcode order.
s3 <- sort(s3)
s4 <- s3[!duplicated(substr(s3, 1L, 12L))]
note("removed_second_sample_of_patient", length(s3) - length(s4))
cat(sprintf(paste0("Samples: %d primary tumors. Removed %d FFPE, %d flagged and %d ",
                   "second samples. %d kept.\n"),
            length(s1), length(s1) - length(s2), length(s2) - length(s3),
            length(s3) - length(s4), length(s4)))
ex <- ex[, s4, drop = FALSE]

# -----------------------------------------------------------------------------
# 2. Sample quality control
# -----------------------------------------------------------------------------
# Back from log2(count + 1) to counts, then the library size of each sample
# relative to the median.
counts <- 2^ex - 1
lib <- colSums(counts)
rel <- lib / median(lib)
keep <- rel >= RD_LIB_MIN_FRAC & rel <= 1 + RD_LIB_MAD_K * mad(rel)
note("lib_min_frac", RD_LIB_MIN_FRAC)
note("lib_mad_k", RD_LIB_MAD_K)
note("removed_library_size", sum(!keep))
cat(sprintf("Library size: %d outliers removed, %d kept.\n", sum(!keep), sum(keep)))
counts <- counts[, keep, drop = FALSE]

# -----------------------------------------------------------------------------
# 3. Normalization
# -----------------------------------------------------------------------------
sc <- colSums(counts) / median(colSums(counts))
logx <- log1p(sweep(counts, 2L, sc, "/"))
note("normalization", "counts / (library size / median library size), then log(1 + x)")

# -----------------------------------------------------------------------------
# 4. Genes
# -----------------------------------------------------------------------------
# The symbols are mapped to the GENCODE v36 identifiers of the data through
# the probe map. A symbol with several identifiers takes the one with the
# largest mean count. A symbol without an identifier is dropped and
# recorded.
map <- RD_GENES
map$id_in_data <- vapply(map$symbol, function(g) {
  cand <- intersect(pm$id[pm$gene == g], rownames(counts))
  if (!length(cand)) return(NA_character_)
  cand[which.max(rowMeans(counts[cand, , drop = FALSE]))]
}, character(1))
map$n_ids <- vapply(map$symbol, function(g) sum(pm$gene == g), integer(1))
missing <- is.na(map$id_in_data)
if (any(missing)) {
  message("Not found in the data, dropped: ", paste(map$symbol[missing], collapse = ", "))
  note("genes_missing", map$symbol[missing])
  map <- map[!missing, , drop = FALSE]
}

# Drop a gene with too many samples at zero.
map$frac_zero <- rowMeans(counts[map$id_in_data, , drop = FALSE] < 0.5)
low <- map$frac_zero > RD_MAX_ZERO_FRAC
if (any(low)) {
  message("More than ", 100 * RD_MAX_ZERO_FRAC, " percent of the samples at zero, dropped: ",
          paste(sprintf("%s (%.0f%%)", map$symbol[low], 100 * map$frac_zero[low]),
                collapse = ", "))
  note("genes_low_expression", map$symbol[low])
  map <- map[!low, , drop = FALSE]
}
if (nrow(map) < 3L) stop("Fewer than three genes of the set are usable.")

# The samples x genes matrix.
X <- t(logx[map$id_in_data, , drop = FALSE])
colnames(X) <- map$symbol
map$source <- "GENCODE v36 probe map"

# -----------------------------------------------------------------------------
# 5. Covariate adjustment
# -----------------------------------------------------------------------------
# The covariates of each sample. Purity and subtype are matched by the first
# 15 characters of the barcode.
smp <- rownames(X)
s15 <- substr(smp, 1L, 15L)
pur1 <- pur[!is.na(pur$purity), , drop = FALSE]
pur1 <- pur1[!duplicated(pur1$array), , drop = FALSE]
sub1 <- subt[!duplicated(subt$sampleID), , drop = FALSE]
cov <- data.frame(
  age = suppressWarnings(as.numeric(clin[smp, "age_at_index.demographic"])),
  sex = as.character(clin[smp, "gender.demographic"]),
  purity = pur1$purity[match(s15, pur1$array)],
  site = substr(smp, 6L, 7L),
  subtype = as.character(sub1$Subtype_Selected[match(s15, sub1$sampleID)]),
  row.names = smp, stringsAsFactors = FALSE)
# A subtype such as "KIRC.NA" means that no subtype was assigned.
cov$subtype[grepl("\\.NA$", cov$subtype)] <- NA_character_

# Sex: remove a sex with fewer than RD_MIN_SEX_FRAC of the samples. Sex is
# a covariate when two sexes remain.
tab_sex <- table(cov$sex)
minor <- names(tab_sex)[tab_sex / sum(tab_sex) < RD_MIN_SEX_FRAC]
drop_sex <- cov$sex %in% minor
note("removed_minority_sex", sum(drop_sex))
use_sex <- length(unique(stats::na.omit(cov$sex[!drop_sex]))) > 1L

# Subtype: a covariate when it is known for at least RD_MIN_SUBTYPE_FRAC of
# the samples.
subtype_known <- mean(!is.na(cov$subtype[!drop_sex]))
use_subtype <- subtype_known >= RD_MIN_SUBTYPE_FRAC
note("subtype_known_frac", round(subtype_known, 3))
note("subtype_used", use_subtype)

# Remove the samples without a value of a covariate that is used.
need <- c("age", "purity", if (use_sex) "sex", if (use_subtype) "subtype")
drop_na <- !drop_sex & !stats::complete.cases(cov[, need, drop = FALSE])
note("removed_missing_covariates", sum(drop_na))
keep_s <- !drop_sex & !drop_na
X <- X[keep_s, , drop = FALSE]
cov <- cov[keep_s, , drop = FALSE]
cat(sprintf("Covariates: %d samples removed for a minority sex, %d for a missing %s. %d kept.\n",
            sum(drop_sex), sum(drop_na), paste(need, collapse = "/"), nrow(X)))

# Pool the tissue source sites with fewer than RD_MIN_SITE samples.
tab_site <- table(cov$site)
cov$site[cov$site %in% names(tab_site)[tab_site < RD_MIN_SITE]] <- "pooled"
note("n_sites", length(unique(cov$site)))

# The design of the regressions.
design <- data.frame(age = cov$age, purity = cov$purity, row.names = rownames(cov))
if (use_sex) design$sex <- factor(cov$sex)
if (length(unique(cov$site)) > 1L) design$site <- factor(cov$site)
if (use_subtype) design$subtype <- factor(cov$subtype)
stopifnot(!anyNA(design))

# Replace each gene by its residual. With na.action = na.exclude a missing
# value would give a missing residual in its own row. There is none,
# because the samples with missing covariates were removed above.
for (g in colnames(X)) {
  X[, g] <- residuals(lm(X[, g] ~ ., data = design, na.action = na.exclude))
}
stopifnot(!anyNA(X))
note("covariates_removed", paste(names(design), collapse = "+"))
cat(sprintf("Covariates removed: %s\n", paste(names(design), collapse = ", ")))
note("n", nrow(X))
note("p", ncol(X))

# -----------------------------------------------------------------------------
# 6. Centering and scaling
# -----------------------------------------------------------------------------
sdv <- apply(X, 2L, stats::sd)
if (any(sdv == 0)) {
  stop("Genes without variation: ", paste(colnames(X)[sdv == 0], collapse = ", "))
}
X <- scale(X, center = TRUE, scale = TRUE)
attr(X, "scaled:center") <- NULL
attr(X, "scaled:scale") <- NULL
note("centered_scaled", "each gene to mean 0 and standard deviation 1")

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------
# The matrix, the gene map and the quality-control summary.
write.csv(data.frame(sample = rownames(X), X, check.names = FALSE),
          file.path(RD_TAB_DIR, sprintf("expression_%s.csv", RD_NETWORK)), row.names = FALSE)
write.csv(map[, c("symbol", "role", "id_in_data", "n_ids", "frac_zero", "source")],
          file.path(RD_TAB_DIR, sprintf("gene_map_%s.csv", RD_NETWORK)), row.names = FALSE)
write.csv(data.frame(key = names(qc),
                     value = vapply(qc, function(v) paste(v, collapse = ";"), character(1))),
          file.path(RD_TAB_DIR, sprintf("qc_summary_%s.csv", RD_NETWORK)), row.names = FALSE)

# Descriptive normality of each gene across patients.
skew <- function(v) mean((v - mean(v))^3) / sd(v)^3
kurt <- function(v) mean((v - mean(v))^4) / sd(v)^4 - 3
sw <- data.frame(symbol = colnames(X),
                 shapiro_p = vapply(colnames(X), function(g) shapiro.test(X[, g])$p.value, numeric(1)),
                 skewness = vapply(colnames(X), function(g) skew(X[, g]), numeric(1)),
                 excess_kurtosis = vapply(colnames(X), function(g) kurt(X[, g]), numeric(1)),
                 stringsAsFactors = FALSE)
write.csv(sw, file.path(RD_TAB_DIR, sprintf("gene_normality_%s.csv", RD_NETWORK)),
          row.names = FALSE)

# Print a summary.
cat(sprintf("\nWritten: expression_%s.csv (%d patients x %d genes)\n",
            RD_NETWORK, nrow(X), ncol(X)))
cat("Normality of each gene across patients (Shapiro-Wilk p-value, skewness, excess kurtosis):\n")
print(sw, row.names = FALSE, digits = 3)
