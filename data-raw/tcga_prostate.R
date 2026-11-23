# Creates the data set tcga_prostate (data/tcga_prostate.rda).
#
# Run this script once from the root folder of the package, before
# devtools::document() and R CMD check:
#   source("data-raw/tcga_prostate.R")
#
# The file data-raw/tcga_prostate.csv holds the preprocessed expression of
# the 8 genes of the androgen-receptor program in 462 primary prostate
# tumors of The Cancer Genome Atlas. It is the matrix that LEWIS is fitted
# to in the real-data analysis of the paper. The first column holds the
# sample barcodes and the other columns hold the genes. The preprocessing
# is described in the help page ?tcga_prostate.

# Read the table. check.names = FALSE keeps the gene name NKX3-1 as it is.
raw <- utils::read.csv("data-raw/tcga_prostate.csv", check.names = FALSE,
                       stringsAsFactors = FALSE)

# Check the columns and the number of samples.
genes <- c("AR", "FOXA1", "HOXB13", "KLK3", "KLK2", "TMPRSS2", "NKX3-1", "FKBP5")
stopifnot(identical(colnames(raw), c("sample", genes)),
          nrow(raw) == 462L,
          !anyDuplicated(raw$sample))

# Keep a numeric matrix of the 8 genes with the sample barcodes as row names.
tcga_prostate <- as.matrix(raw[, genes])
storage.mode(tcga_prostate) <- "double"
rownames(tcga_prostate) <- raw$sample

# Check that the values are finite and that every gene is centered and
# scaled to unit standard deviation.
stopifnot(all(is.finite(tcga_prostate)),
          all(abs(colMeans(tcga_prostate)) < 1e-8),
          all(abs(apply(tcga_prostate, 2L, stats::sd) - 1) < 1e-8))

# Save the data set in data/tcga_prostate.rda.
usethis::use_data(tcga_prostate, overwrite = TRUE, compress = "xz")
