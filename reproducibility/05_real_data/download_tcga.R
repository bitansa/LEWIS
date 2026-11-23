# =============================================================================
# download_tcga.R
#
# Downloads the public TCGA files of the real-data analysis into data/raw/
# and checks each file against the md5 sum in data/raw_files.csv. The md5
# sums are those of the files used for the paper. The files are not in the
# repository or in the Zenodo archive. The address of the subtype file
# redirects to another server. download.file() follows the redirect.
#
# The files (data/raw_files.csv):
#   TCGA-<cohort>.star_counts.tsv.gz   gene counts, log2(count + 1), STAR
#                                      pipeline, GENCODE v36 (GDC data
#                                      release 40), UCSC Xena GDC hub
#   TCGA-<cohort>.clinical.tsv.gz      sample type, preservation, age, sex
#                                      (GDC data release 41), UCSC Xena GDC hub
#   gencode.v36.annotation.gtf.gene.probemap
#                                      gene identifier to gene symbol, UCSC
#                                      Xena GDC hub
#   TCGA_mastercalls.abs_tables_JSedit.fixed.txt
#                                      tumor purity by ABSOLUTE, Pan-Cancer
#                                      Atlas supplemental data (GDC)
#   merged_sample_quality_annotations.tsv
#                                      sample quality flags, Pan-Cancer Atlas
#                                      supplemental data (GDC)
#   TCGASubtype.20170308.tsv.gz        molecular subtypes, UCSC Xena
#                                      Pan-Cancer Atlas hub
# The column needed_for is "reported" for the files of the two reported
# analyses (prostate PRAD and liver LIHC, and the files shared by all
# cohorts) and "screen" for the files that only the screen of the 18 gene
# sets needs (the cohorts BRCA, KIRC, LUAD and STAD).
#
# The UCSC Xena hubs may replace a file by a newer data release. The md5
# check then fails and the script stops. Set TCGA_ALLOW_MD5_MISMATCH=1 to
# keep such a file anyway. The results may then differ from the paper.
#
# Run:  Rscript download_tcga.R          (the files of the reported analyses)
#       Rscript download_tcga.R all      (also the files of the screen)
# The total size is about 115 MB for the reported analyses and 450 MB for
# all files.
# =============================================================================

# ---- Folder of this script --------------------------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})

# ---- Which files --------------------------------------------------------------
cli <- commandArgs(trailingOnly = TRUE)
ALL <- length(cli) >= 1L && identical(cli[[1L]], "all")
ALLOW_MISMATCH <- identical(Sys.getenv("TCGA_ALLOW_MD5_MISMATCH", "0"), "1")
manifest <- read.csv(file.path(.here, "data", "raw_files.csv"), stringsAsFactors = FALSE)
if (!ALL) manifest <- manifest[manifest$needed_for == "reported", ]
RAW_DIR <- file.path(.here, "data", "raw")
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

# Large files need a long timeout. The default of R is 60 seconds.
options(timeout = max(3600, getOption("timeout")))

# The md5 sum of a file.
md5_of <- function(path) unname(tools::md5sum(path))

# Download one file to a temporary name, with up to `tries` attempts, and
# return the temporary path.
fetch <- function(url, dest, tries = 3L) {
  tmp <- paste0(dest, ".part")
  for (k in seq_len(tries)) {
    ok <- tryCatch({
      utils::download.file(url, tmp, mode = "wb", quiet = FALSE)
      TRUE
    }, error = function(e) {
      message(sprintf("Attempt %d of %d failed for %s: %s", k, tries, url,
                      conditionMessage(e)))
      FALSE
    })
    if (ok && file.exists(tmp) && file.size(tmp) > 0) return(tmp)
    Sys.sleep(5 * k)
  }
  stop("Could not download ", url, " after ", tries, " attempts.")
}

# ---- Download and check each file ---------------------------------------------
problems <- character(0)
for (r in seq_len(nrow(manifest))) {
  f <- manifest$file[r]
  dest <- file.path(RAW_DIR, f)
  # Skip a file that is present with the right md5 sum.
  if (file.exists(dest) && identical(md5_of(dest), manifest$md5[r])) {
    cat(sprintf("Present and checked: %s\n", f))
    next
  }
  cat(sprintf("Downloading %s (%.1f MB)\n", f, manifest$bytes[r] / 1e6))
  tmp <- fetch(manifest$url[r], dest)
  # Check the md5 sum before the file takes its final name.
  got <- md5_of(tmp)
  if (!identical(got, manifest$md5[r])) {
    msg <- sprintf(paste0("The md5 sum of %s is %s, not %s. The source may hold a newer ",
                          "data release than the one used in the paper."),
                   f, got, manifest$md5[r])
    if (!ALLOW_MISMATCH) {
      unlink(tmp)
      stop(msg, " Set TCGA_ALLOW_MD5_MISMATCH=1 to keep the file anyway.")
    }
    warning(msg, call. = FALSE)
    problems <- c(problems, f)
  }
  file.rename(tmp, dest)
  cat(sprintf("Saved: %s\n", dest))
}

# ---- Report -------------------------------------------------------------------
cat(sprintf("\n%d files in %s.\n", nrow(manifest), RAW_DIR))
if (length(problems) > 0L) {
  cat("These files differ from the files of the paper:", paste(problems, collapse = ", "), "\n")
} else {
  cat("Every file matches the md5 sum of the paper.\n")
}
