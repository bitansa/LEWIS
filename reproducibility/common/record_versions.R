# =============================================================================
# record_versions.R
#
# Records the versions of R, of the R packages, of the Python packages used
# through reticulate, of the C++ compiler and of MATLAB, and writes them to
#
#     reproducibility/common/versions.txt
#
# Run it in the same environment as the jobs (the same R and Python modules
# and the same RETICULATE_PYTHON). Otherwise reticulate may start another
# Python, and the recorded Python versions would not be those of the runs:
#
#     export RETICULATE_PYTHON="$(command -v python3)"
#     Rscript reproducibility/common/record_versions.R
#
# A file name given as the first argument replaces versions.txt, for
# example for a second machine:
#
#     Rscript reproducibility/common/record_versions.R versions_laptop.txt
# =============================================================================

# ---- Folder of this script --------------------------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})

# ---- Name of the output file ------------------------------------------------
cli <- commandArgs(trailingOnly = TRUE)
out_name <- if (length(cli) >= 1L) cli[[1L]] else "versions.txt"

# ---- R packages and their role ----------------------------------------------
r_pkgs <- c(
  LEWIS      = "the method of the paper",
  Rcpp       = "interface to the C++ code of LEWIS and of Appendix F",
  quantreg   = "least absolute deviation fits of the exact chain (Appendix F)",
  igraph     = "topological orders in the data generation",
  pcalg      = "PC, GES, and the CPDAG of the NOTEARS estimate",
  SID        = "structural intervention distance",
  huge       = "graphical lasso inside the TL method (LEWIS::tl_lingam)",
  dcov       = "distance covariance tests inside the TL method (LEWIS::tl_lingam)",
  reticulate = "interface to the Python packages lingam and notears",
  MASS       = "generalized inverse in the skeleton of Pairwise-LiNGAM",
  data.table = "reading the TCGA files",
  R.utils    = "reading the compressed TCGA files",
  ggplot2    = "figures",
  patchwork  = "figures with several panels",
  scales     = "axis labels of the figures",
  dplyr      = "aggregation",
  tidyr      = "aggregation"
)

# ---- Machine ------------------------------------------------------------------
cpu <- tryCatch({
  info <- readLines("/proc/cpuinfo")
  model <- sub("^model name\\s*:\\s*", "", grep("^model name", info, value = TRUE)[1])
  sprintf("%s (%d logical cores)", model, sum(grepl("^processor", info)))
}, error = function(e) {
  # On macOS, ask sysctl for the processor name.
  out <- tryCatch(system2("sysctl", c("-n", "machdep.cpu.brand_string"),
                          stdout = TRUE, stderr = FALSE),
                  error = function(e) character(0))
  if (length(out) > 0L) out[1] else "unknown"
})

# The date, the machine and R.
lines <- c(sprintf("Recorded: %s", format(Sys.time(), "%Y-%m-%d")),
           sprintf("OS: %s %s", Sys.info()[["sysname"]], Sys.info()[["release"]]),
           sprintf("CPU: %s", cpu),
           sprintf("R: %s", R.version.string),
           sprintf("BLAS: %s", extSoftVersion()[["BLAS"]]),
           "", "R packages:")

# ---- Versions of the R packages ---------------------------------------------
for (pk in names(r_pkgs)) {
  v <- tryCatch(as.character(utils::packageVersion(pk)),
                error = function(e) "not installed")
  lines <- c(lines, sprintf("  %-12s %-12s %s", pk, v, r_pkgs[[pk]]))
}

# ---- Versions of the Python packages ----------------------------------------
lines <- c(lines, "", "Python (through reticulate):")
py <- tryCatch({
  cfg <- reticulate::py_config()
  sprintf("  python       %s", cfg$version)
}, error = function(e) "  reticulate could not start Python")
lines <- c(lines, py)

# The distribution name of each module, used to read its version.
py_mods <- c(lingam = "lingam", notears = "notears", numpy = "numpy",
             scipy = "scipy", sklearn = "scikit-learn")
for (mod in names(py_mods)) {
  v <- tryCatch({
    # Import the module, then read its version from the package metadata.
    reticulate::import(mod, delay_load = FALSE)
    meta <- reticulate::import("importlib.metadata")
    as.character(meta$version(py_mods[[mod]]))
  }, error = function(e) "not installed")
  lines <- c(lines, sprintf("  %-12s %s", py_mods[[mod]], v))
}

# ---- C++ compiler used to build LEWIS ---------------------------------------
lines <- c(lines, "", "C++ compiler:")
cxx <- tryCatch(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX"),
                        stdout = TRUE, stderr = FALSE),
                error = function(e) character(0))
cxxv <- tryCatch({
  drv <- strsplit(trimws(paste(cxx, collapse = " ")), "[[:space:]]+")[[1]][1]
  system2(drv, "--version", stdout = TRUE, stderr = FALSE)[1]
}, error = function(e) NA_character_)
lines <- c(lines,
           sprintf("  CXX          %s", if (length(cxx)) paste(cxx, collapse = " ") else "unknown"),
           sprintf("  version      %s", if (!is.na(cxxv)) cxxv else "unknown"))

# ---- MATLAB (Pairwise-LiNGAM) -----------------------------------------------
lines <- c(lines, "", "MATLAB (Pairwise-LiNGAM):")
mat <- tryCatch(system2("matlab", c("-batch", shQuote("disp(version)")),
                        stdout = TRUE, stderr = FALSE, timeout = 120),
                error = function(e) character(0))
lines <- c(lines, if (length(mat) > 0L) paste0("  ", mat) else
                    "  MATLAB is not on the PATH of this machine")

# ---- Write the file -----------------------------------------------------------
out <- file.path(.here, out_name)
writeLines(lines, out)
cat(paste(lines, collapse = "\n"), "\n\nWritten to ", out, "\n", sep = "")
