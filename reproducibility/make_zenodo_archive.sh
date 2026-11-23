#!/usr/bin/env bash
# =============================================================================
# make_zenodo_archive.sh
#
# Packs the results of each dataset (the raw results of the simulations)
# into one .tar.gz file per study, for deposit on Zenodo. ZENODO.md lists
# what each archive contains.
#
# Usage:
#     bash make_zenodo_archive.sh <reproducibility folder> <output folder>
#
# <reproducibility folder> is a copy of this folder in which the studies
# were run, so that the files listed in ZENODO.md are in place. Only these
# files are read from it.
# <output folder> receives the archives, MD5SUMS.txt (the md5 sum of each
# archive) and FILES_MD5.txt (the md5 sum of each file in the archives).
#
# Each archive holds paths relative to reproducibility/. Unpack it in
# reproducibility/ and every script finds its files.
#
# Users of the repository do not need this script. It documents how the
# archive was made.
#
# The script works on macOS and Linux. It stops at the first error.
# =============================================================================

set -euo pipefail

# ---- Arguments ----------------------------------------------------------------
if [ "$#" -ne 2 ]; then
  echo "Usage: bash make_zenodo_archive.sh <reproducibility folder> <output folder>" >&2
  exit 1
fi
SRC="$(cd "$1" && pwd)"
OUT="$2"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

# A temporary folder for the files to pack. It is removed at the end.
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

# Do not add macOS resource files to the archives.
export COPYFILE_DISABLE=1

# ---- md5 sum of a file, on Linux (md5sum) or macOS (md5) ----------------------
md5_of() {
  if command -v md5sum > /dev/null 2>&1; then
    md5sum "$1" | awk '{print $1}'
  else
    md5 -q "$1"
  fi
}

# ---- Stop if a file is missing --------------------------------------------------
require_file() {
  if [ ! -f "$SRC/$1" ]; then
    echo "Missing file: $1" >&2
    exit 1
  fi
}

# ---- Use the terms of the paper in the header of a staged CSV file -------------
# Column names that contain pcpdag are renamed to rpdag. Only the first
# line can change. The values are not changed.
rename_header() {
  awk 'NR == 1 { gsub(/pcpdag/, "rpdag"); gsub(/PCPDAG/, "RPDAG") } { print }' \
    "$1" > "$1.tmp"
  mv "$1.tmp" "$1"
}

# ---- Copy one CSV file -----------------------------------------------------------
# copy_csv <path relative to the reproducibility folder>
copy_csv() {
  require_file "$1"
  mkdir -p "$(dirname "$STAGE/$1")"
  cp "$SRC/$1" "$STAGE/$1"
  rename_header "$STAGE/$1"
}

# ---- Copy one file unchanged ----------------------------------------------------
copy_file() {
  require_file "$1"
  mkdir -p "$(dirname "$STAGE/$1")"
  cp "$SRC/$1" "$STAGE/$1"
}

# ---- Count the files that match a pattern in a staged folder ---------------------
# check_count <staged folder> <name pattern> <expected number>
check_count() {
  local found
  found=$(find "$STAGE/$1" -type f -name "$2" | wc -l | tr -d ' ')
  if [ "$found" -ne "$3" ]; then
    echo "Expected $3 files $2 in $1, found $found." >&2
    exit 1
  fi
}

# ---- Pack the staged folder of one study, then empty the stage --------------------
pack() {
  local name="$1"
  local archive="$OUT/lewis_raw_${name}.tar.gz"
  # Refuse any file that must not be published: logs, system files and the
  # raw TCGA data.
  if (cd "$STAGE" && find . -type f | grep -E -i '05_real_data|/logs/|\.DS_Store|\.Rhistory|\.RData|data/raw/'); then
    echo "A file that must not be published was staged. Nothing is packed." >&2
    exit 1
  fi
  # Record the md5 sum of each file, then pack.
  (cd "$STAGE" && find . -type f | LC_ALL=C sort | while read -r f; do
     echo "$(md5_of "$f")  ${f#./}"
   done) >> "$OUT/FILES_MD5.txt"
  (cd "$STAGE" && tar -czf "$archive" .)
  echo "$(md5_of "$archive")  $(basename "$archive")" >> "$OUT/MD5SUMS.txt"
  echo "Wrote $archive ($(du -h "$archive" | awk '{print $1}'))"
  rm -rf "${STAGE:?}"/*
}

# Start new md5 lists.
rm -f "$OUT/MD5SUMS.txt" "$OUT/FILES_MD5.txt"

# =============================================================================
# 01_motivation (Figure 1)
# =============================================================================
# The 25 regimes of Figure 1.
REGIMES="all_t5 gaussian t15"
for k in 01 02 03 04 05 06 07 08 09 10 11; do
  REGIMES="$REGIMES count_t5g_k$k count_t5t15_k$k"
done

# Results of the ten methods. Only the rows of the 25 regimes are packed.
# The regime is the first column, in quotes.
require_file "01_motivation/results/motivation_results.csv"
mkdir -p "$STAGE/01_motivation/results"
awk -F',' -v keep="$REGIMES" '
  BEGIN { n = split(keep, r, " "); for (i = 1; i <= n; i++) ok["\"" r[i] "\""] = 1 }
  NR == 1 || ($1 in ok) { print }
' "$SRC/01_motivation/results/motivation_results.csv" \
  > "$STAGE/01_motivation/results/motivation_results.csv"
rename_header "$STAGE/01_motivation/results/motivation_results.csv"

# Likelihood-ratio matrices of Pairwise-LiNGAM for the 25 regimes. The
# pattern starts with the regime name, so "t15" does not match
# "count_t5t15_k01".
mkdir -p "$STAGE/01_motivation/pairwise_lingam/outputs"
for rg in $REGIMES; do
  for n in 100 250 500; do
    for f in "$SRC/01_motivation/pairwise_lingam/outputs/${rg}_n${n}_r"[0-9][0-9][0-9]"_LR.csv"; do
      cp "$f" "$STAGE/01_motivation/pairwise_lingam/outputs/"
    done
  done
done
check_count "01_motivation/pairwise_lingam/outputs" "*_LR.csv" 7500

# Run times of Pairwise-LiNGAM for the 25 regimes (first column, no quotes).
require_file "01_motivation/pairwise_lingam/outputs/runtimes.csv"
awk -F',' -v keep="$REGIMES" '
  BEGIN { n = split(keep, r, " "); for (i = 1; i <= n; i++) ok[r[i]] = 1 }
  NR == 1 || ($1 in ok) { print }
' "$SRC/01_motivation/pairwise_lingam/outputs/runtimes.csv" \
  > "$STAGE/01_motivation/pairwise_lingam/outputs/runtimes.csv"
pack "01_motivation"

# =============================================================================
# 02_comparison (Section 7.2, Appendices E.1, G.1 and G.3)
# =============================================================================
# Results of the ten methods.
for m in lewis_complexity lewis_uniform pc ges lingam direct_lingam rcd \
         tl_zhao notears pairwise_lingam; do
  copy_csv "02_comparison/results/results_${m}.csv"
done
# The list of the datasets exported for MATLAB.
copy_file "02_comparison/pairwise_lingam/manifest.csv"
# Likelihood-ratio matrices of Pairwise-LiNGAM and their run times.
require_file "02_comparison/pairwise_lingam/outputs/runtimes.csv"
mkdir -p "$STAGE/02_comparison/pairwise_lingam/outputs"
find "$SRC/02_comparison/pairwise_lingam/outputs" -maxdepth 1 -type f \
  \( -name '*_LR.csv' -o -name 'runtimes.csv' \) \
  -exec cp {} "$STAGE/02_comparison/pairwise_lingam/outputs/" \;
check_count "02_comparison/pairwise_lingam/outputs" "*_LR.csv" 5400
pack "02_comparison"

# =============================================================================
# 03_calibration (remark on alpha and d in Section 5, Appendix G.2)
# =============================================================================
copy_csv "03_calibration/results/d_sweep_results.csv"
copy_csv "03_calibration/results/alpha_sweep_results.csv"
pack "03_calibration"

# =============================================================================
# 04_uq (Section 7.3 and Appendix H)
# =============================================================================
# The retained draws of the selected dataset.
copy_csv "04_uq/results/tables/samples_per_dag.csv"
# The tables of each of the 10 datasets of the seed rule.
for s in 01 02 03 04 05 06 07 08 09 10; do
  for f in "$SRC/04_uq/results/seed_sweep/seed_${s}/tables/"*.csv; do
    copy_csv "04_uq/results/seed_sweep/seed_${s}/tables/$(basename "$f")"
  done
done
pack "04_uq"

# =============================================================================
# 06_laplace_vs_exact (Appendix F)
# =============================================================================
copy_csv "06_laplace_vs_exact/results/exact_vs_laplace.csv"
pack "06_laplace_vs_exact"

# ---- Report -------------------------------------------------------------------
echo
echo "Archives and md5 sums are in $OUT:"
cat "$OUT/MD5SUMS.txt"
echo "The md5 sum of every packed file is in FILES_MD5.txt."
