#!/bin/bash
# =============================================================================
# comparison_direct_lingam.slurm
#
# Comparative study (Section 7.2), step 1. Runs DirectLiNGAM (Python package
# lingam) on the 100 datasets of each of the 54 settings. Writes
# 02_comparison/results/results_direct_lingam.csv. About 40 CPU hours.
#
# Submit from the root of the repository (the folder that contains
# reproducibility/). Create the folder of the logs once before the first
# job:
#     mkdir -p reproducibility/slurm/logs
#     sbatch reproducibility/slurm/comparison_direct_lingam.slurm
#
# Without SLURM, run the same command from the root of the repository:
#     MC_CORES=8 Rscript reproducibility/02_comparison/run_one_method.R direct_lingam
# MC_CORES is the number of datasets run at the same time. Without it,
# all cores but one are used.
#
# Edit the partition, modules and resources for your cluster. The values
# below are those used for the paper.
# =============================================================================
#SBATCH --job-name=cmp_direct_lingam
#SBATCH --output=reproducibility/slurm/logs/%x_%j.out
#SBATCH --error=reproducibility/slurm/logs/%x_%j.err
#SBATCH --time=21-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=128G
# To choose a partition, remove the first # of the next line and edit it.
##SBATCH --partition=<partition>

# Stop at the first error.
set -eo pipefail

# Load R. Edit the module name for your cluster.
module load R
# Load the Python that has the packages lingam and notears, and point
# reticulate to it. Edit the module name for your cluster.
module load Python
export RETICULATE_PYTHON="${RETICULATE_PYTHON:-$(command -v python3)}"

# Run one dataset per core, with single-threaded linear algebra.
export MC_CORES="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Run the script from the root of the repository.
cd "$SLURM_SUBMIT_DIR"
Rscript reproducibility/02_comparison/run_one_method.R direct_lingam
