#!/bin/bash
# =============================================================================
# lewis_runtime.slurm
#
# Computing times of Appendix E.1. Runs LEWIS (complexity) with the package
# on 1,200 datasets of the comparative study and writes
# 02_comparison/results/lewis_runtime_complexity.csv. It checks the SHD
# against results_lewis_complexity.csv, so run it after
# comparison_lewis_complexity.slurm. About 70 CPU hours.
#
# Submit from the root of the repository (the folder that contains
# reproducibility/). Create the folder of the logs once before the first
# job:
#     mkdir -p reproducibility/slurm/logs
#     sbatch reproducibility/slurm/lewis_runtime.slurm
#
# Without SLURM, run the same command from the root of the repository:
#     MC_CORES=8 Rscript reproducibility/02_comparison/measure_lewis_runtime.R
# MC_CORES is the number of datasets run at the same time. Without it,
# all cores but one are used.
#
# Edit the partition, modules and resources for your cluster. The values
# below are those used for the paper.
# =============================================================================
#SBATCH --job-name=lewis_runtime
#SBATCH --output=reproducibility/slurm/logs/%x_%j.out
#SBATCH --error=reproducibility/slurm/logs/%x_%j.err
#SBATCH --time=7-00:00:00
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

# Run one dataset per core, with single-threaded linear algebra.
export MC_CORES="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Run the script from the root of the repository.
cd "$SLURM_SUBMIT_DIR"
Rscript reproducibility/02_comparison/measure_lewis_runtime.R
