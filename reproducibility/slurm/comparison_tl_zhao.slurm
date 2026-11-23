#!/bin/bash
# =============================================================================
# comparison_tl_zhao.slurm
#
# Comparative study (Section 7.2), step 1. Runs TL on the 100 datasets
# of each of the 54 settings. TL is LEWIS::tl_lingam, a copy of the TL code
# of the R package TransGraph 1.1.0. Writes
# 02_comparison/results/results_tl_zhao.csv. About 80 CPU hours.
#
# Submit from the root of the repository (the folder that contains
# reproducibility/). Create the folder of the logs once before the first
# job:
#     mkdir -p reproducibility/slurm/logs
#     sbatch reproducibility/slurm/comparison_tl_zhao.slurm
#
# Without SLURM, run the same command from the root of the repository:
#     MC_CORES=8 Rscript reproducibility/02_comparison/run_one_method.R tl_zhao
# MC_CORES is the number of datasets run at the same time. Without it,
# all cores but one are used.
#
# Edit the partition, modules and resources for your cluster. The values
# below are those used for the paper.
# =============================================================================
#SBATCH --job-name=cmp_tl_zhao
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
Rscript reproducibility/02_comparison/run_one_method.R tl_zhao
