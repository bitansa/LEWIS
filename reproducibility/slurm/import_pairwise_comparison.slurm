#!/bin/bash
# =============================================================================
# import_pairwise_comparison.slurm
#
# Comparative study, step 4. Pairwise-LiNGAM, stage 3 of 3. Builds and
# evaluates the estimates from the MATLAB outputs. Writes
# 02_comparison/results/results_pairwise_lingam.csv.
#
# Submit from the root of the repository (the folder that contains
# reproducibility/). Create the folder of the logs once before the first
# job:
#     mkdir -p reproducibility/slurm/logs
#     sbatch reproducibility/slurm/import_pairwise_comparison.slurm
#
# Without SLURM, run the same command from the root of the repository:
#     MC_CORES=8 Rscript reproducibility/02_comparison/pairwise_lingam/import_pairwise_results.R
# MC_CORES is the number of datasets run at the same time. Without it,
# all cores but one are used.
#
# Edit the partition, modules and resources for your cluster. The values
# below are those used for the paper.
# =============================================================================
#SBATCH --job-name=cmp_pw_import
#SBATCH --output=reproducibility/slurm/logs/%x_%j.out
#SBATCH --error=reproducibility/slurm/logs/%x_%j.err
#SBATCH --time=1-00:00:00
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
Rscript reproducibility/02_comparison/pairwise_lingam/import_pairwise_results.R
