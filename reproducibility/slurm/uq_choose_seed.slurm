#!/bin/bash
# =============================================================================
# uq_choose_seed.slurm
#
# Uncertainty quantification (Section 7.3 and Appendix H), step 1. Applies
# the rule of Appendix H to the datasets of seeds 1 to 10. Writes
# 04_uq/results/seed_sweep/. A few minutes on one core.
#
# Submit from the root of the repository (the folder that contains
# reproducibility/). Create the folder of the logs once before the first
# job:
#     mkdir -p reproducibility/slurm/logs
#     sbatch reproducibility/slurm/uq_choose_seed.slurm
#
# Without SLURM, run the same command from the root of the repository:
#     Rscript reproducibility/04_uq/choose_seed_for_demo.R
#
# Edit the partition, modules and resources for your cluster. The values
# below are those used for the paper.
# =============================================================================
#SBATCH --job-name=uq_choose_seed
#SBATCH --output=reproducibility/slurm/logs/%x_%j.out
#SBATCH --error=reproducibility/slurm/logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# To choose a partition, remove the first # of the next line and edit it.
##SBATCH --partition=<partition>

# Stop at the first error.
set -eo pipefail

# Load R. Edit the module name for your cluster.
module load R

# Use one core, with single-threaded linear algebra.
export MC_CORES="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Run the script from the root of the repository.
cd "$SLURM_SUBMIT_DIR"
Rscript reproducibility/04_uq/choose_seed_for_demo.R
