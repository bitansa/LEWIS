#!/bin/bash
# =============================================================================
# uq_demo.slurm
#
# Uncertainty quantification (Section 7.3 and Appendix H), step 2. Runs LEWIS
# (complexity) on the selected dataset and writes 04_uq/results/tables/. Run
# it after uq_choose_seed.slurm. Less than a minute on one core.
#
# Submit from the root of the repository (the folder that contains
# reproducibility/). Create the folder of the logs once before the first
# job:
#     mkdir -p reproducibility/slurm/logs
#     sbatch reproducibility/slurm/uq_demo.slurm
#
# Without SLURM, run the same command from the root of the repository:
#     Rscript reproducibility/04_uq/uq_decision_rule_demo.R
#
# Edit the partition, modules and resources for your cluster. The values
# below are those used for the paper.
# =============================================================================
#SBATCH --job-name=uq_demo
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
Rscript reproducibility/04_uq/uq_decision_rule_demo.R
