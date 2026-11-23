#!/bin/bash
# =============================================================================
# calibration_d_sweep.slurm
#
# Calibration study (Appendix G.2), sweep over d with alpha = 0.75: 84
# settings, 10 values of d and 100 datasets. Writes
# 03_calibration/results/d_sweep_results.csv. About 8,500 CPU hours in the
# runs of the paper, that is about 4 days on 100 cores.
#
# Submit from the root of the repository (the folder that contains
# reproducibility/). Create the folder of the logs once before the first
# job:
#     mkdir -p reproducibility/slurm/logs
#     sbatch reproducibility/slurm/calibration_d_sweep.slurm
#
# Without SLURM, run the same command from the root of the repository:
#     MC_CORES=8 Rscript reproducibility/03_calibration/calibration_sweep.R d
# MC_CORES is the number of datasets run at the same time. Without it,
# all cores but one are used.
#
# Edit the partition, modules and resources for your cluster. The values
# below are those used for the paper.
# =============================================================================
#SBATCH --job-name=cal_d
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

# Run one dataset per core, with single-threaded linear algebra.
export MC_CORES="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Run the script from the root of the repository.
cd "$SLURM_SUBMIT_DIR"
Rscript reproducibility/03_calibration/calibration_sweep.R d
