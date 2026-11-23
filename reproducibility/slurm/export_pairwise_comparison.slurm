#!/bin/bash
# =============================================================================
# export_pairwise_comparison.slurm
#
# Comparative study, step 2. Pairwise-LiNGAM, stage 1 of 3. Writes the 5,400
# datasets to 02_comparison/pairwise_lingam/inputs/ and their list to
# manifest.csv, for the MATLAB stage (see 02_comparison/README.md). About 1.9
# GB.
#
# Submit from the root of the repository (the folder that contains
# reproducibility/). Create the folder of the logs once before the first
# job:
#     mkdir -p reproducibility/slurm/logs
#     sbatch reproducibility/slurm/export_pairwise_comparison.slurm
#
# Without SLURM, run the same command from the root of the repository:
#     Rscript reproducibility/02_comparison/pairwise_lingam/export_datasets_for_matlab.R
#
# Edit the partition, modules and resources for your cluster. The values
# below are those used for the paper.
# =============================================================================
#SBATCH --job-name=cmp_pw_export
#SBATCH --output=reproducibility/slurm/logs/%x_%j.out
#SBATCH --error=reproducibility/slurm/logs/%x_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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
Rscript reproducibility/02_comparison/pairwise_lingam/export_datasets_for_matlab.R
