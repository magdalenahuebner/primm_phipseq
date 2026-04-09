#!/bin/bash
#SBATCH --job-name=primm_analysis
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=40G

set -euo pipefail

# Always run from the directory where you submitted the job
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

echo "Working directory: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

module purge
module load R
module load Pandoc/3.9

# Run analysis (project .Rprofile will be picked up from the working directory)
Rscript --no-save --no-restore R/02-analysis.R \
  N_CORES="${SLURM_CPUS_PER_TASK}" \
  MAX_GB=40 \
  LOG=true \
  FORCE=true
  