#!/bin/bash
#SBATCH -J AR_Sims
#SBATCH -A slt_bf_simulation
#SBATCH -p normal_q
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-5
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kmg0122@vt.edu
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err

mkdir -p logs

module load R-bundle-CRAN/2025.11-foss-2025b

export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

case "${SLURM_ARRAY_TASK_ID}" in
  1) n=500 tau_name="medium" tau_sd=0.4472136 start=11  R=12 ;;
  2) n=500 tau_name="medium" tau_sd=0.4472136 start=13  R=14 ;;
  3) n=500 tau_name="medium" tau_sd=0.4472136 start=15  R=16 ;;
  4) n=500 tau_name="medium" tau_sd=0.4472136 start=17  R=18 ;;
  5) n=500 tau_name="medium" tau_sd=0.4472136 start=19  R=20 ;;
esac

echo "Starting task ${SLURM_ARRAY_TASK_ID}: n=${n}, tau=${tau_name}, start=${start}, R=${R}"

Rscript --vanilla simulation_ARM.R "${n}" "${tau_name}" "${tau_sd}" "${start}" "${R}"

echo "Finished task ${SLURM_ARRAY_TASK_ID}"
