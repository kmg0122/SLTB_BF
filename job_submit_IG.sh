#!/bin/bash
#SBATCH -J IG_Sims
#SBATCH -A slt_bf_simulation
#SBATCH -p normal_q
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kmg0122@vt.edu
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err

mkdir -p logs

module load R-bundle-CRAN/2025.11-foss-2025b

export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

case "${SLURM_ARRAY_TASK_ID}" in
  1) n=500 tau_name="medium" tau_sd=0.4472136 start=11 R=15 ;;
  2) n=500 tau_name="medium" tau_sd=0.4472136 start=16 R=19 ;;
  3) n=500 tau_name="medium" tau_sd=0.4472136 start=27  R=30 ;;
  4) n=500 tau_name="medium" tau_sd=0.4472136 start=31  R=35 ;;
  5) n=500 tau_name="medium" tau_sd=0.4472136 start=36  R=40 ;;
  6) n=500 tau_name="medium" tau_sd=0.4472136 start=41  R=45 ;;
  7) n=500 tau_name="medium" tau_sd=0.4472136 start=46  R=50 ;;
  8) n=500 tau_name="big"    tau_sd=0.7071068 start=14  R=19;;
  9) n=500 tau_name="big"    tau_sd=0.7071068 start=30  R=35;;
  10) n=500 tau_name="big"    tau_sd=0.7071068 start=36  R=40;;
  11) n=500 tau_name="big"    tau_sd=0.7071068 start=41  R=45;;
  12) n=500 tau_name="big"    tau_sd=0.7071068 start=46  R=50;;
esac

echo "Starting task ${SLURM_ARRAY_TASK_ID}: n=${n}, tau=${tau_name}, start=${start}, R=${R}" 

Rscript --vanilla simulation_HC.R "${n}" "${tau_name}" "${tau_sd}" "${start}" "${R}"

echo "Finished task ${SLURM_ARRAY_TASK_ID}"
