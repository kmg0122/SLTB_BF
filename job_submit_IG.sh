#!/bin/bash
#SBATCH -J IG_Sims
#SBATCH -A slt_bf_simulation
#SBATCH -p normal_q
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-23
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
  3) n=500 tau_name="medium" tau_sd=0.4472136 start=27  R=29 ;;
  4) n=500 tau_name="medium" tau_sd=0.4472136 start=30  R=32 ;;
  5) n=500 tau_name="medium" tau_sd=0.4472136 start=33  R=35 ;;
  6) n=500 tau_name="medium" tau_sd=0.4472136 start=36  R=38 ;;
  7) n=500 tau_name="medium" tau_sd=0.4472136 start=39  R=40 ;;
  8) n=500 tau_name="medium" tau_sd=0.4472136 start=41  R=43 ;;
  9) n=500 tau_name="medium" tau_sd=0.4472136 start=44  R=45 ;;
  10) n=500 tau_name="medium" tau_sd=0.4472136 start=46  R=48;;
  11) n=500 tau_name="medium" tau_sd=0.4472136 start=49  R=50;;
  12) n=500 tau_name="big"    tau_sd=0.7071068 start=14  R=16;;
  13) n=500 tau_name="big"    tau_sd=0.7071068 start=17  R=19;;
  13) n=500 tau_name="big"    tau_sd=0.7071068 start=20  R=22;;
  14) n=500 tau_name="big"    tau_sd=0.7071068 start=23  R=25;;
  15) n=500 tau_name="big"    tau_sd=0.7071068 start=27  R=29;;
  16) n=500 tau_name="big"    tau_sd=0.7071068 start=30  R=32;;
  17) n=500 tau_name="big"    tau_sd=0.7071068 start=33  R=35;; 
  18) n=500 tau_name="big"    tau_sd=0.7071068 start=36  R=38;;
  19) n=500 tau_name="big"    tau_sd=0.7071068 start=39  R=40;;
  20) n=500 tau_name="big"    tau_sd=0.7071068 start=41  R=42;;
  21) n=500 tau_name="big"    tau_sd=0.7071068 start=43  R=44;;
  22) n=500 tau_name="big"    tau_sd=0.7071068 start=45  R=47;;
  23) n=500 tau_name="big"    tau_sd=0.7071068 start=48  R=50;;
esac

echo "Starting task ${SLURM_ARRAY_TASK_ID}: n=${n}, tau=${tau_name}, start=${start}, R=${R}" 

Rscript --vanilla simulation_IG.R "${n}" "${tau_name}" "${tau_sd}" "${start}" "${R}"

echo "Finished task ${SLURM_ARRAY_TASK_ID}"
