#!/bin/bash
#SBATCH -J AR_Sims
#SBATCH -A slt_bf_simulation
#SBATCH -p normal_q
#SBATCH -t 5-00:01:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-9

module load R-bundle-CRAN/2025.11-foss-2025b

case "${SLURM_ARRAY_TASK_ID}" in
  1) n=50  tau_name="small"  tau_sd=0.1 ;;
  2) n=100 tau_name="small"  tau_sd=0.1 ;;
  3) n=500 tau_name="small"  tau_sd=0.1 ;;
  4) n=50  tau_name="medium" tau_sd=0.4472136 ;;
  5) n=100 tau_name="medium" tau_sd=0.4472136 ;;
  6) n=500 tau_name="medium" tau_sd=0.4472136 ;;
  7) n=50  tau_name="big"    tau_sd=0.7071068 ;;
  8) n=100 tau_name="big"    tau_sd=0.7071068 ;;
  9) n=500 tau_name="big"    tau_sd=0.7071068 ;;
esac

echo "Starting task ${SLURM_ARRAY_TASK_ID}: n=${n}, tau=${tau_name}"
Rscript simulation_ARM.R "${n}" "${tau_name}" "${tau_sd}"
