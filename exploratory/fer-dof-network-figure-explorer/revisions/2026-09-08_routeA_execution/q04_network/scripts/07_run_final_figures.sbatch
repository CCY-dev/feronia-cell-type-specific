#!/bin/bash
# Prototype; not yet carefully vetted.
#SBATCH --job-name=FER_q04_figures
#SBATCH --partition=expansion
#SBATCH --account=nolanlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=/resnick/groups/NolanLab/nolan/FER_resubmission_20260908_RouteA/q04_network/logs/q04_figures_%j.out
#SBATCH --error=/resnick/groups/NolanLab/nolan/FER_resubmission_20260908_RouteA/q04_network/logs/q04_figures_%j.err
set -eo pipefail
TASK_ROOT=/resnick/groups/NolanLab/nolan/FER_resubmission_20260908_RouteA/q04_network
export OMP_NUM_THREADS=2
/resnick/groups/NolanLab/cchen/miniconda3/envs/muscat/bin/Rscript "$TASK_ROOT/scripts/03_network_sensitivity_figures.R" "$TASK_ROOT"
