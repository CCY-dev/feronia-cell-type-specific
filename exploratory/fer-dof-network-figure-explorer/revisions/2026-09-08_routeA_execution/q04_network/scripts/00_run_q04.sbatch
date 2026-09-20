#!/bin/bash
# Prototype; not yet carefully vetted.
#SBATCH --job-name=FER_q04_repair
#SBATCH --partition=expansion
#SBATCH --account=nolanlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/resnick/groups/NolanLab/nolan/FER_resubmission_20260908_RouteA/q04_network/logs/q04_%j.out
#SBATCH --error=/resnick/groups/NolanLab/nolan/FER_resubmission_20260908_RouteA/q04_network/logs/q04_%j.err
set -eo pipefail
TASK_ROOT=/resnick/groups/NolanLab/nolan/FER_resubmission_20260908_RouteA/q04_network
export PATH="/resnick/groups/NolanLab/cchen/miniconda3/bin:$PATH"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
printf 'Started: '; date -Is
hostname
/resnick/groups/NolanLab/cchen/miniconda3/envs/cell_oracle/bin/python "$TASK_ROOT/scripts/01_extract_network_metadata.py" "$TASK_ROOT"
/resnick/groups/NolanLab/cchen/miniconda3/envs/muscat/bin/Rscript "$TASK_ROOT/scripts/02_repair_network_annotations.R" "$TASK_ROOT"
/resnick/groups/NolanLab/cchen/miniconda3/envs/muscat/bin/Rscript "$TASK_ROOT/scripts/03_network_sensitivity_figures.R" "$TASK_ROOT"
printf 'Finished: '; date -Is
