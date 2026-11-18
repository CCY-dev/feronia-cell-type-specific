#!/bin/bash
#SBATCH --job-name=20261020_CO_seeded
#SBATCH --array=0-10  # 0-9 = the 10 cell-type groups; 10 = columella replicate (determinism test vs task 9)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20  # matches n_jobs=20 in knn_imputation
#SBATCH --mem=128G          # the 20260317 run peaked at ~26 GB per task (sacct MaxRSS)
#SBATCH --time=48:00:00     # the 20260317 run took 1.5-6.2 h per task
# NOTE: lab cluster path below -- change it to your own path
#SBATCH --output=/resnick/groups/NolanLab/cchen/Feronia/claude_script/logs/CO_seeded_%A_%a.o
# NOTE: lab cluster path below -- change it to your own path
#SBATCH --error=/resnick/groups/NolanLab/cchen/Feronia/claude_script/logs/CO_seeded_%A_%a.e

# Rerun of the 20260317 CellOracle GRNs with: per-task working folders (no shared ./network_analysis/),
# tissue-only phloem/xylem DE features, a fixed hash seed, and a self-check of every score table.

# NOTE: lab cluster path below -- change it to your own path
export PATH="/resnick/groups/NolanLab/cchen/miniconda3/bin:$PATH"
source activate cell_oracle

# Fixed hash seed: without it the GRN coefficients differ between runs
export PYTHONHASHSEED=0
# Per-task genomepy cache to avoid collisions between parallel array tasks
export XDG_CACHE_HOME="/tmp/genomepy_cache_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

echo "Task $SLURM_ARRAY_TASK_ID started at $(date) on $(hostname)"
# NOTE: lab cluster path below -- change it to your own path
python /resnick/groups/NolanLab/cchen/Feronia/claude_script/20261020_celloracle_grn_seeded_isolated.py $SLURM_ARRAY_TASK_ID \
  # NOTE: lab cluster path below -- change it to your own path
  2> /resnick/groups/NolanLab/cchen/Feronia/claude_script/logs/CO_seeded_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_python.err
exit_code=$?
echo "Task $SLURM_ARRAY_TASK_ID finished with exit code $exit_code at $(date)"
exit $exit_code
