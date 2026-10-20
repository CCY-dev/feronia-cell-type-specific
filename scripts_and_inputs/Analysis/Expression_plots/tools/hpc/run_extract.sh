#SBATCH --account=nolanlab
#SBATCH --partition=expansion
#SBATCH --mem=128G
#SBATCH --cpus-per-task=2
#SBATCH --time=01:30:00
# Read-only export from the shared FER Seurat v3 object; writes only into this job's workdir under /resnick/scratch/nolan.
set -eo pipefail
ENV=/resnick/groups/NolanLab/cchen/miniconda3/envs/scKB_copilot
echo "host: $(hostname)  start: $(date -Is)  env: $ENV"
"$ENV/bin/Rscript" --vanilla extract_si_expression.R 2>&1 | tee si_expression_extract.log
sha256sum si_expression_cells.csv.gz | tee si_expression_cells.csv.gz.sha256
sha256sum extract_si_expression.R >> si_expression_cells.csv.gz.sha256
echo "end: $(date -Is)" | tee -a si_expression_extract.log
sacct -j "${SLURM_JOB_ID:-0}" --format=JobID,Elapsed,MaxRSS,State 2>/dev/null | tee -a si_expression_extract.log || true
ls -lh si_expression_* || true
