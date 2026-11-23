#!/bin/bash
#SBATCH --job-name=brunello10k-score-only
#SBATCH --partition=workers
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=36:00:00
# Supply --output ROOT/score-%j.log when submitting. No array/cohort entry point.
set -euo pipefail
ROOT=${1:?usage: sbatch script ROOT [work|finalize SCORE_JOB]}
MODE=${2:-work}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 TF_NUM_INTRAOP_THREADS=1 TF_NUM_INTEROP_THREADS=1 CUDA_VISIBLE_DEVICES=-1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$ROOT/tmp-$SLURM_JOB_ID"
mkdir "$TMPDIR"
if [[ "$MODE" == work ]]; then
    python3 -B "$ROOT/tooling/supervise.py" --limit-kib 4194304 --interval 0.5 \
        --log "$ROOT/supervisor-$SLURM_JOB_ID.log" --summary "$ROOT/supervisor-$SLURM_JOB_ID.json" \
        -- python3 -B "$ROOT/tooling/scoring_completion.py" work --root "$ROOT"
elif [[ "$MODE" == finalize ]]; then
    # Submit separately with --dependency=afterok:SCORE_JOB and a short time limit.
    python3 -B "$ROOT/tooling/scoring_completion.py" finalize --root "$ROOT" --score-job "${3:?scoring job required}"
else
    echo "Invalid mode: $MODE" >&2; exit 2
fi
# Only remove this invocation's empty private temporary directory. Keep failures.
rmdir "$TMPDIR" 2>/dev/null || true
