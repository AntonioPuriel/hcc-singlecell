#!/bin/bash
#SBATCH --job-name=hcc_pipeline
#SBATCH --partition=standard
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err

# Runs (or resumes) the targets pipeline. Only outdated targets are rebuilt,
# so the same command is used after a crash or after changing code.
# Launch from the project root:  sbatch slurm/10_run_pipeline.sh
# Report: pipeline_report.txt (progress, errors, storage)

set +e
PROJ=/work/aphernand001/hcc-singlecell
MAX_GB=60
cd "$PROJ" || exit 1

REPORT="$PROJ/pipeline_report.txt"
: > "$REPORT"
exec >> "$REPORT" 2>&1
section() { echo; echo "=================== $1 ==================="; date; }

export MAMBA_ROOT_PREFIX="$PROJ/mamba"
export TMPDIR="$PROJ/tmp"
export OMP_NUM_THREADS=1            # avoid BLAS oversubscription with parallel workers
mkdir -p "$TMPDIR" logs
MM="$PROJ/tools/bin/micromamba"

section "0. STORAGE CHECK (limit ${MAX_GB} GB)"
used_gb=$(du -s --block-size=1G "$PROJ" | cut -f1)
echo "Project size: ${used_gb} GB"
if [ "$used_gb" -ge "$MAX_GB" ]; then
  echo "ERROR: project exceeds ${MAX_GB} GB. Clean up _targets/ or tmp/ before running."
  exit 1
fi

section "1. GIT VERSION"
git log -1 --format='%h %ad %s' --date=short 2>/dev/null || echo "(not a git checkout)"

section "2. OUTDATED TARGETS"
"$MM" run -p "$PROJ/env" Rscript -e 'print(targets::tar_outdated())'

section "3. tar_make()"
"$MM" run -p "$PROJ/env" Rscript -e 'targets::tar_make(reporter = "timestamp")'
MAKE_EXIT=$?
echo "exit=$MAKE_EXIT"

section "4. STATUS OF EVERY TARGET"
"$MM" run -p "$PROJ/env" Rscript -e '
m <- targets::tar_meta(fields = c("name","seconds","bytes","error","warnings"))
m$minutes <- round(m$seconds / 60, 1); m$MB <- round(m$bytes / 1e6, 1)
print(as.data.frame(m[, c("name","minutes","MB","warnings","error")]), right = FALSE)'

section "5. RESULT TABLES"
for f in results/tables/*.tsv; do echo "--- $f"; column -t -s $'\t' "$f" | head -40; done

section "6. STORAGE AFTER RUN"
du -sh "$PROJ/_targets" "$PROJ" 2>/dev/null
section "FIN"
exit $MAKE_EXIT
