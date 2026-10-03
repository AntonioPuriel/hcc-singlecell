#!/bin/bash
#SBATCH --job-name=hcc_all
#SBATCH --partition=standard
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=46G
#SBATCH --output=logs/run_all_%j.out
#SBATCH --error=logs/run_all_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=apurielh@gmail.com

# One launch for the whole project:
#   1. checks every R package the pipeline uses and installs any that is missing
#   2. validates the pipeline definition (stops here if _targets.R is broken)
#   3. runs targets::tar_make(); only outdated steps are rebuilt, and a failing
#      step does not stop independent steps (error = "continue" in _targets.R)
#   4. reports the status of every step, the result tables and storage
# The last step writes the HTML report to docs/index.html (GitHub Pages).
# Launch from the project root:  sbatch slurm/run_all.sh
# Report: pipeline_report.txt. Exit status 0 only if no step errored.
# Memory: standard enforces MaxMemPerCPU=2000 MB, so memory comes with CPUs (24 x 2000 MB).

set +e
PROJ=/work/aphernand001/hcc-singlecell
ENV="$PROJ/env"
MAX_GB=60
cd "$PROJ" || exit 1

REPORT="$PROJ/pipeline_report.txt"
: > "$REPORT"
exec >> "$REPORT" 2>&1
section() { echo; echo "=================== $1 ==================="; date; }

export MAMBA_ROOT_PREFIX="$PROJ/mamba"
export TMPDIR="$PROJ/tmp"
export OMP_NUM_THREADS=1            # avoid BLAS oversubscription with parallel workers
export EXPERIMENT_HUB_CACHE="$PROJ/cache/experimenthub"
export ANNOTATION_HUB_CACHE="$PROJ/cache/annotationhub"
mkdir -p "$TMPDIR" logs "$PROJ/cache"
MM="$PROJ/tools/bin/micromamba"
RUN() { "$MM" run -p "$ENV" "$@"; }

section "0. STORAGE CHECK (limit ${MAX_GB} GB)"
used_gb=$(du -s --block-size=1G "$PROJ" | cut -f1)
echo "Project size: ${used_gb} GB"
if [ "$used_gb" -ge "$MAX_GB" ]; then
  echo "ERROR: project exceeds ${MAX_GB} GB. Clean up _targets/ or tmp/ before running."
  exit 1
fi

section "1. GIT VERSION"
git log -1 --format='%h %ad %s' --date=short 2>/dev/null || echo "(not a git checkout)"

section "2. R PACKAGES (missing ones are installed)"
if ! RUN Rscript -e 'quit(status = !rmarkdown::pandoc_available())' 2>/dev/null; then
  echo "pandoc missing: installing from conda-forge (existing packages frozen)"
  "$MM" install -y -p "$ENV" -c conda-forge --freeze-installed pandoc
fi
RUN Rscript -e '
options(Ncpus = 4, timeout = 900)
github <- c(presto = "immunogenomics/presto", liana = "saezlab/liana", msigdbrdata = "igordot/msigdbdf")
needed <- c("Seurat", "SeuratObject", "harmony", "igraph", "scDblFinder", "SingleR", "celldex",
            "infercnv", "DESeq2", "fgsea", "msigdbr", "slingshot", "SingleCellExperiment",
            "targets", "tarchetypes", "qs2", "rmarkdown", "knitr", "data.table", "ggplot2",
            "patchwork", names(github), "CytoTRACE2")
ok <- function(p) requireNamespace(p, quietly = TRUE)
missing <- needed[!vapply(needed, ok, logical(1))]
cat("Missing:", if (length(missing)) missing else "none", "\n")
for (p in intersect(missing, names(github)))
  remotes::install_github(github[[p]], lib = .Library, upgrade = "never", build_vignettes = FALSE)
other <- setdiff(missing, c(names(github), "CytoTRACE2"))
if (length(other)) BiocManager::install(other, lib = .Library, update = FALSE, ask = FALSE)
if ("CytoTRACE2" %in% missing) cat("CytoTRACE2 missing: install with slurm/06_install_cytotrace2.sh (needs the HiClimR workaround)\n")
for (p in needed) cat(sprintf("%-22s %s\n", p, if (ok(p)) as.character(packageVersion(p)) else "MISSING"))
'

section "3. VALIDATE PIPELINE"
if ! RUN Rscript -e 'targets::tar_validate(); cat("pipeline definition OK\n")'; then
  echo "ERROR: _targets.R does not validate; nothing was run."
  exit 1
fi

section "4. OUTDATED TARGETS"
RUN Rscript -e 'print(targets::tar_outdated())'

section "5. tar_make()"
RUN Rscript -e 'targets::tar_make(reporter = "timestamp")'
echo "tar_make exit=$?"

section "6. STATUS OF EVERY TARGET"
RUN Rscript -e '
m <- targets::tar_meta(fields = c("name","seconds","bytes","error","warnings"))
m$minutes <- round(m$seconds / 60, 1); m$MB <- round(m$bytes / 1e6, 1)
print(as.data.frame(m[, c("name","minutes","MB","warnings","error")]), right = FALSE)
err <- m$name[!is.na(m$error)]
cat("\nERRORED TARGETS:", if (length(err)) paste(err, collapse = ", ") else "none", "\n")
writeLines(err, file.path(Sys.getenv("TMPDIR"), "errored.txt"))'
N_ERR=$(grep -c . "$TMPDIR/errored.txt" 2>/dev/null); N_ERR=${N_ERR:-0}

section "7. RESULT TABLES"
for f in results/tables/*.tsv; do echo "--- $f"; column -t -s $'\t' "$f" | head -40; done
ls -l docs/index.html 2>/dev/null || echo "docs/index.html not produced"

section "8. STORAGE AFTER RUN"
du -sh "$PROJ/_targets" "$PROJ" 2>/dev/null
section "FIN"
echo "errored targets: $N_ERR"
[ "$N_ERR" -eq 0 ]
