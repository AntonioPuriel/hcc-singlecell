#!/bin/bash
#SBATCH --job-name=hcc_cytotrace2
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/cytotrace2_%j.out
#SBATCH --error=logs/cytotrace2_%j.err

# Installs CytoTRACE2 (GitHub, R version) and its CRAN dependencies into env/,
# then runs it on a small synthetic matrix to check that it works.
# Does not touch packages used by earlier pipeline steps.
# Launch from the project root:  sbatch slurm/06_install_cytotrace2.sh
# Report: cytotrace2_report.txt

set +e
PROJ=/work/aphernand001/hcc-singlecell
ENV="$PROJ/env"
cd "$PROJ" || exit 1
REPORT="$PROJ/cytotrace2_report.txt"
: > "$REPORT"
exec >> "$REPORT" 2>&1
section() { echo; echo "=================== $1 ==================="; date; }

export MAMBA_ROOT_PREFIX="$PROJ/mamba"
export TMPDIR="$PROJ/tmp"
mkdir -p "$TMPDIR"
MM="$PROJ/tools/bin/micromamba"

section "1. INSTALL"
"$MM" run -p "$ENV" Rscript -e '
options(Ncpus = 4, timeout = 900)
repos <- BiocManager::repositories()
lib <- .Library
check <- function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("%-12s %s %s\n", pkg, if (ok) "OK  " else "FAIL",
              if (ok) as.character(packageVersion(pkg)) else "-"))
  invisible(ok)
}
deps <- c("data.table", "doParallel", "dplyr", "HiClimR", "magrittr", "plyr", "Rfast", "RSpectra", "stringr")
missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
cat("Missing CRAN dependencies:", if (length(missing)) missing else "none", "\n")
if (length(missing)) install.packages(missing, lib = lib, repos = repos)
for (p in deps) check(p)

cat("\n--- CytoTRACE2 (GitHub) ---\n")
remotes::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r", lib = lib, repos = repos,
                        dependencies = c("Depends", "Imports", "LinkingTo"),
                        upgrade = "never", build_vignettes = FALSE)
check("CytoTRACE2")
'

section "2. FUNCTIONAL CHECK (synthetic matrix, 200 cells)"
"$MM" run -p "$ENV" Rscript -e '
suppressPackageStartupMessages(library(CytoTRACE2))
genes <- unique(read.table("data/ref/gene_order_hg38.txt", sep = "\t")$V1)[1:5000]
set.seed(1)
m <- matrix(rpois(length(genes) * 200, 0.5), nrow = length(genes),
            dimnames = list(genes, paste0("c", 1:200)))
res <- tryCatch(cytotrace2(m, species = "human", is_seurat = FALSE, ncores = 2),
                error = function(e) e)
if (inherits(res, "error")) cat("cytotrace2: FAIL -", conditionMessage(res), "\n") else {
  cat("cytotrace2: OK -", nrow(res), "cells; columns:", paste(colnames(res), collapse = ", "), "\n")
  print(summary(res$CytoTRACE2_Score))
}
'

section "3. UPDATE LOCK FILES"
"$MM" env export -p "$ENV" > "$PROJ/environment.lock.yml"
"$MM" run -p "$ENV" Rscript -e '
pk <- c("presto","CellChat","liana","msigdbrdata","NMF","leidenbase","CytoTRACE2")
rows <- lapply(pk, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(NULL)
  d <- packageDescription(p)
  data.frame(package = p, version = d$Version,
             repo = if (!is.null(d$RemoteRepo)) paste0(d$RemoteUsername, "/", d$RemoteRepo) else d$Repository,
             sha  = if (!is.null(d$RemoteSha)) d$RemoteSha else NA)
})
tab <- do.call(rbind, rows)
write.table(tab, "github_packages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)'

rm -rf "$TMPDIR"/*
du -sh "$ENV"
section "FIN"
