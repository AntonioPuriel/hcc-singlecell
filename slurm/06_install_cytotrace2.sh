#!/bin/bash
#SBATCH --job-name=hcc_cytotrace2
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/cytotrace2_%j.out
#SBATCH --error=logs/cytotrace2_%j.err

# Installs CytoTRACE2 (GitHub, R version) into env/ and checks it on a synthetic matrix.
#
# CytoTRACE2 imports HiClimR, and HiClimR imports ncdf4. First attempt (2 Oct) failed:
#  - ncdf4 from CRAN source compiled against the system netCDF (CentOS 7), which is too
#    old (NC_FORMAT_CDF5 undeclared);
#  - r-ncdf4 from conda-forge cannot be solved against the hdf5 already in env/.
# Route A: install only the netCDF C library (libnetcdf) from conda-forge with every
#          existing package frozen, then build ncdf4 from source against env/bin/nc-config.
# Route B (only if A fails): install HiClimR with ncdf4 removed from its imports.
#          ncdf4 is used by HiClimR only to read/write netCDF climate files; CytoTRACE2
#          uses HiClimR for fast correlations, so its results are unaffected.
# The report states which route was used.
#
# Exit status 0 only if CytoTRACE2 loads and runs, so the pipeline can be chained:
#   jid=$(sbatch --parsable slurm/06_install_cytotrace2.sh)
#   sbatch --dependency=afterok:$jid slurm/10_run_pipeline.sh
# Report: cytotrace2_report.txt

set +e
PROJ=/work/aphernand001/hcc-singlecell
ENV="$PROJ/env"
cd "$PROJ" || exit 1
mkdir -p logs
REPORT="$PROJ/cytotrace2_report.txt"
: > "$REPORT"
exec >> "$REPORT" 2>&1
section() { echo; echo "=================== $1 ==================="; date; }

export MAMBA_ROOT_PREFIX="$PROJ/mamba"
export TMPDIR="$PROJ/tmp"
mkdir -p "$TMPDIR"
MM="$PROJ/tools/bin/micromamba"

section "0. hdf5 / netCDF CURRENTLY IN env/"
"$MM" list -p "$ENV" | grep -E "hdf5|netcdf" || echo "(none)"
echo "system nc-config: $(command -v nc-config) $(nc-config --version 2>/dev/null)"

section "1. ROUTE A: libnetcdf FROM CONDA-FORGE (existing packages frozen)"
"$MM" list -p "$ENV" > "$TMPDIR/pkgs_before.txt"
"$MM" install -y -p "$ENV" -c conda-forge --freeze-installed libnetcdf
"$MM" list -p "$ENV" > "$TMPDIR/pkgs_after.txt"
echo "--- packages added or changed in env/:"
diff "$TMPDIR/pkgs_before.txt" "$TMPDIR/pkgs_after.txt" | grep '^[<>]' || echo "(none)"
NC_CONFIG="$ENV/bin/nc-config"
if [ -x "$NC_CONFIG" ]; then echo "env nc-config: $("$NC_CONFIG" --version)"; else echo "no nc-config in env/"; fi
export NC_CONFIG

section "2. ncdf4 / HiClimR / OTHER CRAN DEPENDENCIES"
"$MM" run -p "$ENV" Rscript -e '
options(Ncpus = 4, timeout = 900)
repos <- BiocManager::repositories()
lib <- .Library
ok <- function(p) requireNamespace(p, quietly = TRUE)

# Route A: ncdf4 against the netCDF library in env/
nc <- Sys.getenv("NC_CONFIG")
if (!ok("ncdf4") && file.exists(nc)) {
  install.packages("ncdf4", lib = lib, repos = repos, type = "source",
                   configure.args = paste0("--with-nc-config=", nc))
}
route <- NA
if (ok("ncdf4")) {
  route <- "A (ncdf4 built against conda libnetcdf)"
  if (!ok("HiClimR")) install.packages("HiClimR", lib = lib, repos = repos)
}

# Route B: HiClimR without ncdf4
if (!ok("HiClimR")) {
  route <- "B (HiClimR installed without ncdf4)"
  dir <- file.path(tempdir(), "hiclimr"); dir.create(dir, showWarnings = FALSE)
  tgz <- download.packages("HiClimR", destdir = dir, repos = repos, type = "source")[1, 2]
  untar(tgz, exdir = dir)
  pkg <- file.path(dir, "HiClimR")
  d <- read.dcf(file.path(pkg, "DESCRIPTION"))
  for (f in intersect(c("Depends", "Imports"), colnames(d))) {
    x <- trimws(strsplit(d[1, f], ",")[[1]])
    d[1, f] <- paste(x[!grepl("^ncdf4", x)], collapse = ", ")
  }
  write.dcf(d, file.path(pkg, "DESCRIPTION"))
  ns <- readLines(file.path(pkg, "NAMESPACE"))
  writeLines(ns[!grepl("ncdf4", ns)], file.path(pkg, "NAMESPACE"))
  install.packages(pkg, lib = lib, repos = NULL, type = "source")
}
cat("\nRoute used:", if (ok("HiClimR")) route else "none (HiClimR failed)", "\n\n")

deps <- c("data.table", "doParallel", "dplyr", "magrittr", "plyr", "Rfast", "RSpectra", "stringr")
missing <- deps[!vapply(deps, ok, logical(1))]
if (length(missing)) install.packages(missing, lib = lib, repos = repos)
for (p in c("ncdf4", "HiClimR", deps))
  cat(sprintf("%-12s %s\n", p, if (ok(p)) as.character(packageVersion(p)) else "MISSING"))

cat("\n--- CytoTRACE2 (GitHub) ---\n")
remotes::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r", lib = lib, repos = repos,
                        dependencies = c("Depends", "Imports", "LinkingTo"),
                        upgrade = "never", build_vignettes = FALSE)
cat("CytoTRACE2:", if (ok("CytoTRACE2")) as.character(packageVersion("CytoTRACE2")) else "MISSING", "\n")
'

section "3. FUNCTIONAL CHECK (synthetic matrix, 200 cells)"
"$MM" run -p "$ENV" Rscript -e '
if (!requireNamespace("CytoTRACE2", quietly = TRUE)) { cat("cytotrace2: FAIL - not installed\n"); quit(status = 1) }
suppressPackageStartupMessages(library(CytoTRACE2))
genes <- unique(read.table("data/ref/gene_order_hg38.txt", sep = "\t")$V1)[1:5000]
set.seed(1)
m <- matrix(rpois(length(genes) * 200, 0.5), nrow = length(genes),
            dimnames = list(genes, paste0("c", 1:200)))
res <- tryCatch(cytotrace2(m, species = "human", is_seurat = FALSE, ncores = 2),
                error = function(e) e)
if (inherits(res, "error")) { cat("cytotrace2: FAIL -", conditionMessage(res), "\n"); quit(status = 1) }
cat("cytotrace2: OK -", nrow(res), "cells; columns:", paste(colnames(res), collapse = ", "), "\n")
print(summary(res$CytoTRACE2_Score))
'
STATUS=$?
echo "functional check exit=$STATUS"

section "4. UPDATE LOCK FILES"
"$MM" env export -p "$ENV" > "$PROJ/environment.lock.yml"
"$MM" run -p "$ENV" Rscript -e '
pk <- c("presto","CellChat","liana","msigdbrdata","NMF","leidenbase","HiClimR","CytoTRACE2")
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
exit $STATUS
