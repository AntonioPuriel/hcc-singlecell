#!/bin/bash
#SBATCH --job-name=hcc_fixpkg
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --output=logs/fixpkg_%j.out
#SBATCH --error=logs/fixpkg_%j.err

# Fixes the two packages that failed in 03_install_env.sh:
#   - CellChat: needs NMF >= 0.23 (conda-forge ships 0.21) -> NMF from CRAN source
#   - msigdbrdata: not on r-universe for R 4.5 -> install from GitHub
# Does not touch packages used by the QC/integration steps.
# Launch from the project root:  sbatch slurm/04_fix_packages.sh
# Report: fixpkg_report.txt

set +e
PROJ=/work/aphernand001/hcc-singlecell
ENV="$PROJ/env"
cd "$PROJ" || exit 1
REPORT="$PROJ/fixpkg_report.txt"
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
# Verify by loading, not by trusting install.packages() (which only warns on failure)
check <- function(pkg, min = NULL) {
  ok <- requireNamespace(pkg, quietly = TRUE) &&
        (is.null(min) || packageVersion(pkg) >= min)
  cat(sprintf("%-10s %s %s\n", pkg, if (ok) "OK  " else "FAIL",
              if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else "-"))
  invisible(ok)
}
cat("\n--- NMF (CRAN source, >= 0.23) ---\n")
install.packages("NMF", lib = lib, repos = repos, type = "source")
nmf_ok <- check("NMF", "0.23.0")

cat("\n--- CellChat (GitHub) ---\n")
if (nmf_ok) remotes::install_github("jinworks/CellChat", lib = lib, repos = repos,
                                    dependencies = c("Depends", "Imports", "LinkingTo"),
                                    upgrade = "never", build_vignettes = FALSE)
check("CellChat")

cat("\n--- msigdbrdata (GitHub) ---\n")
remotes::install_github("igordot/msigdbdf", lib = lib, repos = repos,
                        upgrade = "never", build_vignettes = FALSE)
check("msigdbrdata")
'

section "2. FUNCTIONAL CHECKS"
"$MM" run -p "$ENV" Rscript -e '
h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"), error = function(e) e)
if (inherits(h, "error")) cat("msigdbr Hallmark: FAIL -", conditionMessage(h), "\n") else
  cat("msigdbr Hallmark: OK -", length(unique(h$gs_name)), "gene sets\n")
cc <- tryCatch({ suppressMessages(library(CellChat)); nrow(CellChatDB.human$interaction) }, error = function(e) NA)
cat("CellChatDB.human interactions:", cc, "\n")
lr <- tryCatch(nrow(liana::select_resource("Consensus")[[1]]), error = function(e) NA)
cat("LIANA consensus resource interactions:", lr, "\n")
'

section "3. UPDATE LOCK FILES"
"$MM" env export -p "$ENV" > "$PROJ/environment.lock.yml"
"$MM" run -p "$ENV" Rscript -e '
pk <- c("presto","CellChat","liana","msigdbrdata","NMF","leidenbase")
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
