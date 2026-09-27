#!/bin/bash
#SBATCH --job-name=hcc_install
#SBATCH --partition=standard
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=logs/install_%j.out
#SBATCH --error=logs/install_%j.err

# Instala el entorno definitivo de hcc-singlecell en /work/aphernand001/hcc-singlecell/env
# Requiere environment.yml en el mismo directorio.
# Lanzar desde /work/aphernand001/hcc-singlecell:  sbatch slurm/03_install_env.sh
# Resultado: install_report.txt, environment.lock.yml, github_packages.tsv

set +e
PROJ=/work/aphernand001/hcc-singlecell
ENV="$PROJ/env"
cd "$PROJ" || exit 1

REPORT="$PROJ/install_report.txt"
: > "$REPORT"
exec >> "$REPORT" 2>&1
section() { echo; echo "=================== $1 ==================="; date; }

export MAMBA_ROOT_PREFIX="$PROJ/mamba"
export CONDA_PKGS_DIRS="$PROJ/conda_pkgs"
export TMPDIR="$PROJ/tmp"
export R_LIBS_USER="$ENV/lib/R/library"   # nada en el home
MM="$PROJ/tools/bin/micromamba"
mkdir -p "$TMPDIR"

# ------------------------------------------------------------------ 0. LIMPIEZA
section "0. LIMPIEZA DE LAS PRUEBAS"
# Ruta Singularity descartada (conda resuelve) y matriz normalizada de GEO no usada
rm -rf "$PROJ/sing_cache" "$PROJ/containers" "$PROJ/tmp/rlib"
rm -f  "$PROJ/data/raw/GSE149614_HCC.scRNAseq.S71915.normalized.txt.gz"
du -sh "$PROJ"

# ------------------------------------------------------------------ 1. ENTORNO BASE
section "1. CREAR ENTORNO DESDE environment.yml"
"$MM" create -y --no-rc -p "$ENV" -f "$PROJ/environment.yml"
BASE_EXIT=$?
echo "exit=$BASE_EXIT"
if [ "$BASE_EXIT" -ne 0 ]; then echo "ERROR: el entorno base no se creó. Abortando."; exit 1; fi

# ------------------------------------------------------------------ 2. EXTRAS CONDA
section "2. EXTRAS DE CONDA (uno a uno, no bloqueantes)"
# Dependencias de CellChat y LIANA disponibles en conda
EXTRAS="bioconductor-singlecellexperiment bioconductor-scuttle bioconductor-scran \
bioconductor-omnipathr bioconductor-complexheatmap bioconductor-biocneighbors \
r-circlize r-nmf r-ggalluvial r-rspectra r-svglite r-reticulate r-sna r-ggnetwork \
r-future r-leidenbase"
for p in $EXTRAS; do
  "$MM" install -y --no-rc -p "$ENV" -c conda-forge -c bioconda "$p" > "$TMPDIR/extra_$p.log" 2>&1 \
    && echo "OK    $p" || { echo "FAIL  $p"; tail -5 "$TMPDIR/extra_$p.log" | sed 's/^/      /'; }
done

# ------------------------------------------------------------------ 3. PAQUETES R FUERA DE CONDA
section "3. PAQUETES DE GITHUB / R-UNIVERSE (no bloqueantes)"
"$MM" run -p "$ENV" Rscript -e '
options(Ncpus = 8, timeout = 600)
repos <- BiocManager::repositories()
lib <- .Library
try_install <- function(label, expr) {
  res <- tryCatch({ expr; "OK" }, error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(sprintf("%-10s %s\n", label, res))
}
if (!requireNamespace("leidenbase", quietly = TRUE))
  try_install("leidenbase", install.packages("leidenbase", lib = lib, repos = repos))
try_install("msigdbdf", install.packages("msigdbdf", lib = lib,
            repos = c("https://igordot.r-universe.dev", repos)))
gh <- c(presto = "immunogenomics/presto",
        CellChat = "jinworks/CellChat",
        liana = "saezlab/liana")
for (n in names(gh))
  try_install(n, remotes::install_github(gh[[n]], lib = lib, repos = repos,
              dependencies = c("Depends", "Imports", "LinkingTo"),
              upgrade = "never", build_vignettes = FALSE))
'

# ------------------------------------------------------------------ 4. VERIFICACIÓN
section "4. VERIFICACIÓN: ¿CARGA CADA PAQUETE?"
"$MM" run -p "$ENV" Rscript -e '
pk <- c("Seurat","SeuratObject","harmony","scDblFinder","SingleR","celldex","infercnv",
        "DESeq2","fgsea","msigdbr","msigdbdf","slingshot","monocle3","igraph","leidenbase",
        "liana","CellChat","presto","OmnipathR","targets","tarchetypes","crew","qs2",
        "data.table","dplyr","ggplot2","patchwork","rmarkdown","quarto")
ok <- vapply(pk, function(p) requireNamespace(p, quietly = TRUE), logical(1))
ver <- vapply(pk, function(p) if (ok[[p]]) as.character(packageVersion(p)) else "-", character(1))
print(data.frame(paquete = pk, carga = ok, version = ver), row.names = FALSE)
cat("\n", R.version.string, "\n")
cat("Quarto CLI:", Sys.which("quarto"), "\n")
'

# ------------------------------------------------------------------ 5. LOCKS
section "5. EXPORTAR VERSIONES EXACTAS"
"$MM" env export -p "$ENV" > "$PROJ/environment.lock.yml"
echo "environment.lock.yml: $(wc -l < "$PROJ/environment.lock.yml") líneas"
"$MM" run -p "$ENV" Rscript -e '
pk <- c("presto","CellChat","liana","msigdbdf","leidenbase")
rows <- lapply(pk, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(NULL)
  d <- packageDescription(p)
  data.frame(package = p, version = d$Version,
             repo = if (!is.null(d$RemoteRepo)) paste0(d$RemoteUsername, "/", d$RemoteRepo) else d$Repository,
             sha  = if (!is.null(d$RemoteSha)) d$RemoteSha else NA)
})
tab <- do.call(rbind, rows)
write.table(tab, "github_packages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)
'

# ------------------------------------------------------------------ 6. LIMPIEZA FINAL
section "6. LIMPIAR CACHÉ DE PAQUETES Y ESPACIO FINAL"
"$MM" clean -a -y > /dev/null 2>&1
rm -rf "$TMPDIR"/*
du -sh "$ENV" "$PROJ/conda_pkgs" "$PROJ/data" 2>/dev/null
du -sh "$PROJ"

section "FIN"
