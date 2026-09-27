#!/bin/bash
#SBATCH --job-name=hcc_test_conda
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/envtest_conda_%j.out
#SBATCH --error=logs/envtest_conda_%j.err

# Ruta A: micromamba propio + conda-forge/bioconda. Solo dry-run, no instala.
# Lanzar desde /work/aphernand001/hcc-singlecell:  sbatch slurm/01_test_conda.sh
# Resultado: envtest_conda_report.txt (+ envtest_dryrun_full.txt)

set +e
PROJ=/work/aphernand001/hcc-singlecell
cd "$PROJ" || exit 1
mkdir -p tools conda_pkgs mamba tmp

REPORT="$PROJ/envtest_conda_report.txt"
: > "$REPORT"
exec >> "$REPORT" 2>&1
section() { echo; echo "=================== $1 ==================="; date; }

export MAMBA_ROOT_PREFIX="$PROJ/mamba"
export CONDA_PKGS_DIRS="$PROJ/conda_pkgs"
export TMPDIR="$PROJ/tmp"

section "0. SISTEMA"
hostname
ldd --version | head -1

section "A1. MICROMAMBA"
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$PROJ/tools" bin/micromamba
MM="$PROJ/tools/bin/micromamba"
"$MM" --version
"$MM" info

CH="--override-channels -c conda-forge -c bioconda"
PKGS="r-base r-seurat r-harmony bioconductor-scdblfinder bioconductor-singler bioconductor-celldex \
bioconductor-deseq2 bioconductor-fgsea r-msigdbr bioconductor-slingshot bioconductor-infercnv \
r-targets r-tarchetypes r-crew r-qs2 r-igraph r-remotes r-rmarkdown r-quarto"

section "A2. DRY-RUN DEL ENTORNO COMPLETO"
"$MM" create -n hcc-test --dry-run -y $CH $PKGS > "$PROJ/envtest_dryrun_full.txt" 2>&1
FULL_EXIT=$?
echo "exit=$FULL_EXIT"
echo "--- versiones elegidas ---"
grep -E '^\s*\+\s+(r-base|r-seurat|r-seuratobject|r-harmony|bioconductor-deseq2|bioconductor-infercnv|bioconductor-singler|r-targets|r-crew|r-qs2)\s' "$PROJ/envtest_dryrun_full.txt"
echo "--- últimas 40 líneas ---"
tail -40 "$PROJ/envtest_dryrun_full.txt"

if [ "$FULL_EXIT" -ne 0 ]; then
  section "A3. PAQUETE A PAQUETE"
  for p in $PKGS; do
    [ "$p" = "r-base" ] && continue
    out=$("$MM" create -n t --dry-run -y $CH r-base "$p" 2>&1)
    if [ $? -eq 0 ]; then
      echo "OK    $p  $(echo "$out" | grep -E "^\s*\+\s+$p\s" | awk '{print $3}')"
    else
      echo "FAIL  $p"; echo "$out" | tail -8 | sed 's/^/      /'
    fi
  done
fi

section "A4. monocle3 (opcional)"
"$MM" create -n t --dry-run -y $CH r-base r-monocle3 > /dev/null 2>&1 && echo "r-monocle3: OK" || echo "r-monocle3: FAIL"

section "ESPACIO"
du -sh "$PROJ"/tools "$PROJ"/mamba "$PROJ"/conda_pkgs 2>/dev/null
section "FIN"
