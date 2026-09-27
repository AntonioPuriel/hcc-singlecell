#!/bin/bash
#SBATCH --job-name=hcc_envcheck
#SBATCH --partition=short
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/envcheck_%j.out
#SBATCH --error=logs/envcheck_%j.err

# Diagnóstico del entorno de Pyrene para el proyecto hcc-singlecell.
# Lanzar desde /work/aphernand001/hcc-singlecell:  sbatch slurm/00_check_env.sh
# Resultado: envcheck_report.txt (+ modules_full.txt) en el mismo directorio.

set +e
cd /work/aphernand001/hcc-singlecell || exit 1
REPORT="envcheck_report.txt"
: > "$REPORT"
exec >> "$REPORT" 2>&1

section() { echo; echo "=================== $1 ==================="; }

section "SISTEMA"
date
hostname
grep -E '^(NAME|VERSION)=' /etc/os-release
echo "CPUs visibles: $(nproc)"
free -g

section "MÓDULOS (lista completa en modules_full.txt)"
module avail > modules_full.txt 2>&1
module avail 2>&1 | tr ' ' '\n' | grep -i -E '^(R|r-|rstudio|seurat|singularity|apptainer|python|miniconda|anaconda|mamba|gcc|hdf5|jags|gsl|glpk|cmake|quarto|pandoc)([/-]|$)' | sort -u

section "INTENTO: module load R (versión por defecto)"
if module load R 2>/dev/null; then
  which Rscript
  Rscript -e 'cat(R.version.string, "\n"); pk <- c("Seurat","SeuratObject","harmony","scDblFinder","SingleR","celldex","DESeq2","fgsea","msigdbr","slingshot","monocle3","infercnv","copykat","liana","CellChat","targets","crew","crew.cluster","qs","igraph","renv","quarto"); inst <- rownames(installed.packages()); print(data.frame(paquete=pk, instalado=pk %in% inst))'
  module unload R
else
  echo "No hay módulo 'R' por defecto (revisar modules_full.txt por nombres alternativos)."
fi

section "CONTENEDORES"
for c in singularity apptainer; do
  if command -v "$c" >/dev/null 2>&1; then echo "$c: $($c --version)"; else echo "$c: no encontrado en PATH"; fi
done

section "CONDA COMPARTIDO"
module load bioconda-tools/3
conda --version
command -v mamba && mamba --version | head -1
echo "--- entornos compartidos con R/single-cell en el nombre ---"
ls /softs/contrib/apps/anaconda/3/envs 2>/dev/null | grep -i -E '(^r|r-|seurat|single|scrna|bioc)'
echo "--- conda env list ---"
conda env list
echo "--- canales configurados ---"
conda config --show channels
module unload bioconda-tools/3

section "ACCESO A INTERNET DESDE NODO DE CÁLCULO"
for u in https://ftp.ncbi.nlm.nih.gov/geo/ https://cloud.r-project.org/ https://bioconductor.org/ https://github.com/ https://conda.anaconda.org/ https://cf.10xgenomics.com/; do
  code=$(curl -s -o /dev/null -m 20 -w '%{http_code}' "$u")
  echo "$u -> HTTP $code"
done

section "ESPACIO EN DISCO"
df -h /work /home 2>/dev/null
echo "--- home ---"
du -sh "$HOME" 2>/dev/null
echo "--- /work/aphernand001 por subdirectorio (timeout 20 min) ---"
timeout 1200 du -sh /work/aphernand001/* 2>/dev/null | sort -h
echo "--- total ---"
timeout 600 du -sh /work/aphernand001 2>/dev/null

section "SLURM"
sinfo -o "%P %a %l %c %m %D" 2>&1
echo "--- QOS / asociaciones ---"
sacctmgr -nP show assoc user="$USER" format=Account,Partition,QOS,GrpTRES,MaxTRES 2>&1
echo "--- ¿se puede lanzar sbatch desde un job? ---"
command -v sbatch
echo "--- mis jobs en cola ---"
squeue -u "$USER" -o "%i %P %j %T %M %l %C %m" 2>&1

section "FIN"
date
