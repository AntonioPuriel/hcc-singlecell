#!/bin/bash
#SBATCH --job-name=hcc_download
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err

# Descarga GSE149614 (Lu et al., 2022) desde GEO: archivos suplementarios + series matrix.
# Lanzar desde /work/aphernand001/hcc-singlecell:  sbatch slurm/02_download.sh
# Resultado: data/raw/ + download_report.txt

set +e
PROJ=/work/aphernand001/hcc-singlecell
RAW="$PROJ/data/raw"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614"
mkdir -p "$RAW"
cd "$PROJ" || exit 1

REPORT="$PROJ/download_report.txt"
: > "$REPORT"
exec >> "$REPORT" 2>&1
section() { echo; echo "=================== $1 ==================="; date; }

section "LISTADO REMOTO"
for d in suppl matrix; do
  echo "--- $d ---"
  curl -s "$BASE/$d/" | grep -o 'href="[^"]*"' | sed 's/href="//; s/"$//' | grep -v -E '^(\?|/|\.\./)'
done

section "DESCARGA"
for d in suppl matrix; do
  wget -q -r -np -nd -nH -N -R "index.html*" -e robots=off -P "$RAW" "$BASE/$d/"
  echo "$d: exit=$?"
done

section "ARCHIVOS"
ls -lh "$RAW"

section "INTEGRIDAD (gzip -t) Y CHECKSUMS"
for f in "$RAW"/*.gz; do
  gzip -t "$f" && echo "OK   $(basename "$f")" || echo "FAIL $(basename "$f")"
done
(cd "$RAW" && md5sum * > "$PROJ/data/raw_md5.txt")
cat "$PROJ/data/raw_md5.txt"

section "ESPACIO"
du -sh "$RAW"
section "FIN"
