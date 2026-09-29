#!/bin/bash
#SBATCH --job-name=hcc_refs
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/refs_%j.out
#SBATCH --error=logs/refs_%j.err

# Gene positions (hg38) needed by inferCNV: gene<TAB>chr<TAB>start<TAB>end
# 1st choice: the inferCNV/CTAT gene-order file from the Broad.
# Fallback: built from the NCBI RefSeq GRCh38.p14 annotation (NCBI is reachable
# from the compute nodes; gene_id in that GTF is the HGNC symbol, as in the matrix).
# Launch from the project root:  sbatch slurm/05_download_refs.sh
# Report: refs_report.txt

set +e
PROJ=/work/aphernand001/hcc-singlecell
REF="$PROJ/data/ref"; OUT="$REF/gene_order_hg38.txt"
mkdir -p "$REF" "$PROJ/logs"; cd "$PROJ" || exit 1
REPORT="$PROJ/refs_report.txt"; : > "$REPORT"; exec >> "$REPORT" 2>&1
section() { echo; echo "=================== $1 ==================="; date; }

section "1. CTAT gene order (Broad)"
curl -sSfL -m 300 -o "$REF/ctat.txt" https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt
if [ $? -eq 0 ] && [ -s "$REF/ctat.txt" ]; then
  mv "$REF/ctat.txt" "$OUT"; echo "source: CTAT hg38_gencode_v27" > "$REF/gene_order_source.txt"
else
  echo "CTAT not reachable -> NCBI fallback"; rm -f "$REF/ctat.txt"
  section "2. NCBI RefSeq GTF fallback"
  GTF=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
  curl -sSfL -m 1800 "$GTF" | gzip -dc | awk -F'\t' '
    $3 == "gene" && $1 ~ /^NC_0000(0[1-9]|1[0-9]|2[0-4])\./ {
      n = substr($1, 8, 2) + 0
      chr = (n == 23) ? "chrX" : (n == 24) ? "chrY" : "chr" n
      match($9, /gene_id "[^"]+"/); g = substr($9, RSTART + 9, RLENGTH - 10)
      if (!(g in seen)) { seen[g] = 1; key = (n < 10 ? "0" n : n); print key "\t" $4 "\t" g "\t" chr "\t" $4 "\t" $5 }
    }' | sort -k1,1 -k2,2n | cut -f3- > "$OUT"
  echo "source: NCBI RefSeq GRCh38.p14 (GCF_000001405.40)" > "$REF/gene_order_source.txt"
fi

section "3. CHECK"
cat "$REF/gene_order_source.txt"
wc -l "$OUT"; head -3 "$OUT"; cut -f2 "$OUT" | sort | uniq -c | sort -k2,2V | head -30
md5sum "$OUT"
section "FIN"
