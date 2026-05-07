#!/bin/bash
set -euo pipefail

REFS="/sietch_colab/scamison/talapas/wolbachia/strains/wPip-pk1-refs.fasta"
CONS_DIR="/sietch_colab/scamison/talapas/wolbachia/consensus"
IDLIST="/sietch_colab/scamison/talapas/wolbachia/strains/ids-to-blast-Wb.txt"
OUTDIR="/sietch_colab/scamison/talapas/wolbachia/strains/3blast_results_pk1_definitivo"
SUMMARY="${OUTDIR}/pk1_hits_summary.tsv"

mkdir -p "$OUTDIR"

# header
echo -e "sample\tqseqid\tsseqid\tpident\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tmismatch\tgaps\tqcovhsp\tsstrand" > "$SUMMARY"

while read -r id; do
  fasta="${CONS_DIR}/${id}_Wb_consensus.fasta"
  if [[ -f "$fasta" ]]; then
    echo "[INFO] Procesando muestra: $id"

    makeblastdb -in "$fasta" -dbtype nucl -out "${OUTDIR}/${id}_db" >/dev/null 2>&1

    blastn -db "${OUTDIR}/${id}_db" -query "$REFS" \
      -evalue 1e-20 \
      -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore mismatch gaps qcovhsp sstrand' \
      > "${OUTDIR}/${id}.tsv"

    awk -v samp="$id" 'BEGIN{OFS="\t"} {print samp, $0}' "${OUTDIR}/${id}.tsv" >> "$SUMMARY"
  else
    echo "[WARN] No existe archivo ${fasta}" >&2
  fi
done < "$IDLIST"

echo "[INFO] Resumen combinado: $SUMMARY"
echo "[INFO] Total de hits: $(($(wc -l < "$SUMMARY") - 1))"