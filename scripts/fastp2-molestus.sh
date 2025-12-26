#!/bin/bash
#bash fastp2-molestus.sh ./molestus > fastp2-molestus.out 2> fastp2-molestus.err
# Directory containing the files
data_dir=${1} #4.enero

# Put bash in strict mode
set -e
set -o pipefail

# Output and log directories (adjust or create as needed)
output_dir="$data_dir/fastp"
#mkdir -p "$output_di

# Loop through all *_1.fastq files
for r1 in "$data_dir/to-map/primer-fastp"/*.R1.fastp.fq.gz; do
    # Get base name by removing _1.fastq
    sample=$(basename "$r1" .R1.fastp.fq.gz)
    r2="$data_dir/to-map/primer-fastp/${sample}.R2.fastp.fq.gz"

    # Make sure pair exists
    if [[ ! -f "$r2" ]]; then
        echo "Warning: missing pair for $sample"
        continue
    fi

    echo "Processing $sample"

    # Run fastp
    fastp \
        -i "$r1" -I "$r2" \
        -o "$data_dir/${sample}.R1.fastp2.fq.gz" -O "$data_dir/${sample}.R2.fastp2.fq.gz" \
        -h "$output_dir/${sample}.fastp2.html" -j "$output_dir/${sample}.fastp2.json" \
         --detect_adapter_for_pe \
        --adapter_fasta /sietch_colab/scamison/talapas/data/adapters-fastp2-molestus.fasta