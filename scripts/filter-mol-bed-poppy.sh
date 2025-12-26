#!/bin/bash

#set strict mode
set -euo pipefail

chrom=${1}  # Specify the chromosome to filter, e.g., "chr2"

# This script filters a BCF file using a BED file for SNPs on chromosome $chrom
echo "Filtering SNPs for chromosome ${chrom}..."
bcftools view \
 -R prueba-chr1.bed \
 molestus_mcbride_${chrom}_pop.bcf \
 -Ob -o prueba-bed-ch1.bcf \
 --threads 2
echo "Filtering complete for chromosome ${chrom}."

# index ouput BCF file
bcftools index molestus_mcbride_${chrom}_snps.bcf
echo "Indexing complete for molestus_mcbride_${chrom}_snps.bcf"