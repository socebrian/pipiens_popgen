#!/bin/bash
set -e
set -o pipefail
THR=4

# General paths
work=/sietch_colab/data_share/culex/divergence/molpip_div
out_dir=$work/phylogenies

# Input VCF
invcf=$work/pipiens1234_all.N90.miss0.0.maf0.05.vcf.gz

# Output
# For 2p
target="2p"
chrom=2
start=1
stop=112539609

# For 2q
target="2q"
chrom=2
start=113453184
stop=213114244

# Window size in KB, for prunning
window=1

# Allele frequency, for filtering (both min and max)
maf=0.05
xaf=$(awk "BEGIN {print 1-$maf}")

# Output VCF
outvcf=$out_dir/${target}-pipiens1234_all.N90.miss0.0.maf${maf}.${window}kb_prunned.vcf.gz

# Run BCFtools to:
# 1. Subset the variants to only the target region
# 2. Keep only biallelic SNPs
# 3. Filter sites based on maf (for NREF too rare or too common)
# 3. Prune based on a window size
bcftools view \
        --output-type u \
        --threads $THR \
        --regions "${chrom}:${start}-${stop}" \
        $invcf | \
    bcftools view \
        --output-type u \
        --type 'snps' \
        --min-alleles 2 \
        --max-alleles 2 | \
    bcftools view \
        --output-type u \
        --min-af $maf \
        --max-af $xaf | \
    bcftools +prune \
        --output-type v \
        --nsites-per-win 1 \
        --window ${window}kb | \
    bgzip -c > $outvcf
# Index the VCF
tabix $outvcf
