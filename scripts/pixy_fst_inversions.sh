#!/bin/bash

THR=8
NPR=4
work=/sietch_colab/data_share/culex/divergence
invcf=$work/merged/pipiens1234_all.filtered.vcf.gz
#invcf=$work/merged/pipiens1234_all.filtered.maf0.03.vcf.gz
outdir=$work/pixy_outs

# The three inversions to process
invs=(
    inv1
    # inv3.1
    # inv3.2
)

# The statistics to run
stats=(
    pi
    dxy
    fst
)

# Window size for the calculations (in bp)
win=1

# Loop over the three inversions
for inv in "${invs[@]}"; do
    # Select the correct popmap
    popmap=$work/info/popmap_${inv}.tsv
    # Loop over the three chromosomes
    # for chr in {1..3}; do
    for chr in {2..2}; do
        # Loop over the statistics
        for stat in "${stats[@]}"; do
            # Generate output prefix
            prefix=pipiens1234_all_${inv}_chr${chr}_win${win}bp
            # Pixy command
            cmd=(
                pixy
                --stats $stat
                --vcf $invcf
                --populations $popmap
                --window_size $win
                --chromosomes "'${chr}'"
                --n_cores $THR
                --output_folder $outdir
                --output_prefix $prefix
                #--include_multiallelic_snps
            )
            echo "${cmd[@]}"
        done
    done
done | parallel -j $NPR