#!/bin/bash

THR=8
NPR=4
work=/path/to/workdir
invcf=$work/vcf
outdir=$work/pixy_outs/
sites=$work/prefix

# The statistics to run
stats=(
    pi
    dxy
    fst
)

# Window size for the calculations (in bp)
windows=(
    5000
    10000
)

# The popmaps to run
popmaps=(
    ecot95
    ecot60
)

# The chromosomes to run
chrs=(
    1
    2
    3
)

# Loop over the three chromosomes
for chr in "${chrs[@]}"; do
    # Get the sites for that chromosome
    chr_sites=$work/chr${chr}_sites.tsv
    cat $sites | awk -v c="$chr" '$1 == c {print $0}' > $chr_sites
    # Loop over the statistics
    for stat in "${stats[@]}"; do
        # Loop over the window sizes
        for win in "${windows[@]}"; do
            # Loop over the population assignments
            for pop in "${popmaps[@]}"; do
                # Select the popmap to run
                popmap=$work/popmap_${pop}.tsv
                # Generate output prefix
                prefix=pipiens1234_N90_${pop}_chr${chr}_win${win}bp
                # Pixy command
                cmd=(
                    pixy
                    --stats $stat
                    --vcf $invcf
                    --populations $popmap
                    --window_size $win
                    --chromosomes "'${chr}'"
                    --n_cores $THR
                    # --sites_file $chr_sites
                    --output_folder $outdir
                    --output_prefix $prefix
                )
                echo "${cmd[@]}"
            done
        done
    done
done | parallel -j $NPR