#!/bin/bash

NPR=4
src=/path/to/files
work=$src/path/to/outputs

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
    # Loop over the window sizes
    for win in "${windows[@]}"; do
        # Loop over the population assignments
        for pop in "${popmaps[@]}"; do

            # Generate the prefix for the target files
            prefix=pipiens1234_N90_${pop}_chr${chr}_win${win}bp

            # Input pixy files
            in_pi=$work/${prefix}_pi.txt
            in_dxy=$work/${prefix}_dxy.txt

            # command
            cmd=(
                python3
                $src/calculate_Da_pixy.py
                --pi $in_pi
                --dxy $in_dxy
                --out-dir $work
                --out-prefix $prefix
            )
            echo "${cmd[@]}"
        done
    done
done | parallel -j $NPR