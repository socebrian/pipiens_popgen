#!/bin/bash

#this script runs qualimap over all *merged.bam files in a directory and the runs multibam qc over all the qualimap outputs
#module load miniconda3
#conda activate culex

dir=${1} #3.repetidas
export TMPDIR="./tmpr"
subdir=${2}

# Paths and variables
input_dir="./${dir}/${subdir}"
output_dir="./${dir}/qualimaps"

# Limit Java memory and GC threads for Qualimap
export QUALIMAP_OPTS="-Xmx8g -XX:ParallelGCThreads=3 -XX:ConcGCThreads=3"

## Loop through each merged BAM file in the input directory
for bam_file in "$input_dir"/*merged.sorted.bam; do
    # Check if the BAM file exists to avoid errors
    if [[ -f "$bam_file" ]]; then
        # Extract sample ID from filename
        sample_id=$(basename "$bam_file" .merged.sorted.bam)

        echo "Running qualimap for sample: ${sample_id}"
        qualimap bamqc -bam "$bam_file" -outdir "${output_dir}/${sample_id}" -outfile "${sample_id}.merged.bamqc.pdf" --java-mem-size=8G -outformat PDF:HTML -nt 6
    else
        echo "No BAM files found in ${input_dir}"
    fi
done

echo "done"
# sbatch --array=1-5 -n 1 --cpus-per-task=1 --mem-per-cpu=2G -t 02:00:00 --job-name=3rep.qualimap  qualimap.sh 3.repetidas
