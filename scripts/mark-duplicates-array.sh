#!/bin/bash
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL
#SBATCH --output=mark_dup_%A_%a.out
#SBATCH --error=mark_dup_%A_%a.err
#SBATCH --array=1-50

#this script takes a list of bam files and remove duplicates

# Load any necessary modules
module load miniconda3
conda activate vcf
module load racs-eb/1
module load samtools
module load picard

#variables
dir=${1} #3.repetidas
list=${2} #list of bam files

# Specify the path to the config file
input="/path/to/${dir}/data/${list}"
bam_dir="/path/to/${dir}/bams"
out_dir="/path/to/${dir}/mark-duplicates"

IFS=$'\t' read -r -a samples <<< $(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$input")
# Extract parameters
id=${samples[1]}
lane=${samples[2]}
r=${samples[3]}

# Read the sample ID for the current array task
sample=$(awk -v taskid="$SLURM_ARRAY_TASK_ID" '$1 == taskid {print $2}' "$sample_list")
in_bam="$bam_dir/${id}.${lane}.PR.sorted.rg.bam"
# Define output filenames
output_bam="$out_dir/${id}.${lane}.md.bam"
output_metrics="$out_dir/${id}.${lane}_marked_dup_metrics.txt"


# Debug: Print the sample variable
echo "Sample: ${is}.${lane}"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=$in_bam \
    O=$output_bam \
    M=$output_metrics \
    REMOVE_DUPLICATES=true

##sbatch --array=1-5 -n 1 --cpus-per-task=10 --mem-per-cpu=20G -t 12:00:00 --job-name=3rep.markdup  mark-duplicates-array.sh 3.repetidas list-mark-duplicates.txt