#!/bin/bash
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL
#SBATCH --output=bwa_align_%A_%a.out
#SBATCH --error=bwa_align_%A_%a.err
#SBATCH --array=1-50


module load miniconda3
conda activate culex
module load bwa
module load samtools
module load java/17
module load picard

# Put bash in strict mode
set -e
set -o pipefail


#variables
dir=${1} #3.repetidas
list=${2} #3.list-fastp-output.txt  (print el 3 campo)
ref=${3} # idCulPipi1.1.primary.fa

input="/path/to/${dir}/data/${list}"

IFS=$'\t' read -r -a samples <<< $(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$input")
# Extract parameters
id=${samples[1]}
lane=${samples[2]}
r=${samples[3]}

###############################################
# ALLIGN SAMPLES
###############################################
#create output dir
#echo "Creating output directory: /path/to/${dir}/bams"
#mkdir /path/to/${dir}/bams


#paired alignment
bwa mem \
    /path/to/data/${ref}\
    /path/to/${dir}/data/${id}.${lane}.R1.fastp.fq.gz \
    /path/to/${dir}/data/${id}.${lane}.R2.fastp.fq.gz \
    -t 15|
samtools view -hbS - -o /path/to/${dir}/bams/${id}.${lane}.PR.bam -@ 15


###############################################
# INDEX SAMPLES
###############################################


# Debug information
echo "Sample List: $output_path"
echo "Batch: $dir"
echo "ID: $id"
echo "Lane: $lane"
echo "R: $r"
echo "Executing: sbatch index-bam.sh $id $lane $r $dir"

# File paths
input_bam="/path/to/${dir}/bams/${id}.${lane}.${r}.bam"
sorted_bam="/path/to/${dir}/bams/${id}.${lane}.${r}.sorted.bam"
sorted_rg_bam="/path/to/${dir}/bams/${id}.${lane}.${r}.sorted.rg.bam"

# Check if input BAM file exists
if [[ ! -f "$input_bam" ]]; then
    echo "Error: Input BAM file not found at $input_bam"
    exit 1
fi

# Sort BAM file
echo "Sorting BAM file: $input_bam"
samtools sort "$input_bam" -o "$sorted_bam"

# Add or replace read group
#id flag in picard should be only the mosquito ID so different lanes can be merged afterwards wihout problems messing with the metadata
echo "Adding or replacing read groups in BAM file: $sorted_bam"
picard AddOrReplaceReadGroups \
    I="$sorted_bam" \
    O="$sorted_rg_bam" \
    RGID="${id}" \
    RGLB="${dir}" \
    RGPL=Illumina \
    RGPU="${id}.${lane}.${r}" \
    RGSM="${id}" \
    VALIDATION_STRINGENCY=SILENT

# Index the sorted BAM file with read groups
echo "Indexing BAM file: $sorted_rg_bam"
samtools index "$sorted_rg_bam"

###############################################
# GENERATE ALLIGNMENT STATS
###############################################

echo "Generating alignment stats for: $sorted_rg_bam, output: ${id}.${lane}.stats"
samtools flagstat "${sorted_rg_bam}" > "${id}.${lane}.stats"


#sbatch --array=1-168 -n 1 --cpus-per-task=15 --mem-per-cpu=20G -t 12:00:00 allign-bams-index.sh 4.enero list-fastp-output.txt idCulPipi1.1.primary.fa