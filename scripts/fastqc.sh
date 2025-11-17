#!/bin/bash
#SBATCH --mem-per-cpu=20G
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err

#this script will list all the files in a directory and write the output to a file
#specificly designed to list the files to perform fastp
#loasd modules
module load miniconda3
conda activate fastp

# Directory containing the files
dir=${1} #4.enero

# Put bash in strict mode
set -e
set -o pipefail

#run fastqc and multiqc
for file in /path/to/${dir}/data/*; do
    echo $file
    fastqc ${file}  -o "/path/to/${dir}/fastqc"
done

echo "fastqc done, stariting multiqc in ${dir}/fastqc"
multiqc "${dir}/fastqc" --outdir "${dir}/fastqc"

#sbatch -n 1 --cpus-per-task=10 --mem-per-cpu=20G -t 12:00:00 fastqc.sh 4.enero