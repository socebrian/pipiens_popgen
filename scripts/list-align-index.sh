#!/bin/bash
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=list_align_%j.out
#SBATCH --error=list_align_%j.err


module load miniconda3
conda activate culex

# Put bash in strict mode
set -e
set -o pipefail


#variables
dir=${1} #3.repetidas

###############################################
# CREATE LIST of SAMPLES
###############################################
echo "Creating list of samples to align and index: /path/to/${dir}/data/sample-list.txt"
output_path="/path/to/${dir}/data/sample-list.txt"
touch $output_path

# Write the header line to the output file
echo -e "taskid\tid\tlane\tread" > $output_path

# Initialize a counter for the ID column
taskid=1

# List all files in the directory
for filename in "/path/to/${dir}/data"/*.R1.fastp.fq.gz; do
    # Check if the file ends with any of the specified suffixes
    if [[ "$filename" =~ \.(R1.fastp.fq.gz)$ ]]; then
        # Extract the base name of the file (without the directory path)
        base_filename=$(basename "$filename")
        # Split the filename into components using the dot (.) as a delimiter
        IFS='.' read -r id lane read ext <<< "$base_filename"
        # Write the extracted parts to the output file
        echo -e "$taskid\t$id\t$lane\tPR\t$ext" >> "$output_path"
        # Increment the task ID counter
        taskid=$((taskid + 1))
    fi
done
