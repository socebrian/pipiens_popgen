#!/bin/bash
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL


module load miniconda3
conda activate culex

# Exit immediately on error
set -e
set -o pipefail

# Create header for the output file
echo -e "ID\tTotalSeq\tPairedSeq\tPctMapped\tProperlyPaired\tPctProperlyPaired" > summary.txt

# Loop through all .stats files
for file in *.stats
do
  # Extract sample ID from filename (strip off .stats)
  ID=$(basename "$file" .stats)

  # 1) Total reads: first line -> "NNNNNNNN + 0 in total ..."
  totalSeq=$(awk 'NR==1 {print $1}' "$file")

  # 2) Paired in sequencing: line containing "... paired in sequencing"
  pairedSeq=$(awk '/paired in sequencing/ {print $1}' "$file")

  # 3) Mapped %: line containing "... mapped (XX.XX% ...)"
  #    We capture the numeric percentage inside parentheses
  pctMapped=$(awk '/ mapped \(/ {
      match($0,/\(([0-9]+\.[0-9]+)%/,m);
      if (m[1] != "") {
          print m[1];
          exit;
      }
    }' "$file")

  # 4) Properly paired reads: line containing "... properly paired (XX.XX% ...)"
  properlyPaired=$(awk '/properly paired \(/ {print $1}' "$file")

  # 5) Properly paired %: same line as above, capture the percentage
  pctProperlyPaired=$(awk '/properly paired \(/ {
      match($0,/\(([0-9]+\.[0-9]+)%/,m);
      if (m[1] != "") {
          print m[1];
          exit;
      }
    }' "$file")

  # Append one row per file
  echo -e "${ID}\t${totalSeq}\t${pairedSeq}\t${pctMapped}\t${properlyPaired}\t${pctProperlyPaired}" >> summary.txt
done

echo "Summary written to summary.txt"
