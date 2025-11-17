#!/bin/bash

#this script will take bcf and filter by missingness and MAF. 
###############################################################################
# USAGE: filter_variants_MISyMAF3.sh  <input.vcf.gz>  <basename_out>  | & tee job.log
###############################################################################

input=${1} 
output=${2} 

#activate environmnet before
source activate vcf
echo "Starting"

# Put bash in strict mode
set -eo pipefail

echo "Filtering sites with Missing genotypes>=20%"
bcftools filter --soft-filter MISS --mode +  -e 'F_MISSING>=0.2' --threads 8  "${input}" | bcftools view -f .,PASS -Ob -o "${output}_BDQM.bcf"
bcftools +counts "${output}_BDQM.bcf"
################################################################
echo "Filtering sites with MAF < 3%"
bcftools filter -e 'MAF < 0.03' --threads 8 -Ob -o "${output}_BDQMF.bcf" "${output}_BDQM.bcf"
vcftools --bcf  ${output}_BDQMF.bcf --out ${output}_BDQMF --missing-indv
echo " Remaining sites."
bcftools +counts "${output}_BDQMF.bcf"
echo "done with MAF filtering"
echo "done with filtering"