#!/bin/bash
#this script will take vcf files and filter out invariant sites and and sites marked as repeates, and generate depth and missingness statistics for variant sites only
source activate vcf

# Get the chromosome number from the command line argument
bcf=${1}
chr=${2}
input=${3}

repeat_masked_bcf="variant_filtering/${input}"
variant_bcf="variant_filtering/${bcf}_${chr}_variant.bcf"
invariant_bcf="variant_filtering/${bcf}_${chr}_invariant.bcf"

#filter out the repeat masked sites
echo "Filter variants sites for chrom  in ${repeat_masked_bcf}"
bcftools view -f PASS,. --min-ac 1 --threads 4 -Ob -o $variant_bcf $repeat_masked_bcf
echo "done generating filtered bcf ${variant_bcf}"

echo "number of sites in ${variant_bcf}. SNPS after filtering repeats"
bcftools +counts ${variant_bcf}

#site depth total
vcftools --bcf ${variant_bcf} --site-depth --out variant_filtering/${bcf}_${chr}_variant
#mean site depth
vcftools --bcf ${variant_bcf} --site-mean-depth --out variant_filtering/${bcf}_${chr}_variant
#mean depth individual
vcftools --bcf ${variant_bcf} --depth --out variant_filtering/${bcf}_${chr}_variant
#missing sites per individuals
vcftools --bcf ${variant_bcf} --missing-indv --out variant_filtering/${bcf}_${chr}_variant

echo "done with depth and missing checks for variant sites"
