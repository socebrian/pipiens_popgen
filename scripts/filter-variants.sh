
#!/bin/bash

#this script will take vcf files, and apply some basic filters to the putput bcf files. vcf files must be indexed.
#the goal of the filters is to clean data to explore missing data and distributions
###############################################################################
# USAGE:  bash filter-variants.sh  <input.vcf.gz>  <basename_out>  <fdp>  <fgq> | & tee job.log
###############################################################################

input=${1} #pipiens123_1_filter1n_qd.vcf.gz
output=${2} #pipiens123_f1n
fdp=${3}
fgq=${4}

#activate environment before
source activate vcf
echo "Starting"

# Put bash in strict mode
set -e
set -o pipefail


###############################################################################
# helper: count sites with F_MISSING > 0.20
###############################################################################
count_fail20 () {
    local vcf=$1
    local expr='F_MISSING>=0.20'

    if bcftools +counts -h &>/dev/null ; then               # fast path
        bcftools +counts -i "${expr}" "${vcf}" |
            awk '/Number of sites/ {print $NF}'
    else                                                    # fallback
        bcftools view -i "${expr}" -H "${vcf}" | wc -l
    fi
}

###############################################################################
# helper: total number of sites in a VCF/BCF
###############################################################################
total_sites () {
    local vcf=$1
    if bcftools +counts -h &>/dev/null ; then
        bcftools +counts "${vcf}" | awk '/Number of sites/ {print $NF}'
    else
        bcftools view --threads 4 -H "${vcf}" | wc -l
    fi
}
###############################################################################

total_input=$(total_sites "$input")
echo "Total input sites: ${total_input}"

###############################################################

echo "Filtering biallelic SNPS"
bcftools view -v snps -m2 -M2 --threads 4 -Ob -o "${output}_biSNPs.bcf" "$input"
echo " Remaining $(bcftools query -f "%CHROM %POS" "${output}_biSNPs.bcf" | wc -l) sites."

total_bi=$(total_sites "${output}_biSNPs.bcf")
base_fail20=$(count_fail20 "${output}_biSNPs.bcf")
keep_base=$(( total_bi - base_fail20 ))
echo "${output}_biSNPs.bcf -- Sites remaining: ${total_bi}. Sites **passing** 20 % rule before masking: ${keep_base}"

################################################################

echo "Filtering sites with DP > 2 * mean DP"
meandp=$(bcftools query -f '%DP\n' "${output}_biSNPs.bcf" | awk '{{sum+=$1; count+=1}} END {{print sum/count}}')
echo "mean DP: $meandp"
bcftools filter --soft-filter highDP -e "INFO/DP > 2 * ${meandp}" --threads 4 "${output}_biSNPs.bcf" | bcftools view -f .,PASS -Ob -o ${output}_BinfoDP.bcf
vcftools --bcf ${output}_BinfoDP.bcf --missing-indv --out "${output}_BinfoDP"

total_highdp=$(total_sites "${output}_BinfoDP.bcf")
echo "${output}_BinfoDP.bcf -- Sites remaining: ${total_highdp}"

################################################################

echo "Filtering sites with FORMAT/DP<${fdp}"
bcftools filter -e "FORMAT/DP<${fdp}" --threads 4 ${output}_BinfoDP.bcf -S . -Ob -o ${output}_BDP.bcf
#vcftools --bcf  ${output}_BDP.bcf --out ${output}_BDP --missing-indv
echo " Remaining $(bcftools query -f "%CHROM %POS" "${output}_BDP.bcf" | wc -l) sites."

total_dp=$(total_sites "${output}_BDP.bcf")
dp_fail20=$(count_fail20 "${output}_BDP.bcf")
keep_dp=$(( total_dp - dp_fail20 ))
echo "  ${output}_BDP.bcf - Sites **passing** 20 % rule after DP mask: ${keep_dp}"

################################################################

echo "Filtering sites with FORMAT/GQ<${fgq}"
bcftools filter -e "FORMAT/GQ<${fgq}" --threads 4 ${output}_BDP.bcf -S . -Ob -o ${output}_BDQ.bcf
#vcftools --bcf  ${output}_BDQ.bcf --out ${output}_BDQ --missing-indv
echo " Remaining $(bcftools query -f "%CHROM %POS" "${output}_BDQ.bcf" | wc -l) sites."

total_gq=$(total_sites "${output}_BDQ.bcf")
gq_fail20=$(count_fail20 "${output}_BDQ.bcf")
keep_gq=$(( total_gq - gq_fail20 ))
echo "${output}_BDQ.bcf -- Sites **passing** 20 % rule after GQ mask: ${keep_gq}"
################################################################

echo "Filtering sites with Missing genotypes>=50%"
bcftools filter --soft-filter MISS --mode +  -e 'F_MISSING>=0.5' --threads 4  "${output}_BDQ.bcf" | bcftools view -f .,PASS -Ob -o "${output}_BDQM5.bcf"
vcftools --bcf  ${output}_BDQM5.bcf --out ${output}_BDQM5 --missing-indv
echo " Remaining $(bcftools query -f "%CHROM %POS" "${output}_BDQM5.bcf" | wc -l) sites."

total_miss=$(total_sites "${output}_BDQM5.bcf")
echo "${output}_BDQM.bcf -- Sites remaining after filtering out those with >50% missing genotypes: ${total_miss}"
echo "Do not need to match any of the above cause is a combination of the opriginally missing and the ones not passing DP and GT filters"

################################################################



echo "done with filtering"

#I often use intermediate files when I am unsure of the chosen filtering, you can remove them after they are use by uncommenting the following lines
#echo "removing intermediate SNPS"
#rm ${output}_BinfoDP.bcf
#rm ${output}_BDP.bcf
#rm ${output}_BDQ.bcf
#rm ${output}_BDQM.bcf

echo "done"

#Notes
# --mode + to add filter tags
# --soft-filter to add filter tags without removing the sites
# -S. OPTION PARA MARCAR Los GT AS MISSING
