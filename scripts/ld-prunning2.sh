#!/bin/bash

#script to perform LD pruning using plink, recode vcf and perform pca
#load modules
module load miniconda3
conda activate culex
module load plink

#set variables
input=${1} # file.bcf
prefix=${2} #pipiens123_ch1

# perform linkage pruning - i.e. identify prune sites
echo "performing prunning on ${input}"

plink --vcf  $input --double-id \
--set-missing-var-ids @:# \
--indep-pairwise 200 10 0.2 --out $prefix

echo "prunning complete"
echo "recoding vcf"

plink --vcf $input --set-missing-var-ids @:# \
  --extract ${prefix}.prune.in \
  --recode vcf-iid \
  --out ${prefix}_pruned

echo "performing pca from recoded vcf"

plink --vcf ${prefix}_pruned.vcf \
      --double-id \
      --set-missing-var-ids @:# \
      --make-bed \
      --pca \
      --out ${prefix}_PCA