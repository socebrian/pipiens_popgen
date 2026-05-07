#!/bin/bash
for bcf in *inv*.vcf.gz; do
echo "Running PCA for $bcf"
plink --vcf "$bcf" --double-id \
  --set-missing-var-ids @:# \
  --make-bed \
  --pca \
  --out "${bcf%.vcf.gz}"
done