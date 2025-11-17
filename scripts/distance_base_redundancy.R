## Sonia Cebrián Camisón
## Proyecto Culex SNPs
# Script to run PCA distance based redundancy analyses (dbRDA)

# Load packages
library(vegan)
library(ape)
library(dplyr)
library(tidyr)
library(reshape2) 
library(ggplot2)

#----Load data and transform as needed----
geo_dist <- read.csv("geo-distances-pipiens1234.csv" , header = TRUE) 
gen_dist <- read.table("pipiens1234_1-IBS.mdist")
gen_dist_ids <- read.table("pipiens1234_1-IBS.mdist.id")
metadata <- read.table("metadata-pipiens-217samples.txt", stringsAsFactors = TRUE, header=T)


#--- genetic distance ---# transform to matrix
gen_mat <- as.matrix(gen_dist)
gen_dist <- as.dist(gen_mat) 
str(gen_mat)

#---geographic data-----# transform to vectors
#transform in matrix
geo_matrix <- dcast(geo_dist, InputID ~ TargetID, value.var="Distance")

#fix format and colnames
rownames(geo_matrix) <- geo_matrix$InputID
geo_matrix$InputID <- NULL 
geo_matrix <- as.matrix(geo_matrix)
str(geo_matrix)

# Copy upper triangle to lower triangle to make symmetric
geo_matrix[lower.tri(geo_matrix)] <- t(geo_matrix)[lower.tri(geo_matrix)]
all.equal(geo_matrix, t(geo_matrix))  # should return TRUE
#convert to proper dist object
geo_dist_matrix <- as.dist(geo_matrix)
str(geo_dist_matrix)

#Convert geographic distance matrix to PCNM spatial variables
pcnm_geo <- pcnm(geo_dist_matrix)
scores(pcnm_geo)

#----metadata ----# create table with 2 explanatory variables same order as matrix
#  extract sample IDs from rownames
sample_ids <- rownames(gen_mat)

# Subset and reorder metadata accordingly
meta <- metadata[metadata$ID %in% sample_ids, ]
meta <- meta[match(sample_ids, meta$ID), ]

# Quick check
if (!all(meta$ID == sample_ids)) {
  warning("Mismatch in sample ID order between metadata and gen_mat!")
} else {
  message("✔ Sample order matches.")
}

expl <- meta[, c("habitat", "pQ1", "Wolbachia")]
rownames(expl) <- meta$ID
str(expl)

####--------------- model with geo block ------------------####
geo_block <- scores(pcnm_geo, choices = 1:5)
str(geo_block)
model_geo <- capscale(gen_dist ~ geo_block, data = expl_full2)
anova(model_geo, permutations = 999)
model_geo

model_geo2 <- capscale(gen_dist ~ habitat + pQ1 + Wolbachia + geo_block, data = expl_full2)
anova_geo2 <- anova(model_geo2, permutations = 999)

summary(model_geo2)
anova_geo2
anova(model_geo2, by = "margin", permutations = 999)

#----------- save model results-------------------#
saveRDS(model_geo2, file = "dbrda_pipiens1234_Wolbachia_geoBlock.rds")
sink("dbRDA_pipiens1234_Wolbachia_geoBlock_summary.txt")
print(summary(model_geo2))
print(anova(model_geo2, permutations = 999))
print(anova(model_geo2, by = "margin", permutations = 999))
sink()
site_scores <- scores(model_geo2, display = "sites")
write.csv(site_scores, "dbRDA_pipiens1234_Wolbachia_geoBlock_site_scores.csv")
eig_vals <- summary(model_geo2)$CCA$eig
var_expl <- round(100 * eig_vals / sum(eig_vals), 2)
write.csv(
  data.frame(CAP = names(eig_vals), Eigenvalue = eig_vals, VariancePct = var_expl),
  "dbRDA_pipiens1234_Wolbachia_geoBlock_variance.csv",
  row.names = FALSE
)

#------------------check colinearity................##
# 1. Correlation among geo_block PCNMs
cor(geo_block)  # Look for values >0.7 or <−0.7 → potential collinearity

# 2. Check association between geo_block (PCNMs) and Wolbachia/habitat
# Boxplots of PCNM values by habitat or Wolbachia infection status 
par(mfrow = c(2, 3))
for (i in 1:5) {
  boxplot(geo_block[, i] ~ expl_full2$pQ1,
          main = paste("PCNM", i, "~ pQ1"))
}
# ANOVA per PCNM ~ admixture
apply(geo_block, 2, function(pcnm) summary(aov(pcnm ~ expl_full2$pQ1)))

# ANOVA per PCNM ~ habitat
apply(geo_block, 2, function(pcnm) summary(aov(pcnm ~ expl_full2$habitat)))

# ANOVA per PCNM ~ Wolbachia
apply(geo_block, 2, function(pcnm) summary(aov(pcnm ~ expl_full2$Wolbachia)))


# -------------------------------
# 3. Variance Inflation Factor (VIF)
# -------------------------------
library(vegan)
vif.results <- vif.cca(model_geo2)  # VIF >10 suggests problematic collinearity
vif.results

# -------------------------------
# Interpretation notes:
# -------------------------------
# - Strong correlation between predictors (numeric or factor) can inflate variance.
# - High VIF → consider removing or conditioning on that variable.
# - PCNMs highly correlated with ecotype/habitat → consider using `Condition()` in capscale.
