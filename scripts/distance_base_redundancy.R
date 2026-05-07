##plots fro SEG2025 about cx pipines stuctures.
#take into account that to work with this i willl have to be changing session dir
#from on plot to the nex

###INCLUDES 3 MODELS
## 1) MODEL WITH OCNM 1 AND 2
## 2) MODEL WITH PCNM 1-5
## 3) MODEL WITH PCNMs As a block to test geography in general
## 4) mantel-test separado por ecotipos

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

expl <- meta[, c("habitat", "pQ1", "Wolbachia", "location")]
rownames(expl) <- meta$ID
str(expl)

#----  1) dbRDA with PCNM 1 and 2--------------#
#before doing it, select only geographic vectors that are significant, to reduce noise
mod <- capscale(gen_dist ~ scores(pcnm_geo))
anov <- anova(mod, by = "axis")
anov
#none of them are significative, but as i have previous evidence of a weak relation
#given by mantel test, i will stil include first 2 vectors

geo_sig <- as.data.frame(scores(pcnm_geo, choices = 1:2))

#check same order in ids
geo_sig <- geo_sig[rownames(expl), ] #reorder so they match
all(rownames(expl) == rownames(geo_sig))  # should be TRUE

expl_full <- cbind(expl, geo_sig) # bind tables

model <- capscale(gen_dist ~ habitat + pQ1 + Wolbachia + PCNM1 + PCNM2, data = expl_full)

summary(model)

# Test full model
anova(model, permutations = 999)

# Test each variable (habitat, admixture, PCNM1, PCNM2)
anova(model, by = "margin", permutations = 999)

plot(model, display = "sites", type = "n")
points(model, display = "sites", col = expl_full$admixture, pch = 19)
legend("topright", legend = levels(expl_full$admixture), col = 1:3, pch = 19)


####--------------- 2) try a model with more PCNM------------------####
geo_sig2 <- as.data.frame(scores(pcnm_geo, choices = 1:5))
scores(pcnm_geo)
#check same order in ids
geo_sig2 <- geo_sig2[rownames(expl), ] #reorder so they match
all(rownames(expl) == rownames(geo_sig2))  # should be TRUE

expl_full2 <- cbind(expl, geo_sig2) # bind tables

model2 <- capscale(gen_dist ~ habitat + pQ1 + Wolbachia + PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5, data = expl_full2)

summary(model2)

# Test full model
anova(model2, permutations = 999)

# Test each variable (habitat, admixture, PCNM1, PCNM2)
anova(model2, by = "margin", permutations = 999)

#-------------plot new model--------------------
# Extract CAP axes (site scores)
cap_scores <- scores(model2, display = "sites", choices = 1:2)

# Turn into a dataframe
cap_df <- as.data.frame(cap_scores)
cap_df$ID <- rownames(cap_scores)

# Merge with metadata (expl_full2 must have rownames = IDs)
cap_df <- cbind(cap_df, expl_full2[rownames(cap_df), ])

#get % variance explained
summary(model2)$CCA
eig_vals <- summary(model2)$CCA$eig
cap1_pct <- round(100 * eig_vals[1] / sum(eig_vals), 1)
cap2_pct <- round(100 * eig_vals[2] / sum(eig_vals), 1)

#plot
p2 <- ggplot(cap_df, aes(x = CAP1, y = CAP2, color = pQ1)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "dbRDA – Genetic structure (CAP1 vs CAP2)",
    x = paste0("CAP1 (", cap1_pct, "%)"),
    y = paste0("CAP2 (", cap2_pct, "%)")
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(
    legend.title = element_blank(),
    panel.grid = element_line(color = "gray90")
  )


#----------- save model results-------------------#
ggsave("dbRDA-pipiens1234_Wolbachia.png", device = "png", units = "in", plot=p2, width=6, height=5)
saveRDS(model2, file = "dbrda_pipiens1234_Wolbachia.rds")
#anova reuslts
sink("dbRDA_pipiens1234_Wolbachia_summary.txt")
summary(model2)
sink()
sink("dbRDA_pipiens1234_Wolbachia_summary.tx")
print(summary(model2))
print(anova(model2, permutations = 999))
print(anova(model2, by = "margin", permutations = 999))
sink()

#key data for plotting
site_scores <- scores(model2, display = "sites")
write.csv(site_scores, "dbRDA_pipiens1234_Wolbachia_site_scores.csv")
#eigenvals and explained variance
eig_vals <- summary(model2)$CCA$eig
var_expl <- round(100 * eig_vals / sum(eig_vals), 2)
write.csv(data.frame(CAP = names(eig_vals), Eigenvalue = eig_vals, VariancePct = var_expl),
          "dbRDA_pipiens1234_Wolbachia_variance.csv", row.names = FALSE)

####--------------- 3) model with geo block------------------####
geo_sig2 <- as.data.frame(scores(pcnm_geo, choices = 1:5))
scores(pcnm_geo)
#check same order in ids
geo_sig2 <- geo_sig2[rownames(expl), ] #reorder so they match
all(rownames(expl) == rownames(geo_sig2))  # should be TRUE

expl_full2 <- cbind(expl, geo_sig2) # bind tables


geo_block <- scores(pcnm_geo, choices = 1:5)
str(geo_block)
model_geo <- capscale(gen_dist ~ geo_block, data = expl_full2)
anova(model_geo, permutations = 9999)
model_geo

model_geo2 <- capscale(gen_dist ~ habitat + pQ1 + Wolbachia + geo_block, data = expl_full2)
anova_geo2 <- anova(model_geo2, permutations = 9999)

summary(model_geo2)
anova_geo2
anova(model_geo2, by = "margin", permutations = 9999)


#----------- save model results-------------------#
saveRDS(model_geo2, file = "dbrda_pipiens1234_Wolbachia_geoBlock_10kperm.rds")
sink("dbRDA_pipiens1234_Wolbachia_geoBlock_10kperm_summary.txt")
print(summary(model_geo2))
print(anova(model_geo2, permutations = 9999))
print(anova(model_geo2, by = "margin", permutations = 9999))
sink()
site_scores <- scores(model_geo2, display = "sites")
write.csv(site_scores, "dbRDA_pipiens1234_Wolbachia_geoBlock_10kperm_site_scores.csv")
eig_vals <- summary(model_geo2)$CCA$eig
var_expl <- round(100 * eig_vals / sum(eig_vals), 2)
write.csv(
  data.frame(CAP = names(eig_vals), Eigenvalue = eig_vals, VariancePct = var_expl),
  "dbRDA_pipiens1234_Wolbachia_geoBlock_10kperm_variance.csv",
  row.names = FALSE
)

#### 2) Partial dbRDA: habitat effect controlling for pQ1 + space ####
model_partial <- capscale(
  gen_dist ~ habitat + Condition(pQ1 + Wolbachia + geo_block),
  data = expl_full2
)

summary(model_partial)

anova(model_partial, permutations = 9999)
anova(model_partial, by = "margin", permutations = 9999)

#----------- save model results-------------------#
saveRDS(model_partial, file = "partialdbrda_pipiens1234_habitat_10kperm.rds")
sink("partialdbRDA_pipiens1234_habitat_10kperm_summary.txt")
print(summary(model_partial))
print(anova(model_partial, permutations = 9999))
print(anova(model_partial, by = "margin", permutations = 9999))
sink()
site_scores <- scores(model_partial, display = "sites")
write.csv(site_scores, "partialdbRDA_pipiens1234_habitar_10kperm_site_scores.csv")
eig_vals <- summary(model_partial)$CCA$eig
var_expl <- round(100 * eig_vals / sum(eig_vals), 2)
write.csv(
  data.frame(CAP = names(eig_vals), Eigenvalue = eig_vals, VariancePct = var_expl),
  "partialdbRDA_pipiens1234_habitat_10kperm_variance.csv",
  row.names = FALSE
)


#####------------------check colinearity................##
# -------------------------------
# 1. Correlation among geo_block PCNMs
# -------------------------------
cor(geo_block)  # Look for values >0.7 or <−0.7 → potential collinearity

# -------------------------------
# 2. Check association between geo_block (PCNMs) and admixture/habitat
# -------------------------------

# Boxplots of PCNM values by ecotype or habitat
par(mfrow = c(2, 3))
for (i in 1:5) {
  boxplot(geo_block[, i] ~ expl_full2$habitat,
          main = paste("PCNM", i, "~explanatory"))
}

par(mfrow = c(1, 1))  # reset layout

# ANOVA per PCNM ~ admixture
apply(geo_block, 2, function(pcnm) summary(aov(pcnm ~ expl_full2$pQ1)))

# ANOVA per PCNM ~ habitat
apply(geo_block, 2, function(pcnm) summary(aov(pcnm ~ expl_full2$habitat)))

apply(geo_block, 2, function(pcnm) summary(aov(pcnm ~ expl_full2$Wolbachia)))

?aov

# -------------------------------
# 3. Variance Inflation Factor (VIF)
# -------------------------------
# Run this after fitting your capscale model (e.g. model_geo2)
library(vegan)
vif.results <- vif.cca(model_geo2)  # VIF >10 suggests problematic collinearity
vif.results
?vif.cca
# -------------------------------
# Interpretation notes:
# -------------------------------
# - Strong correlation between predictors (numeric or factor) can inflate variance.
# - High VIF → consider removing or conditioning on that variable.
# - PCNMs highly correlated with ecotype/habitat → consider using `Condition()` in capscale.


##===============mante test por ecotipos =================#
metadata <- read.table("metadata-pipiens-217samples.txt", sep="\t", header=T)
rownames(metadata) <- metadata$ID

# Seleccionar muestras pipiens
pipiens_ids <- rownames(metadata)[metadata$ecot60 == "pipiens"]

geo_pipiens <- geo_matrix[pipiens_ids, pipiens_ids]
gen_pipiens <- gen_dist[pipiens_ids, pipiens_ids]

# Seleccionar muestras molestus
molestus_ids <- rownames(metadata)[metadata$ecot60 == "molestus"]

geo_molestus <- geo_matrix[molestus_ids, molestus_ids]
gen_molestus <- gen_dist[molestus_ids, molestus_ids]

# Mantel test pipiens
mantel_pipiens <- mantel(geo_pipiens, gen_pipiens, method="pearson", permutations=9999)
mantel_pipiens
# Mantel test molestus
mantel_molestus <- mantel(geo_molestus, gen_molestus, method="pearson", permutations=9999)
mantel_molestus


### -------------------------------------------###
##    modelo para ver si urbano tiene + ancestria molestus ##
### -------------------------------------------##
# Histogram with kernel density

ggplot(meta, aes(x = pQ1)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density()

qqnorm(meta$pQ1, main = "pQ1")
qqline(meta$pQ1)

shapiro.test(meta_mol$pQ2)
# W = 0.76616, p-value < 2.2e-16
#data is normal when p>0.05

#check same order in ids
#### MODEL #######
#expl_full2 is defined above in this same script

#install.packages("glmmTMB") 
library(glmmTMB)

mod_Q <- glmmTMB(
  pQ1 ~ habitat + (1 | location),
  data = expl_full2,
  family = beta_family(link = "logit")
)

summary(mod_Q)

#diagnostics
#install.packages("DHARMa") 
library(DHARMa)

sim <- simulateResiduals(mod_Q)
plot(sim)

#testDispersion(sim) dont save
#testUniformity(sim) dont save

#get coeeficients in orginal scale
#now are in logit scales, if i want proportions:
#install.packages("emmeans")
library(emmeans)

emmeans(mod_Q, ~ habitat, type = "response")


#### ================= SAVE results ================= ####

#### 1) SAVE MODEL OBJECT ####
saveRDS(mod_Q, file = "mod_pQ1_beta_glmm.rds")

#### 2) SAVE MODEL SUMMARY ####
sink("summary_mod_pQ1.txt")
summary(mod_Q)
sink()

#### 3) SAVE CLEAN COEFFICIENT TABLE ####
#install.packages("broom.mixed") # run once if needed
library(broom.mixed)

tidy_mod <- tidy(mod_Q)
write.csv(tidy_mod,
          file = "model_coefficients_mod_pQ1.csv",
          row.names = FALSE)

#### 4) SAVE EMMEANS (ON RESPONSE SCALE) ####
emm_Q <- emmeans(mod_Q, ~ habitat, type = "response")
emm_df <- as.data.frame(emm_Q)

write.csv(emm_df,
          file = "emmeans_mod_pQ1.csv",
          row.names = FALSE)

saveRDS(emm_Q, file = "emmeans_mod_pQ1.rds")

#### 5) SAVE DHARMa OBJECT ####
saveRDS(sim, file = "DHARMa_sim_mod_pQ1.rds")

#### 6) SAVE DHARMa RESIDUAL PLOT #### done by habnd from the display
