
#### Purpose: evalute the SE model predictions

library(dplyr)
library(data.table)
library(ggplot2)
library(pROC)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(grid)
library(gridExtra)
library(binom)
library(glue)
source("src/00_functions.R")

# ====================== #
# Load and process data
# ====================== #

data_SE = fread("results/SE_STAN_AllPairs_SD=3.75_OMULT=0.6_TestPred.csv")
data_SE_similar = fread("results/SE_STAN_PhenoSimilar_SD=3.75_OMULT=0.6_TestPred.csv")
data_SE_dissimilar = fread("results/SE_STAN_PhenoDissimilar_SD=3.75_OMULT=0.6_TestPred.csv")

diseases_names = fread("data/gwas_matching_between_sources_manual_final_2.txt")
diseases_names = diseases_names[,c("mesh_id","mesh_name", "trait_category")]
data_SE = cbind(
  stringr::str_split_fixed(data_SE$V1, pattern = ":", n = 2),
  data_SE[,c("p","label")]
) %>%
  dplyr::select(drug = V1, phenotype = V2, PredProb = p, is_side_effect = label) %>%
  left_join(diseases_names, by = c("phenotype" = "mesh_id")) %>%
  dplyr::select(drug, phenotype, phenotype_name = mesh_name, phenotype_category = trait_category, PredProb, is_side_effect)
data_SE_similar = cbind(
  stringr::str_split_fixed(data_SE_similar$V1, pattern = ":", n = 2),
  data_SE_similar[,c("p","label")]
) %>%
  dplyr::select(drug = V1, phenotype = V2, PredProb = p, is_side_effect = label) %>%
  left_join(diseases_names, by = c("phenotype" = "mesh_id")) %>%
  dplyr::select(drug, phenotype, phenotype_name = mesh_name, phenotype_category = trait_category, PredProb, is_side_effect)
data_SE_dissimilar = cbind(
  stringr::str_split_fixed(data_SE_dissimilar$V1, pattern = ":", n = 2),
  data_SE_dissimilar[,c("p","label")]
) %>%
  dplyr::select(drug = V1, phenotype = V2, PredProb = p, is_side_effect = label) %>%
  left_join(diseases_names, by = c("phenotype" = "mesh_id")) %>%
  dplyr::select(drug, phenotype, phenotype_name = mesh_name, phenotype_category = trait_category, PredProb, is_side_effect)

#### evaluation ####

# ====================================== #
# AUROC on test set (held-out drugs)
# ====================================== #

# both
roc(is_side_effect ~ PredProb, data_SE) # 0.5818
# similar
roc(is_side_effect ~ PredProb, data_SE_similar) # 0.6023
# dissimilar
roc(is_side_effect ~ PredProb, data_SE_dissimilar) # 0.5483

# ======================================= #
# Per phenotype AUROC on test set
# ======================================= #

per_pheno_SE_AUROC = data.frame(
  phenotype = unique(data_SE_dissimilar$phenotype_name),
  phenotype_category = NA,
  AUROC_test = NA,
  Nr_ind_drugs = NA
)
for (i in 1:nrow(per_pheno_SE_AUROC)) {
  pheno = per_pheno_SE_AUROC[i, "phenotype"]
  pheno_data = data_SE_dissimilar %>% 
    filter(phenotype_name == pheno) %>%
    dplyr::select(phenotype_name, phenotype_category, drug, is_side_effect, PredProb) %>%
    distinct()
  pheno_data$is_side_effect = as.numeric(as.character(pheno_data$is_side_effect))
  if (sum(pheno_data$is_side_effect) == 0) {
    per_pheno_SE_AUROC[i, "AUROC_test"] = NA
    per_pheno_SE_AUROC[i, "Nr_ind_drugs"] = 0
    next
  }
  per_pheno_SE_AUROC[i, "AUROC_test"] = roc(pheno_data$is_side_effect, pheno_data$PredProb, levels = c("0", "1"), direction = "<")$auc
  per_pheno_SE_AUROC[i, "phenotype_category"] = unique(pheno_data$phenotype_category)
  per_pheno_SE_AUROC[i, "Nr_ind_drugs"] = sum(pheno_data$is_side_effect)
} ; rm(i, pheno, pheno_data)

per_pheno_SE_AUROC = per_pheno_SE_AUROC %>% arrange(phenotype_category, AUROC_test)
per_pheno_SE_AUROC = na.omit(per_pheno_SE_AUROC) ; rownames(per_pheno_SE_AUROC) = NULL
per_pheno_SE_AUROC$phenotype = factor(per_pheno_SE_AUROC$phenotype, levels = per_pheno_SE_AUROC$phenotype, labels = per_pheno_SE_AUROC$phenotype)
per_pheno_SE_AUROC$phenotype_category = factor(per_pheno_SE_AUROC$phenotype_category, levels = per_pheno_SE_AUROC$phenotype_category, labels = per_pheno_SE_AUROC$phenotype_category)

my_colors = c(
  rep("cornflowerblue", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[1])),
  rep("#FB9A99", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[2])),
  rep("orange", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[3])),
  rep("red2", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[4])),
  rep("darkturquoise", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[5])),
  rep("tan4", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[6])),
  rep("coral1", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[7])),
  rep("green4", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[8])),
  rep("magenta1", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[9])),
  rep("tan2", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[10])),
  rep("dodgerblue2", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[11])),
  rep("brown", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[12])),
  rep("cyan3", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[13])),
  rep("black", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[14])),
  rep("orchid", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[15])),
  rep("darkgoldenrod2", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[16])),
  rep("wheat4", as.numeric(table(per_pheno_SE_AUROC$phenotype_category)[17]))
)

ggplot(per_pheno_SE_AUROC, aes(x = AUROC_test, y = phenotype, size = Nr_ind_drugs)) +
  geom_point(color = my_colors) +
  geom_vline(xintercept = 0.5, color = "black", linetype = "solid", linewidth = 0.7) +
  labs(
    title = "SIDE EFFECTS model",
    subtitle = "AUROC per phenotype",
    x = "AUROC\n(held-out drugs)",
    y = "",
    size = "# side effect drugs"
  ) +
  scale_x_continuous(
    breaks = seq(0, 2, 0.1)
  ) +
  coord_cartesian(xlim = c(0.05, 1), clip = "off") +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 24, color = "black", face = "bold"),
    plot.subtitle = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 9, color = my_colors),
    legend.title = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 18, color = "black"),
    plot.margin = margin(5.5, 5.5, 5.5, 200)   # extra left margin so labels show
  ) +
  annotate("text", x = -0.46, y = 1.5, label = "Blood disorders", size = 4.5, color = "cornflowerblue", hjust = "right") +
  annotate("text", x = -0.46, y = 14.5, label = "Cardiovasular diseases", size = 4.5, color = "#FB9A99", hjust = "right") +
  annotate("text", x = -0.46, y = 34.5, label = "Cell proliferation diseases", size = 4.5, color = "orange", hjust = "right") +
  annotate("text", x = -0.46, y = 44.5, label = "Ear diseases", size = 4.5, color = "red2", hjust = "right") +
  annotate("text", x = -0.46, y = 49, label = "Visual system diseases", size = 4.5, color = "darkturquoise", hjust = "right") +
  annotate("text", x = -0.46, y = 53, label = "Endocrine system diseases", size = 4.5, color = "tan4", hjust = "right") +
  annotate("text", x = -0.46, y = 67.5, label = "Gastrointestinal diseases", size = 4.5, color = "coral1", hjust = "right") +
  annotate("text", x = -0.46, y = 81, label = "Immune system disease", size = 4.5, color = "green4", hjust = "right") +
  annotate("text", x = -0.46, y = 87.5, label = "Integumentary system diseases", size = 4.5, color = "magenta1", hjust = "right") +
  annotate("text", x = -0.46, y = 94.5, label = "Metabolic diseases", size = 4.5, color = "tan2", hjust = "right") +
  annotate("text", x = -0.46, y = 100.5, label = "Musculoskeletal diseases", size = 4.5, color = "dodgerblue2", hjust = "right") +
  annotate("text", x = -0.46, y = 112.5, label = "Nervous system diseases", size = 4.5, color = "brown", hjust = "right") +
  annotate("text", x = -0.46, y = 123, label = "Pancreas diseases", size = 4.5, color = "cyan3", hjust = "right") +
  annotate("text", x = -0.46, y = 125, label = "Psychiatric disease", size = 4.5, color = "black", hjust = "right") +
  annotate("text", x = -0.46, y = 126.5, label = "Reproductive/Breast diseases", size = 4.5, color = "orchid", hjust = "right") +
  annotate("text", x = -0.46, y = 132.5, label = "Respiratory diseases", size = 5, color = "darkgoldenrod2", hjust = "right") +
  annotate("text", x = -0.46, y = 141, label = "Urinary system diseases", size = 5, color = "wheat4", hjust = "right")
ggsave("figures/Figure_5a.png",
       device = "png", dpi = 600, width = 14, height = 18)
sum(per_pheno_SE_AUROC$AUROC_test > 0.5) / nrow(per_pheno_SE_AUROC)

# =============================================== #
# SIDER data enrichment for our predictions
# =============================================== #

## SIDER data
sider = fread("data/SIDER_processed.txt", sep="\t")
# filter for drugs and phenotypes in our data
sider = sider %>% 
  filter(
    drug %in% data_SE_dissimilar$drug,
    mesh_id %in% data_SE_dissimilar$phenotype
  ) %>%
  mutate(sider_se = 1)
# combine our data with sider data
ultimate = data_SE_dissimilar %>%
  dplyr::select(drug, phenotype, is_side_effect, PredProb) %>%
  filter(phenotype %in% sider$mesh_id, drug %in% sider$drug) %>%
  left_join(sider, by = c("drug", "phenotype" = "mesh_id"))
ultimate$sider_se = ifelse(is.na(ultimate$sider_se), 0, 1)
# check overlap between onSIDES and SIDER --> expected to be high
table(onSIDES = ultimate$is_side_effect, SIDER = ultimate$sider_se)
fisher.test(table(onSIDES = ultimate$is_side_effect, SIDER = ultimate$sider_se))
# keep only drug-phenotype pairs that are reported as side effects in SIDER and not in onSIDES (for more fair evaluation of our predictions)
ultimate = ultimate %>% filter(is_side_effect == 0 )

## calculate enrichment for likely true side effects (odds ratio) across different probability thresholds (compare to pairs with probability lower than the baseline frequency)
baseline_prob = prop.table(table(data_SE_dissimilar$is_side_effect))["1"]
baseline_prob

SE_ORs = data.frame(
  threshold = c(baseline_prob, 0.1, 0.2, 0.3, 0.4, 0.5),
  OR = NA,
  OR_ci95_lower = NA,
  OR_ci95_upper = NA,
  nr_predicted_pairs = NA
)
for (i in 1:nrow(SE_ORs)) {
  threshold = SE_ORs[i,"threshold"]
  temp = ultimate %>%
    mutate(predicted = if_else(PredProb > threshold, 1, 2)) %>% # 0: not predicted | 1: predicted
    mutate(predicted = if_else(PredProb < baseline_prob, 0, predicted)) %>%
    filter(predicted %in% c(0,1))
  fisher_or = fisher.test(table(temp$predicted, temp$sider_se))
  SE_ORs[i, "OR"] = fisher_or$estimate
  SE_ORs[i, "OR_ci95_lower"] = fisher_or$conf.int[1]
  SE_ORs[i, "OR_ci95_upper"] = fisher_or$conf.int[2]
  SE_ORs[i, "nr_predicted_pairs"] = sum(temp$predicted == 1)
}
# plot
SE_ORs$threshold = round(SE_ORs$threshold, 3)
SE_ORs$threshold = gsub(round(baseline_prob, 3), paste(round(baseline_prob,3), "\n(baseline probability)", sep = ""), SE_ORs$threshold)
SE_ORs$threshold = factor(SE_ORs$threshold, levels = unique(SE_ORs$threshold), labels = paste(">", unique(SE_ORs$threshold), sep = ""))
ggplot(SE_ORs, aes(x = threshold, y = OR, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = nr_predicted_pairs), show.legend = FALSE) +
  geom_ribbon(aes(ymin = OR_ci95_lower, ymax = OR_ci95_upper), alpha = 0.2, color = NA) +
  # geom_text(aes(label = round(OR,2)), vjust = -1, size = 10, data = SE_ORs %>% filter(OR_ci95_lower > 1)) +
  geom_text(aes(label = scales::comma(nr_predicted_pairs)), nudge_y = 0.1, size = 6) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 1) +
  labs(
    title = "", 
    x = "Predicted probability",
    y = "Odds ratio"
  ) +
  theme_bw(base_size = 20, base_family = "Arial") +
  theme(
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(size = 35, family = "Arial", color = "black", margin = margin(r = 20)),
    axis.title.x = element_text(size = 35, family = "Arial", color = "black", margin = margin(t = 20)),
    axis.text = element_text(size = 20, family = "Arial", color = "black")
  )
ggsave("figures/Figure_5b.png",
       device = "png", dpi = 600, width = 12, height = 6)


# =============================================== #
# offSIDES data enrichment for our predictions
# =============================================== #

## offSIDES data
# this data is provided to us only for the purpose of our analysis and we cannot upload them to our GitHub repository
# for access to this data, please contact Nicholas Tatonetti (Nicholas.Tatonetti@cshs.org)
data_offsides = fread("data/offSIDES_ALL_processed.txt") # won't run
data_offsides = data_offsides %>%
  # filter for drugs and phenotypes in our sample
  filter(chembl_id %in% data_SE_dissimilar$drug) %>%
  filter(mesh_id %in% data_SE_dissimilar$phenotype) %>%
  distinct() %>%
  # filter for drug-side effect pairs with PRR > 1
  filter(PRR > 1) %>%
  # remove drug-phenotype pairs with duplicates (duplicates occur when I match original offSIDES identifiers to ChEMBL for drugs and MeSH for side effects)
  group_by(mesh_id, chembl_id) %>%
  filter(n() == 1) %>%
  ungroup()

## Add predicted probabilities from the the two side effect models
ultimate = data_SE_dissimilar %>%
  dplyr::select(drug, phenotype, is_side_effect, PredProb) %>%
  left_join(data_offsides, by = c("drug" = "chembl_id", "phenotype" = "mesh_id")) %>%
  na.omit() %>% 
  mutate(sig = if_else(p.adjust(pvalue, method = "bonferroni") < 0.05, 1, 0))
# keep only pairs that are not known side effects
table(is_SE = ultimate$is_side_effect, offSIDES_sig = ultimate$sig)
fisher.test(table(is_SE = ultimate$is_side_effect, offSIDES_sig = ultimate$sig))
ultimate = ultimate %>% filter(is_side_effect == 0)

## Calculating enrichment for likely true side effects (odds ratio) across different probability thresholds (compare to pairs with probability lower than the baseline frequency of being side effect)
SE_ORs = data.frame(
  threshold = c(baseline_prob, 0.1, 0.2,  0.3, 0.4, 0.5),
  OR = NA,
  OR_ci95_lower = NA,
  OR_ci95_upper = NA,
  nr_predicted_pairs = NA
)
for (i in 1:nrow(SE_ORs)) {
  threshold = SE_ORs[i,"threshold"]
  temp = ultimate %>%
    mutate(predicted = if_else(PredProb > threshold, 1, 2)) %>% # 0: not predicted | 1: predicted
    mutate(predicted = if_else(PredProb < baseline_prob, 0, predicted)) %>%
    filter(predicted %in% c(0,1))
  fisher_or = fisher.test(table(temp$predicted, temp$sig))
  SE_ORs[i, "OR"] = fisher_or$estimate
  SE_ORs[i, "OR_ci95_lower"] = fisher_or$conf.int[1]
  SE_ORs[i, "OR_ci95_upper"] = fisher_or$conf.int[2]
  SE_ORs[i, "nr_predicted_pairs"] = sum(temp$predicted == 1)
}
# plot
SE_ORs$threshold = round(SE_ORs$threshold, 3)
SE_ORs$threshold = gsub(round(baseline_prob, 3), paste(round(baseline_prob,3), "\n(baseline probability)", sep = ""), SE_ORs$threshold)
SE_ORs$threshold = factor(SE_ORs$threshold, levels = unique(SE_ORs$threshold), labels = paste(">", unique(SE_ORs$threshold), sep = ""))
ggplot(SE_ORs, aes(x = threshold, y = OR, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = nr_predicted_pairs), show.legend = FALSE) +
  geom_ribbon(aes(ymin = OR_ci95_lower, ymax = OR_ci95_upper), alpha = 0.2, color = NA) +
  geom_text(aes(label = scales::comma(nr_predicted_pairs)), nudge_y = 0.1, size = 6) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 1) +
  labs(
    title = "", 
    x = "Predicted probability",
    y = "Odds ratio"
  ) +
  theme_bw(base_size = 20, base_family = "Arial") +
  theme(
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(size = 35, family = "Arial", color = "black", margin = margin(r = 20)),
    axis.title.x = element_text(size = 35, family = "Arial", color = "black", margin = margin(t = 20)),
    axis.text = element_text(size = 20, family = "Arial", color = "black"),
    legend.position = "none"
  )
ggsave("figures/Supplementary_figure_9.png",
       device = "png", dpi = 600, width = 12, height = 6)

##### compare predicted probabilities between "same" and "different" models

ultimate = full_join(
  data_SE_similar,
  data_SE_dissimilar,
  by = c("drug", "phenotype", "phenotype_name", "phenotype_category", "is_side_effect")
)
ultimate = na.omit(ultimate)
rownames(ultimate) = NULL
colnames(ultimate)[c(5,7)] = c("PredProb_similar", "PredProb_dissimilar")

cor.test(ultimate$PredProb_similar, ultimate$PredProb_dissimilar, method = "pearson") # 0.05
cor.test(ultimate$PredProb_similar, ultimate$PredProb_dissimilar, method = "pearson")$p.value # 1.71e-53

cor.test(ultimate$PredProb_similar, ultimate$PredProb_dissimilar, method = "spearman") # 0.06
cor.test(ultimate$PredProb_similar, ultimate$PredProb_dissimilar, method = "spearman")$p.value # 1.62e-77

roc(ultimate$is_side_effect ~ ultimate$PredProb_similar) # 0.5867 - similar
roc(ultimate$is_side_effect ~ ultimate$PredProb_dissimilar) # 0.5489 - dissimilar
