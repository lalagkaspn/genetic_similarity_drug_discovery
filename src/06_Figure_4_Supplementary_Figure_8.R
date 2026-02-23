
#### Purpose: evaluate the IND model predictions

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(grid)
library(gridExtra)
library(pROC)
library(binom)
library(glue)
source("src/00_functions.R")

# ============ #
# Load data
# ============ #
data_I = fread("results/IND_STAN_AllPairs_SD=2.5_OMULT=0.7_TestPred.csv")
data_I_similar = fread("results/IND_STAN_PhenoSimilar_SD=2.5_OMULT=0.7_TestPred.csv")
data_I_dissimilar = fread("results/IND_STAN_PhenoDissimilar_SD=2.5_OMULT=0.7_TestPred.csv")

diseases_names = fread("data/gwas_matching_between_sources_manual_final_2.txt")
diseases_names = diseases_names[,c("mesh_id","mesh_name", "trait_category")]
data_I = cbind(
  stringr::str_split_fixed(data_I$V1, pattern = ":", n = 2),
  data_I[,c("p","label")]
) %>%
  dplyr::select(drug = V1, phenotype = V2, PredProb = p, is_indicated = label) %>%
  left_join(diseases_names, by = c("phenotype" = "mesh_id")) %>%
  dplyr::select(drug, phenotype, phenotype_name = mesh_name, phenotype_category = trait_category, PredProb, is_indicated)
data_I_similar = cbind(
  stringr::str_split_fixed(data_I_similar$V1, pattern = ":", n = 2),
  data_I_similar[,c("p","label")]
) %>%
  dplyr::select(drug = V1, phenotype = V2, PredProb = p, is_indicated = label) %>%
  left_join(diseases_names, by = c("phenotype" = "mesh_id")) %>%
  dplyr::select(drug, phenotype, phenotype_name = mesh_name, phenotype_category = trait_category, PredProb, is_indicated)
data_I_dissimilar = cbind(
  stringr::str_split_fixed(data_I_dissimilar$V1, pattern = ":", n = 2),
  data_I_dissimilar[,c("p","label")]
) %>%
  dplyr::select(drug = V1, phenotype = V2, PredProb = p, is_indicated = label) %>%
  left_join(diseases_names, by = c("phenotype" = "mesh_id")) %>%
  dplyr::select(drug, phenotype, phenotype_name = mesh_name, phenotype_category = trait_category, PredProb, is_indicated)

# ====================================== #
# AUROC on test set (held-out drugs)
# ====================================== #

# both
roc(is_indicated ~ PredProb, data_I) # 0.7468
# similar
roc(is_indicated ~ PredProb, data_I_similar) # 0.72
# dissimilar
roc(is_indicated ~ PredProb, data_I_dissimilar) # 0.5594

## For now on, I use the predicted probabilities calculated based on genetic similarity between drug indication - disease pairs that are phenotypically dissimilar

# ======================================= #
# Per phenotype AUROC on test set
# ======================================= #
per_pheno_I_AUROC = data.frame(
  phenotype = unique(data_I_dissimilar$phenotype_name),
  phenotype_category = NA,
  AUROC_test = NA,
  Nr_ind_drugs = NA
)
for (i in 1:nrow(per_pheno_I_AUROC)) {
  pheno = per_pheno_I_AUROC[i, "phenotype"]
  pheno_data = data_I_dissimilar %>% 
    filter(phenotype_name == pheno) %>%
    dplyr::select(phenotype_name, phenotype_category, drug, is_indicated, PredProb) %>%
    distinct()
  if (sum(pheno_data$is_indicated) == 0) {
    per_pheno_I_AUROC[i, "AUROC_test"] = NA
    per_pheno_I_AUROC[i, "Nr_ind_drugs"] = 0
    next
  }
  per_pheno_I_AUROC[i, "AUROC_test"] = roc(is_indicated ~ PredProb, pheno_data, levels = c("0", "1"), direction = "<")$auc
  per_pheno_I_AUROC[i, "phenotype_category"] = unique(pheno_data$phenotype_category)
  per_pheno_I_AUROC[i, "Nr_ind_drugs"] = sum(pheno_data$is_indicated)
} ; rm(i, pheno, pheno_data)

per_pheno_I_AUROC = per_pheno_I_AUROC %>% arrange(phenotype_category, AUROC_test)
per_pheno_I_AUROC = na.omit(per_pheno_I_AUROC) ; rownames(per_pheno_I_AUROC) = NULL
per_pheno_I_AUROC$phenotype = factor(per_pheno_I_AUROC$phenotype, levels = per_pheno_I_AUROC$phenotype, labels = per_pheno_I_AUROC$phenotype)
per_pheno_I_AUROC$phenotype_category = factor(per_pheno_I_AUROC$phenotype_category, levels = per_pheno_I_AUROC$phenotype_category, labels = per_pheno_I_AUROC$phenotype_category)

my_colors = c(
  rep("cornflowerblue", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[1])),
  rep("#FB9A99", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[2])),
  rep("orange", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[3])),
  rep("red2", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[4])),
  rep("darkturquoise", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[5])),
  rep("coral1", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[6])),
  rep("green4", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[7])),
  rep("magenta1", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[8])),
  rep("dodgerblue2", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[9])),
  rep("brown", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[10])),
  rep("cyan3", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[11])),
  rep("black", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[12])),
  rep("orchid", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[13])),
  rep("darkgoldenrod2", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[14])),
  rep("wheat4", as.numeric(table(per_pheno_I_AUROC$phenotype_category)[15]))
)

ggplot(per_pheno_I_AUROC, aes(x = AUROC_test, y = phenotype, size = Nr_ind_drugs)) +
  geom_point(color = my_colors) +
  geom_vline(xintercept = 0.5, color = "black", linetype = "solid", linewidth = 0.7) +
  labs(
    title = "INDICATIONS model",
    subtitle = "AUROC per phenotype",
    x = "AUROC\n(held-out drugs)",
    y = "",
    size = "# indicated drugs"
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
    axis.text.y = element_text(size = 12, color = my_colors),
    legend.title = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 18, color = "black"),
    plot.margin = margin(5.5, 5.5, 5.5, 200) # extra left margin so labels show
  ) +
  annotate("text", x = -0.46, y = 1.75, label = "Blood disorders", size = 4.5, color = "cornflowerblue", hjust = "right") +
  annotate("text", x = -0.46, y = 10.5, label = "Cardiovasular diseases", size = 4.5, color = "#FB9A99", hjust = "right") +
  annotate("text", x = -0.46, y = 25.5, label = "Cell proliferation diseases", size = 4.5, color = "orange", hjust = "right") +
  annotate("text", x = -0.46, y = 34, label = "Ear diseases", size = 4.5, color = "red2", hjust = "right") +
  annotate("text", x = -0.46, y = 35.5, label = "Visual system diseases", size = 4.5, color = "darkturquoise", hjust = "right") +
  annotate("text", x = -0.46, y = 43, label = "Gastrointestinal diseases", size = 4.5, color = "coral1", hjust = "right") +
  annotate("text", x = -0.46, y = 51, label = "Immune system diseases", size = 4.5, color = "green4", hjust = "right") +
  annotate("text", x = -0.46, y = 55.5, label = "Integumentary system disease", size = 4.5, color = "magenta1", hjust = "right") +
  annotate("text", x = -0.46, y = 63.5, label = "Musculoskeletal diseases", size = 4.5, color = "dodgerblue2", hjust = "right") +
  annotate("text", x = -0.46, y = 71.5, label = "Nervous system diseases", size = 4.5, color = "brown", hjust = "right") +
  annotate("text", x = -0.46, y = 77, label = "Pancreas disease", size = 4.5, color = "cyan3", hjust = "right") +
  annotate("text", x = -0.46, y = 78, label = "Psychiatric diseases", size = 4.5, color = "black", hjust = "right") +
  annotate("text", x = -0.46, y = 79, label = "Reproducive/Breast diseases", size = 4.5, color = "orchid", hjust = "right") +
  annotate("text", x = -0.46, y = 81.5, label = "Respiratory disease", size = 4.5, color = "darkgoldenrod2", hjust = "right") +
  annotate("text", x = -0.46, y = 84.5, label = "Urinary system diseases", size = 4.5, color = "wheat4", hjust = "right")
ggsave("figures/Figure_4A.png",
       device = "png", dpi = 600, width = 14, height = 18)
sum(per_pheno_I_AUROC$AUROC_test > 0.5) / nrow(per_pheno_I_AUROC)

# ==================================================================================================================================== #
# Clinical trial phase-specific enrichment for our predictions
# ==================================================================================================================================== #

## clinical trials data
drug_ct = fread("data/DrugDisease_PerCT.txt")
# remove pairs with "Not Applicable" phase as these are clinical trials of devices or behavioral interventions (according to ClinicalTrials.gov)
drug_ct = drug_ct %>% filter(max_phase_for_ind != "Not Applicable")

data_I_dissimilar_ct = left_join(data_I_dissimilar, drug_ct, by = c("drug", "phenotype"))
# for drugs and phenotypes that do not exist in the clinical trials file, we don't know if they have ever been in clinical trials --> set "status unknown"
shared_drugs = intersect(drug_ct$drug, data_I_dissimilar_ct$drug) # we have clinical trials data for 1,327 drugs
shared_phenos = intersect(drug_ct$phenotype, data_I_dissimilar_ct$phenotype) # and 173 phenotypes
data_I_dissimilar_ct$max_phase_for_ind = ifelse(!data_I_dissimilar_ct$drug %in% shared_drugs | !data_I_dissimilar_ct$phenotype %in% shared_phenos, "Status unknown", data_I_dissimilar_ct$max_phase_for_ind)
# NA means that the pair has never been in clinical trials --> assign 0 (we do have CT information about these drugs/phenos and no evidence that have been tested in any phase)
data_I_dissimilar_ct$max_phase_for_ind = ifelse(is.na(data_I_dissimilar_ct$max_phase_for_ind), 0, data_I_dissimilar_ct$max_phase_for_ind) 
# "" means that the pair has been in clinical trials but has unknown phase --> remove it
data_I_dissimilar_ct$max_phase_for_ind = ifelse(data_I_dissimilar_ct$max_phase_for_ind == "", "Unknown phase", data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = gsub("0", "", data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = gsub("1", "Phase 1", data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = gsub("2", "Phase 2", data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = gsub("3", "Phase 3", data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = gsub("4", "Approved", data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = ifelse(data_I_dissimilar_ct$is_indicated == 1, "Approved", data_I_dissimilar_ct$max_phase_for_ind)

data_I_dissimilar_ct = data_I_dissimilar_ct %>% filter(max_phase_for_ind != "Status unknown") # pairs of drugs and phenotypes that don't exist in clinical trials

# baseline probability that a drug treats a phenotype
baseline_prob = prop.table(table(data_I_dissimilar$is_indicated))["1"]
baseline_prob

IND_ORs = data.frame(
  threshold = c(baseline_prob, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
  phase = rep(c("Phase 1","Phase 2","Phase 3","Approved"), each = 8),
  OR = NA,
  OR_ci95_lower = NA,
  OR_ci95_upper = NA,
  nr_predicted_pairs = NA
)
for (i in 1:nrow(IND_ORs)) {
  phase = IND_ORs[i,"phase"]
  threshold = IND_ORs[i,"threshold"]
  temp = data_I_dissimilar_ct %>%
    mutate(predicted = if_else(PredProb > threshold, 1, 2)) %>% # 0: not predicted | 1: predicted
    mutate(predicted = if_else(PredProb < baseline_prob, 0, predicted)) %>% # 0: not predicted | 1: predicted
    filter(predicted %in% c(0,1)) %>%
    filter(max_phase_for_ind %in% c("",phase)) # "" means the drug-phenotype pair has never been in clinical trials
  fisher_or = fisher.test(table(temp$predicted, temp$max_phase_for_ind), alternative = "two.sided")
  IND_ORs[i, "OR"] = fisher_or$estimate
  IND_ORs[i, "OR_ci95_lower"] = fisher_or$conf.int[1]
  IND_ORs[i, "OR_ci95_upper"] = fisher_or$conf.int[2]
  IND_ORs[i, "nr_predicted_pairs"] = sum(temp$predicted == 1)
} ; rm(phase, threshold, fisher_or, temp)
IND_ORs$threshold = round(IND_ORs$threshold, 3)
IND_ORs$threshold = gsub(round(baseline_prob, 3), paste(round(baseline_prob,3), "\n(baseline probability)", sep = ""), IND_ORs$threshold)
IND_ORs$threshold = factor(IND_ORs$threshold, levels = unique(IND_ORs$threshold), labels = paste("> ", unique(IND_ORs$threshold), sep = ""))
IND_ORs$phase = factor(IND_ORs$phase, levels = c("Phase 1","Phase 2","Phase 3","Approved"), labels = c("Phase 1","Phase 2","Phase 3","Approved"))
IND_ORs$facets = ifelse(IND_ORs$phase == "Approved", "Approved", "Phase 1, 2, 3")
IND_ORs$facets = factor(IND_ORs$facets, levels = c("Approved", "Phase 1, 2, 3"), labels = c("Approved", "Phase 1, 2, 3"))
ggplot(IND_ORs, aes(x = threshold, y = OR, group = phase, color = phase, fill = phase)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = nr_predicted_pairs), show.legend = FALSE) +
  geom_ribbon(aes(ymin = OR_ci95_lower, ymax = OR_ci95_upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 1) +
  geom_text_repel(aes(label = round(OR,2), color = phase), data = IND_ORs %>% filter(OR_ci95_lower > 1),
                  nudge_y = 0.2, box.padding = 0.5,
                  size = 7,
                  show.legend = FALSE) +
  labs(
    title = "", 
    x = "Predicted probability",
    y = "Odds ratio",
    color = "",
    fill = ""
  ) +
  scale_color_manual(values = c("Phase 1" = "#1B9E77",
                                "Phase 2" = "#7570B3", 
                                "Phase 3" = "#D95F02", 
                                "Approved" = "#E7298A")) +
  scale_fill_manual(values = c("Phase 1" = "#1B9E77",
                               "Phase 2" = "#7570B3",
                               "Phase 3" = "#D95F02",
                               "Approved" = "#E7298A")) +
  facet_wrap(~facets, nrow = 2, scales = "free_y") +
  theme_bw(base_size = 16, base_family = "Arial") +
  theme(
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(size = 25, family = "Arial", color = "black", margin = margin(r = 20)),
    axis.title.x = element_text(size = 25, family = "Arial", color = "black", margin = margin(t = 20)),
    axis.text = element_text(size = 20, family = "Arial", color = "black"),
    legend.text = element_text(size = 25, family = "Arial", color = "black"),
    legend.position = "top",
    strip.text = element_text(size = 25, family = "Arial", color = "black", face = "bold")
  )
ggsave("figures/Figure_4b.png",
       device = "png", dpi = 600, width = 20, height = 11)

# any clinical trial phase enrichment
IND_ORs_anyphase = data.frame(
  threshold = c(baseline_prob, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
  OR = NA,
  OR_ci95_lower = NA,
  OR_ci95_upper = NA,
  nr_predicted_pairs = NA
)
for (i in 1:nrow(IND_ORs_anyphase)) {
  threshold = IND_ORs_anyphase[i,"threshold"]
  temp = data_I_dissimilar_ct %>%
    mutate(predicted = if_else(PredProb > threshold, 1, 2)) %>% # 0: not predicted | 1: predicted
    mutate(predicted = if_else(PredProb < baseline_prob, 0, predicted)) %>% # 0: not predicted | 1: predicted
    filter(predicted %in% c(0,1)) %>%
    filter(max_phase_for_ind != "Approved") %>%
    mutate(in_ct = if_else(max_phase_for_ind %in% c("Phase 1", "Phase 2", "Phase 3", "Unknown phase"), 1, 0))
  fisher_or = fisher.test(table(temp$predicted, temp$in_ct), alternative = "two.sided")
  IND_ORs_anyphase[i, "OR"] = fisher_or$estimate
  IND_ORs_anyphase[i, "OR_ci95_lower"] = fisher_or$conf.int[1]
  IND_ORs_anyphase[i, "OR_ci95_upper"] = fisher_or$conf.int[2]
  IND_ORs_anyphase[i, "nr_predicted_pairs"] = sum(temp$predicted == 1)
} ; rm(threshold, fisher_or, temp)
IND_ORs_anyphase$threshold = round(IND_ORs_anyphase$threshold, 3)
IND_ORs_anyphase$threshold = gsub(round(baseline_prob, 3), paste(round(baseline_prob,3), "\n(baseline probability)", sep = ""), IND_ORs_anyphase$threshold)
IND_ORs_anyphase$threshold = factor(IND_ORs_anyphase$threshold, levels = unique(IND_ORs_anyphase$threshold), labels = paste("> ", unique(IND_ORs_anyphase$threshold), sep = ""))
ggplot(IND_ORs_anyphase, aes(x = threshold, y = OR)) +
  geom_line(linewidth = 1, group = 1) +
  geom_point(aes(size = nr_predicted_pairs), show.legend = FALSE) +
  geom_ribbon(aes(ymin = OR_ci95_lower, ymax = OR_ci95_upper), alpha = 0.2, group = 1) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 1) +
  geom_text_repel(aes(label = round(OR,2)), data = IND_ORs_anyphase %>% filter(OR_ci95_lower > 1),
                  nudge_y = 0.2, box.padding = 0.5,
                  size = 7,
                  show.legend = FALSE) +
  labs(
    title = "", 
    x = "Predicted probability",
    y = "Odds ratio",
    color = "",
    fill = ""
  ) +
  theme_bw(base_size = 16, base_family = "Arial") +
  theme(
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(size = 25, family = "Arial", color = "black", margin = margin(r = 20)),
    axis.title.x = element_text(size = 25, family = "Arial", color = "black", margin = margin(t = 20)),
    axis.text = element_text(size = 20, family = "Arial", color = "black"),
    legend.text = element_text(size = 25, family = "Arial", color = "black"),
    legend.position = "top",
    strip.text = element_text(size = 25, family = "Arial", color = "black", face = "bold")
  )


# ==================================================================================================================================== #
# Relative success for advancing to the next clinical trial phase
# ==================================================================================================================================== #

RS = data.frame(
  pred_threshold = c(baseline_prob, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25),
  # Phase I --> Phase II
  PhaseI_to_PhaseII_RS = NA,
  PhaseI_to_PhaseII_RS_95ci_lower = NA,
  PhaseI_to_PhaseII_RS_95ci_upper = NA,
  PhaseI_to_PhaseII_RS_pairs_compared = NA,
  
  # Phase II --> Phase III
  PhaseII_to_PhaseIII_RS = NA,
  PhaseII_to_PhaseIII_RS_95ci_lower = NA,
  PhaseII_to_PhaseIII_RS_95ci_upper = NA,
  PhaseII_to_PhaseIII_RS_pairs_compared = NA,
  
  # Phase III --> Approval
  PhaseIII_to_Approval_RS = NA,
  PhaseIII_to_Approval_RS_95ci_lower = NA,
  PhaseIII_to_Approval_RS_95ci_upper = NA,
  PhaseIII_to_Approval_RS_pairs_compared = NA,
  
  # Phase I --> Approval
  PhaseI_to_Approval_RS = NA,
  PhaseI_to_Approval_RS_95ci_lower = NA,
  PhaseI_to_Approval_RS_95ci_upper = NA,
  PhaseI_to_Approval_RS_pairs_compared = NA
)

data_I_dissimilar_ct = left_join(data_I_dissimilar, drug_ct, by = c("drug", "phenotype"))
# for drugs and phenotypes that do not exist in the clinical trials file, we don't know if they have ever been in clinical trials --> set "status unknown"
shared_drugs = intersect(drug_ct$drug, data_I_dissimilar_ct$drug) # we have clinical trials data for 1,327 drugs
shared_phenos = intersect(drug_ct$phenotype, data_I_dissimilar_ct$phenotype) # and 173 phenotypes
data_I_dissimilar_ct$max_phase_for_ind = ifelse(!data_I_dissimilar_ct$drug %in% shared_drugs | !data_I_dissimilar_ct$phenotype %in% shared_phenos, "Status unknown", data_I_dissimilar_ct$max_phase_for_ind)
# NA means that the pair has never been in clinical trials --> assign 0 (we do have CT information about these drugs/phenos and no evidence that have been tested in any phase)
data_I_dissimilar_ct$max_phase_for_ind = ifelse(is.na(data_I_dissimilar_ct$max_phase_for_ind), 0, data_I_dissimilar_ct$max_phase_for_ind) 
# "" means that the pair has been in clinical trials but has unknown phase --> remove it
data_I_dissimilar_ct$max_phase_for_ind = ifelse(data_I_dissimilar_ct$max_phase_for_ind == "", "Unknown phase", data_I_dissimilar_ct$max_phase_for_ind)
## dataset with clinical trials data
data_I_dissimilar_ct$max_phase_for_ind = gsub("Phase 1", 1, data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = gsub("Phase 2", 2, data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = gsub("Phase 3", 3, data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = gsub("Approved", 4, data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct$max_phase_for_ind = gsub("Unknown phase", "", data_I_dissimilar_ct$max_phase_for_ind)
data_I_dissimilar_ct = data_I_dissimilar_ct %>% filter(max_phase_for_ind != "Status unknown") # pairs of drugs and phenotypes that don't exist in clinical trials
## keep pairs with clinical trials data
data_I_dissimilar_ct = data_I_dissimilar_ct[!is.na(data_I_dissimilar_ct$max_phase_for_ind), ]
data_I_dissimilar_ct$max_phase_for_ind = as.numeric(data_I_dissimilar_ct$max_phase_for_ind)

# calculate relative success for each phase advancement by using different percentile thresholds to define low and high predicted probability groups
for (z in 1:nrow(RS)) {
  pred_threshold = as.numeric(RS[z, "pred_threshold"])
  
  ## compare drug-phenotype pairs with pred_prob > threshold to pairs with pred_prob below the baseline frequency of being indication
  data_I_dissimilar_ct = data_I_dissimilar_ct %>%
    mutate(PredProb_cat = if_else(PredProb > pred_threshold, "high", "medium")) %>%
    mutate(PredProb_cat = if_else(PredProb < baseline_prob, "low", PredProb_cat))
  # keep only needed columns
  data_I_dissimilar_ct = data_I_dissimilar_ct %>% 
    dplyr::select(drug, phenotype, is_indicated, max_phase_for_ind, PredProb, PredProb_cat) %>%
    distinct()
  # convert the phase_CT_reached into binary success for each phase
  data_I_ct_long = data_I_dissimilar_ct %>%
    mutate(
      succ_I   = ifelse(max_phase_for_ind >= 2, 1, 0), # Succeed in Phase 1: pairs with max_phase ≥ 2 (reached phase II or higher)
      succ_II  = ifelse(max_phase_for_ind >= 3, 1, 0), # Succeed in Phase 2: pairs with max_phase ≥ 3 (reached phase III or higher)
      succ_III = ifelse(max_phase_for_ind >= 4, 1, 0)  # Succeed in Phase 3: pairs with max_phase ≥ 4 (approved)
    ) %>%
    tidyr::pivot_longer(succ_I:succ_III, names_to = "phase", values_to = "success") %>%
    mutate(phase = gsub("succ_", "", phase))
  
  ## Calculate relative success for phase progression
  # Phase I --> Phase II
  phaseI_to_phaseII = data_I_ct_long %>% 
    # all pairs that have been tested in Phase I
    filter(max_phase_for_ind >= 1) %>%
    # all pairs that have information about success/failure in Phase I
    filter(phase == "I") %>%
    group_by(PredProb_cat, phase) %>%
    dplyr::summarize(
      x = sum(success), 
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      binom = binom.confint(x, n, conf.level = 0.95, method = "wilson")[, c("mean", "lower", "upper")],
      mean = binom$mean,
      l = binom$lower,
      u = binom$upper
    ) %>%
    select(PredProb_cat, phase, x, n, mean, l, u) %>%
    tidyr::pivot_wider(names_from = PredProb_cat, values_from = c(x, n, mean, l, u)) %>%
    mutate(
      RS_high_low = mean_high / mean_low,
      fraction_high_low = glue("{x_high}/{n_high} vs {x_low}/{n_low}")
    )
  
  # Phase II --> Phase III
  phaseII_to_phaseIII = data_I_ct_long %>% 
    # all pairs that have been tested in Phase II
    filter(max_phase_for_ind >= 2) %>%
    # all pairs that have information about success/failure in Phase II
    filter(phase == "II") %>%
    group_by(PredProb_cat, phase) %>%
    dplyr::summarize(
      x = sum(success),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      binom = binom.confint(x, n, conf.level = 0.95, method = "wilson")[, c("mean", "lower", "upper")],
      mean = binom$mean,
      l = binom$lower,
      u = binom$upper
    ) %>%
    select(PredProb_cat, phase, x, n, mean, l, u) %>%
    tidyr::pivot_wider(names_from = PredProb_cat, values_from = c(x, n, mean, l, u)) %>%
    mutate(
      RS_high_low = mean_high / mean_low,
      fraction_high_low = glue("{x_high}/{n_high} vs {x_low}/{n_low}")
    )
  
  # Phase III --> Approval
  phaseIII_to_phaseIV = data_I_ct_long %>% 
    # all pairs that have been tested in Phase III
    filter(max_phase_for_ind >= 3) %>%
    # all pairs that have information about success/failure in Phase III
    filter(phase == "III") %>%
    group_by(PredProb_cat, phase) %>%
    dplyr::summarize(
      x = sum(success),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      binom = binom.confint(x, n, conf.level = 0.95, method = "wilson")[, c("mean", "lower", "upper")],
      mean = binom$mean,
      l = binom$lower,
      u = binom$upper
    ) %>%
    select(PredProb_cat, phase, x, n, mean, l, u) %>%
    tidyr::pivot_wider(names_from = PredProb_cat, values_from = c(x, n, mean, l, u)) %>%
    mutate(
      RS_high_low = mean_high / mean_low,
      fraction_high_low = glue("{x_high}/{n_high} vs {x_low}/{n_low}")
    )
  
  # Phase I --> Approval
  phaseI_to_Approval = data_I_ct_long %>% 
    # all pairs that have been tested in Phase I
    filter(max_phase_for_ind >= 1) %>%
    # all pairs that have information about success/failure in Phase III
    filter(phase == "III") %>%
    group_by(PredProb_cat, phase) %>%
    dplyr::summarize(
      x = sum(success),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      binom = binom.confint(x, n, conf.level = 0.95, method = "wilson")[, c("mean", "lower", "upper")],
      mean = binom$mean,
      l = binom$lower,
      u = binom$upper
    ) %>%
    select(PredProb_cat, phase, x, n, mean, l, u) %>%
    tidyr::pivot_wider(names_from = PredProb_cat, values_from = c(x, n, mean, l, u)) %>%
    mutate(
      RS_high_low = mean_high / mean_low,
      fraction_high_low = glue("{x_high}/{n_high} vs {x_low}/{n_low}")
    )
  
  # combine all together
  temp = rbind(
    phaseI_to_phaseII,
    phaseII_to_phaseIII,
    phaseIII_to_phaseIV,
    phaseI_to_Approval
  )
  temp$phase = c(
    "Phase I to Phase II", 
    "Phase II to Phase III",
    "Phase III to Phase IV",
    "Phase I to Approval"
  )
  
  for (i in 1:nrow(temp)) {
    ## High - Low
    ci_temp = RS_ci95(
      a = as.numeric(temp[i, "x_high"]),
      b = as.numeric(temp[i, "n_high"]) - as.numeric(temp[i, "x_high"]),
      c = as.numeric(temp[i, "x_low"]),
      d = as.numeric(temp[i, "n_low"]) - as.numeric(temp[i, "x_low"])
    )
    temp[i, "RS_high_low_95ci_lower"] = ci_temp[2]
    temp[i, "RS_high_low_95ci_upper"] = ci_temp[3]
  }
  
  # Phase I --> Phase II
  RS[z, "PhaseI_to_PhaseII_RS"] = temp[1, "RS_high_low"]
  RS[z, "PhaseI_to_PhaseII_RS_95ci_lower"] = temp[1, "RS_high_low_95ci_lower"]
  RS[z, "PhaseI_to_PhaseII_RS_95ci_upper"] = temp[1, "RS_high_low_95ci_upper"]
  RS[z, "PhaseI_to_PhaseII_RS_pairs_compared"] = temp[1, "fraction_high_low"]
  # Phase II --> Phase III
  RS[z, "PhaseII_to_PhaseIII_RS"] = temp[2, "RS_high_low"]
  RS[z, "PhaseII_to_PhaseIII_RS_95ci_lower"] = temp[2, "RS_high_low_95ci_lower"]
  RS[z, "PhaseII_to_PhaseIII_RS_95ci_upper"] = temp[2, "RS_high_low_95ci_upper"]
  RS[z, "PhaseII_to_PhaseIII_RS_pairs_compared"] = temp[2, "fraction_high_low"]
  # Phase III --> Approval
  RS[z, "PhaseIII_to_Approval_RS"] = temp[3, "RS_high_low"]
  RS[z, "PhaseIII_to_Approval_RS_95ci_lower"] = temp[3, "RS_high_low_95ci_lower"]
  RS[z, "PhaseIII_to_Approval_RS_95ci_upper"] = temp[3, "RS_high_low_95ci_upper"]
  RS[z, "PhaseIII_to_Approval_RS_pairs_compared"] = temp[3, "fraction_high_low"]
  # Phase I --> Approval
  RS[z, "PhaseI_to_Approval_RS"] = temp[4, "RS_high_low"]
  RS[z, "PhaseI_to_Approval_RS_95ci_lower"] = temp[4, "RS_high_low_95ci_lower"]
  RS[z, "PhaseI_to_Approval_RS_95ci_upper"] = temp[4, "RS_high_low_95ci_upper"]
  RS[z, "PhaseI_to_Approval_RS_pairs_compared"] = temp[4, "fraction_high_low"]
}
rm(data_I_ct_long, phaseI_to_Approval, phaseI_to_phaseII, phaseII_to_phaseIII, phaseIII_to_phaseIV, ci_temp, i, z, pred_threshold, temp)
# set x-axis names
RS$x_axis_names = paste(">", round(RS$pred_threshold, 4))
RS$x_axis_names = gsub(paste(">", round(baseline_prob, 4)), paste("> ", round(baseline_prob,3), "\n(baseline probability)", sep=""), RS$x_axis_names)
RS$x_axis_names = factor(RS$x_axis_names, levels = RS$x_axis_names, labels = RS$x_axis_names)

# Phase I --> Phase II
ggplot(RS, aes(x = x_axis_names, y = PhaseI_to_PhaseII_RS, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = PhaseI_to_PhaseII_RS_95ci_lower, ymax = PhaseI_to_PhaseII_RS_95ci_upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 1) +
  # text showing relative success above points
  geom_text(
    data = RS %>% filter(PhaseI_to_PhaseII_RS_95ci_lower > 1),
    aes(label = round(PhaseI_to_PhaseII_RS, 2)),
    vjust = -1, size = 7, show.legend = FALSE, fontface = "bold") +
  # text showing pairs compared below points
  geom_label_repel(aes(label = PhaseI_to_PhaseII_RS_pairs_compared), direction = "y",
                   vjust = 0, size = 4, show.legend = FALSE) +
  theme_bw(base_size = 14, base_family = "Arial") +
  labs(title = "Phase I → Phase II",
       x = "\nPredicted probability",
       y = "Relative success\n") +
  theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 18, color = "black"),
    axis.title = element_text(size = 20, family = "Arial", color = "black"),
    axis.text = element_text(size = 16, family = "Arial", color = "black")
  )
ggsave("figures/Supplementary_figure_8a.png",
       device = "png", dpi = 600, width = 12, height = 6)

# Phase II --> Phase III
ggplot(RS, aes(x = x_axis_names, y = PhaseII_to_PhaseIII_RS, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = PhaseII_to_PhaseIII_RS_95ci_lower, ymax = PhaseII_to_PhaseIII_RS_95ci_upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 1) +
  # text showing relative success above points
  geom_text(
    data = RS %>% filter(PhaseII_to_PhaseIII_RS_95ci_lower > 1),
    aes(label = round(PhaseII_to_PhaseIII_RS, 2)),
    vjust = -1, size = 7, show.legend = FALSE, fontface = "bold") +
  # text showing pairs compared below points
  geom_label_repel(aes(label = PhaseII_to_PhaseIII_RS_pairs_compared), direction = "y",
                   vjust = 0, size = 4, show.legend = FALSE) +
  theme_bw(base_size = 14, base_family = "Arial") +
  labs(title = "Phase II → Phase III",
       x = "\nPredicted probability",
       y = "Relative success\n") +
  theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 18, color = "black"),
    axis.title = element_text(size = 20, family = "Arial", color = "black"),
    axis.text = element_text(size = 16, family = "Arial", color = "black")
  )
ggsave("figures/Supplementary_figure_8b.png",
       device = "png", dpi = 600, width = 12, height = 6)

# Phase III --> Approval
ggplot(RS, aes(x = x_axis_names, y = PhaseIII_to_Approval_RS, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = PhaseIII_to_Approval_RS_95ci_lower, ymax = PhaseIII_to_Approval_RS_95ci_upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 1) +
  # text showing relative success above points
  geom_text(
    data = RS %>% filter(PhaseIII_to_Approval_RS_95ci_lower > 1),
    aes(label = round(PhaseIII_to_Approval_RS, 2)),
    vjust = -1, size = 7, show.legend = FALSE, fontface = "bold") +
  # text showing pairs compared below points
  geom_label_repel(aes(label = PhaseIII_to_Approval_RS_pairs_compared), direction = "y",
                   vjust = 0, size = 4, show.legend = FALSE) +
  theme_bw(base_size = 14, base_family = "Arial") +
  labs(title = "Phase III → Approval",
       x = "\nPredicted probability",
       y = "Relative success\n") +
  theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 18, color = "black"),
    axis.title = element_text(size = 20, family = "Arial", color = "black"),
    axis.text = element_text(size = 16, family = "Arial", color = "black")
  )
ggsave("figures/Supplementary_figure_8c.png",
       device = "png", dpi = 600, width = 12, height = 6)

# Phase I --> Approval
ggplot(RS, aes(x = x_axis_names, y = PhaseI_to_Approval_RS, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = PhaseI_to_Approval_RS_95ci_lower, ymax = PhaseI_to_Approval_RS_95ci_upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 1) +
  # text showing relative success above points
  geom_text(
    data = RS %>% filter(PhaseI_to_Approval_RS_95ci_lower > 1),
    aes(label = round(PhaseI_to_Approval_RS, 2)),
    vjust = -1, size = 7, show.legend = FALSE, fontface = "bold") +
  # text showing pairs compared below points
  geom_label_repel(aes(label = PhaseI_to_Approval_RS_pairs_compared), direction = "y",
                   vjust = 0, size = 4, show.legend = FALSE) +
  theme_bw(base_size = 14, base_family = "Arial") +
  labs(title = "Phase I → Approval",
       x = "\nPredicted probability",
       y = "Relative success\n") +
  theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 18, color = "black"),
    axis.title = element_text(size = 20, family = "Arial", color = "black"),
    axis.text = element_text(size = 16, family = "Arial", color = "black")
  )
ggsave("figures/Figure_4c.png",
       device = "png", dpi = 600, width = 16, height = 6)


#################################################################################
##### compare predicted probabilities between "same" and "different" models #####

ultimate = full_join(
  data_I_similar,
  data_I_dissimilar,
  by = c("drug", "phenotype", "phenotype_name", "phenotype_category", "is_indicated")
)
ultimate = na.omit(ultimate)
rownames(ultimate) = NULL
colnames(ultimate)[c(5,7)] = c("PredProb_similar", "PredProb_dissimilar")

cor.test(ultimate$PredProb_similar, ultimate$PredProb_dissimilar, method = "pearson") # 0.15
cor.test(ultimate$PredProb_similar, ultimate$PredProb_dissimilar, method = "pearson")$p.value # 2.66-153

cor.test(ultimate$PredProb_similar, ultimate$PredProb_dissimilar, method = "spearman") # 0.08
cor.test(ultimate$PredProb_similar, ultimate$PredProb_dissimilar, method = "spearman")$p.value # 2.78e-45


roc(ultimate$is_indicated ~ ultimate$PredProb_similar) # 0.6744 - similar
roc(ultimate$is_indicated ~ ultimate$PredProb_dissimilar) # 0.5681 - dissimilar
