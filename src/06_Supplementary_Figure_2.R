
#### Purpose: to test whether phenotype pairs sharing more genetics also share more drugs
## Using disease labels from Open Targets, ICD10 and UK HRCS to group diseases into same and different body systems

library(data.table)
library(dplyr)
library(DescTools)
library(ggplot2)
source("src/00_functions.R")

## ============ ##
## Load data    ##
## ============ ##

data_AllPhenoPairs = fread("data/AllPhenoPairs_DrugOverlap_GeneticSimilarity.txt", data.table = FALSE)

### add phenotype annotations
# UK HRCS & ICD10
trait_categories = fread("data/phenotype_categories.txt") %>%
  dplyr::select(mesh_id, trait_area_HRCS, trait_ICD10_top)
data_AllPhenoPairs = data_AllPhenoPairs %>%
  left_join(trait_categories, by = c("phenotype_1" = "mesh_id")) %>%
  left_join(trait_categories, by = c("phenotype_2" = "mesh_id")) %>%
  dplyr::select(phenotype_1, phenotype_2, drug_se_overlap_coef, drug_ind_overlap_coef, ldsc_score, magma_cor, smultixcan_cor, l2g_cor, coloc_overcoef,
                phenotype_1_category_OT = phenotype_1_category, phenotype_2_category_OT = phenotype_2_category,
                phenotype_1_category_HRCS = trait_area_HRCS.x, phenotype_2_category_HRCS = trait_area_HRCS.y,
                phenotype_1_category_ICD10 = trait_ICD10_top.x, phenotype_2_category_ICD10 = trait_ICD10_top.y) %>%
  mutate(
    same_body_system_consensus_OT_HRCS_ICD10 = if_else(
      phenotype_1_category_OT == phenotype_2_category_OT | phenotype_1_category_HRCS == phenotype_2_category_HRCS | phenotype_1_category_ICD10 == phenotype_2_category_ICD10,
      "Same", "Different")) %>%
  dplyr::select(phenotype_1, phenotype_2, same_body_system_consensus_OT_HRCS_ICD10, drug_se_overlap_coef:coloc_overcoef)

## using disease annotations only (OT, UK HRCS, ICD10)
data_AllPhenoPairs$phenotype_category = factor(
  data_AllPhenoPairs$same_body_system_consensus_OT_HRCS_ICD10,
  levels = c("Same","Different"),
  labels = c("Same","Different")
)

# ============================================================= #
# Disease pairs that share more genetics also share more drugs  #
# ============================================================= #

## create data frame with spearman correlation values and 95% confidence intervals between genetic similarity and drug overlap metrics
corr_DrugSharing_GenSim = data.frame(
  feature = rep(c("LDSC", "MAGMA", "S-MultiXcan", "L2G", "COLOC"), each = 6),
  drug_effect = rep(c("Indications", "Side effects"), each = 3),
  body_system = c("All", "Same", "Different"),
  corr_spearman = c(
    # ldsc - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "Different", conf_int = 0.95)[1],
    # ldsc - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "Different", conf_int = 0.95)[1],
    # magma - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "Different", conf_int = 0.95)[1],
    # magma - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "Different", conf_int = 0.95)[1],
    # smultixcan - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "Different", conf_int = 0.95)[1],
    # smultixcan - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "Different", conf_int = 0.95)[1],
    # l2g - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "Different", conf_int = 0.95)[1],
    # l2g - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "Different", conf_int = 0.95)[1],
    # coloc - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "Different", conf_int = 0.95)[1],
    # coloc - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "All", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "Same", conf_int = 0.95)[1],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "Different", conf_int = 0.95)[1]
  ),
  corr_pvalue = c(
    # ldsc - indication
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "Different"),
    # ldsc - side effect
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "Different"),
    # magma - indication
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "Different"),
    # magma - side effect
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "Different"),
    # smultixcan - indication
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "Different"),
    # smultixcan - side effect
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "Different"),
    # l2g - indication
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "Different"),
    # l2g - side effect
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "Different"),
    # coloc - indication
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "Different"),
    # coloc - side effect
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "All"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "Same"),
    spearman_cor_pvalue(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "Different")
  ),
  corr_95ci_lower = c(
    # ldsc - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "Different", conf_int = 0.95)[2],
    # ldsc - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "Different", conf_int = 0.95)[2],
    # magma - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "Different", conf_int = 0.95)[2],
    # magma - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "Different", conf_int = 0.95)[2],
    # smultixcan - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "Different", conf_int = 0.95)[2],
    # smultixcan - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "Different", conf_int = 0.95)[2],
    # l2g - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "Different", conf_int = 0.95)[2],
    # l2g - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "Different", conf_int = 0.95)[2],
    # coloc - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "Different", conf_int = 0.95)[2],
    # coloc - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "All", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "Same", conf_int = 0.95)[2],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "Different", conf_int = 0.95)[2]
  ),
  corr_95ci_upper = c(
    # ldsc - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "ldsc", body_system = "Different", conf_int = 0.95)[3],
    # ldsc - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "ldsc", body_system = "Different", conf_int = 0.95)[3],
    # magma - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "magma", body_system = "Different", conf_int = 0.95)[3],
    # magma - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "magma", body_system = "Different", conf_int = 0.95)[3],
    # smultixcan - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "smultixcan", body_system = "Different", conf_int = 0.95)[3],
    # smultixcan - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "smultixcan", body_system = "Different", conf_int = 0.95)[3],
    # l2g - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "l2g", body_system = "Different", conf_int = 0.95)[3],
    # l2g - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "l2g", body_system = "Different", conf_int = 0.95)[3],
    # coloc - indication
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "indication", metric = "coloc", body_system = "Different", conf_int = 0.95)[3],
    # coloc - side effect
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "All", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "Same", conf_int = 0.95)[3],
    spearman_cor_ci(data = data_AllPhenoPairs, drug_effect = "side effect", metric = "coloc", body_system = "Different", conf_int = 0.95)[3]
  )
)

## Add number of phenotype pairs tested in each case
corr_DrugSharing_GenSim[1, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "ldsc", drug_effect = "indication")
corr_DrugSharing_GenSim[2, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "ldsc", drug_effect = "indication")
corr_DrugSharing_GenSim[3, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "ldsc", drug_effect = "indication")
corr_DrugSharing_GenSim[4, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "ldsc", drug_effect = "side effect")
corr_DrugSharing_GenSim[5, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "ldsc", drug_effect = "side effect")
corr_DrugSharing_GenSim[6, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "ldsc", drug_effect = "side effect")

corr_DrugSharing_GenSim[7, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "magma", drug_effect = "indication")
corr_DrugSharing_GenSim[8, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "magma", drug_effect = "indication")
corr_DrugSharing_GenSim[9, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "magma", drug_effect = "indication")
corr_DrugSharing_GenSim[10, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "magma", drug_effect = "side effect")
corr_DrugSharing_GenSim[11, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "magma", drug_effect = "side effect")
corr_DrugSharing_GenSim[12, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "magma", drug_effect = "side effect")

corr_DrugSharing_GenSim[13, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "smultixcan", drug_effect = "indication")
corr_DrugSharing_GenSim[14, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "smultixcan", drug_effect = "indication")
corr_DrugSharing_GenSim[15, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "smultixcan", drug_effect = "indication")
corr_DrugSharing_GenSim[16, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "smultixcan", drug_effect = "side effect")
corr_DrugSharing_GenSim[17, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "smultixcan", drug_effect = "side effect")
corr_DrugSharing_GenSim[18, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "smultixcan", drug_effect = "side effect")

corr_DrugSharing_GenSim[19, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "l2g", drug_effect = "indication")
corr_DrugSharing_GenSim[20, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "l2g", drug_effect = "indication")
corr_DrugSharing_GenSim[21, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "l2g", drug_effect = "indication")
corr_DrugSharing_GenSim[22, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "l2g", drug_effect = "side effect")
corr_DrugSharing_GenSim[23, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "l2g", drug_effect = "side effect")
corr_DrugSharing_GenSim[24, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "l2g", drug_effect = "side effect")

corr_DrugSharing_GenSim[25, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "coloc", drug_effect = "indication")
corr_DrugSharing_GenSim[26, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "coloc", drug_effect = "indication")
corr_DrugSharing_GenSim[27, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "coloc", drug_effect = "indication")
corr_DrugSharing_GenSim[28, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "All", metric = "coloc", drug_effect = "side effect")
corr_DrugSharing_GenSim[29, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Same", metric = "coloc", drug_effect = "side effect")
corr_DrugSharing_GenSim[30, "nr_pheno_pairs"] = nr_pheno_pairs_tested(data = data_AllPhenoPairs, body_system = "Different", metric = "coloc", drug_effect = "side effect")

## create factors
corr_DrugSharing_GenSim$feature = factor(
  corr_DrugSharing_GenSim$feature, 
  levels = corr_DrugSharing_GenSim$feature, 
  labels = corr_DrugSharing_GenSim$feature
)
corr_DrugSharing_GenSim$body_system = factor(
  corr_DrugSharing_GenSim$body_system, 
  levels = corr_DrugSharing_GenSim$body_system, 
  labels = corr_DrugSharing_GenSim$body_system
)
corr_DrugSharing_GenSim$drug_effect = factor(
  corr_DrugSharing_GenSim$drug_effect, 
  levels = corr_DrugSharing_GenSim$drug_effect,
  labels = corr_DrugSharing_GenSim$drug_effect
)

## Plot
# shaded areas to separate the genetic similarity metrics on the y-axis
bg_df = corr_DrugSharing_GenSim %>%
  distinct(feature) %>%
  arrange(feature) %>%
  mutate(y_id = as.numeric(feature),
         ymin = y_id - 0.5,
         ymax = y_id + 0.5,
         fill = rep(c(NA, "gray88"), length.out = n()))
ggplot(corr_DrugSharing_GenSim, aes(x = corr_spearman, y = feature, color = body_system)) +
  geom_rect(data = bg_df, inherit.aes = FALSE,
            aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = fill),
            alpha = 0.3) +
  scale_fill_identity() +
  geom_errorbar(aes(xmin = corr_95ci_lower, xmax = corr_95ci_upper), linewidth = 1.2, width = 0, position = position_dodge(width = -0.6)) +
  geom_point(position = position_dodge(width = -0.6), aes(size = nr_pheno_pairs), stroke = 0.2) +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed", color = "grey40") +
  labs(color = "Body system", size = "# tested disease pairs") +
  xlab("Genetic similarity ~ Drug sharing\n(spearman correlation)") +
  ylab("") +
  scale_x_continuous(breaks = seq(-0.15, 0.6, 0.15), limits=c(-0.15, 0.6)) +
  scale_color_viridis_d(option = "C", end = 0.9) +
  scale_size(range = c(3, 9)) +  # more visual separation in sample sizes
  scale_y_discrete(limits=rev) +
  facet_wrap(~drug_effect, ncol = 2) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 30, color = "black"),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.title.x = element_text(size = 30, color = "black", margin = margin(t = 20)),
    strip.text = element_text(size = 30, face = "bold", color = "black"),
    strip.background = element_rect(color = "black", fill = "grey90"),
    legend.title = element_text(size = 30, face = "bold", color = "black"),
    legend.text = element_text(size = 25, color = "black"),
    panel.spacing = unit(1, "cm")
  ) + 
  guides(color = guide_legend(override.aes = list(size = 4, linewidth = 1)))
ggsave("figures/Supplementary_figure_2.png",
       device = "png", dpi = 600, width = 22, height = 10)
