
#### Purpose: to test whether disease pairs sharing more genetics also share more drugs (across different ClinGraph phenotypic similarity thresholds)

library(data.table)
library(dplyr)
library(DescTools)
library(ggplot2)
library(viridisLite)
source("src/00_functions.R")

## ============ ##
## Load data    ##
## ============ ##

data_AllPhenoPairs = fread("data/AllPhenoPairs_DrugOverlap_GeneticSimilarity.txt", data.table = FALSE)

### add clingraph phenotypic similarity (ClinGraph)
clingraph_cosine_df = fread("data/ClinGraph_cosine_similarity.txt")
data_AllPhenoPairs = data_AllPhenoPairs %>%
  left_join(clingraph_cosine_df, by = c("phenotype_1" = "id1_mesh", "phenotype_2" = "id2_mesh")) %>%
  dplyr::select(phenotype_1, phenotype_2, drug_se_overlap_coef, drug_ind_overlap_coef, ldsc_score, magma_cor, smultixcan_cor, l2g_cor, coloc_overcoef,
                clingraph_sim = cosine_similarity_mean) %>%
  dplyr::select(phenotype_1, phenotype_2, clingraph_sim, drug_se_overlap_coef:coloc_overcoef)

## create data frame with spearman correlation values and 95% confidence intervals between genetic similarity and drug overlap metrics
## across varying ClinGraph phenotypic similarity thresholds
corr_ultimate = data.frame()

for (threshold in c(-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
  data_AllPhenoPairs$phenotype_category = ifelse(
    data_AllPhenoPairs$clingraph_sim >= threshold,
    "Same", "Different"
  )
  
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
  
  corr_DrugSharing_GenSim$clingraph_threshold = threshold
  corr_ultimate = rbind(corr_ultimate, corr_DrugSharing_GenSim)
  cat(threshold, "\n")
}

## create factors
corr_ultimate$feature = factor(
  corr_ultimate$feature, 
  levels = corr_ultimate$feature, 
  labels = corr_ultimate$feature
)
corr_ultimate$body_system = factor(
  corr_ultimate$body_system, 
  levels = c("All", "Same", "Different"), 
  labels = c("All", "≥ threshold", "< threshold")
)
corr_ultimate$drug_effect = factor(
  corr_ultimate$drug_effect, 
  levels = corr_ultimate$drug_effect,
  labels = corr_ultimate$drug_effect
)

## Plot
corr_ultimate_all = corr_ultimate %>% filter(body_system == "All")
corr_ultimate_ind = corr_ultimate %>% filter(drug_effect == "Indications")
corr_ultimate = corr_ultimate %>% 
  group_by(feature, drug_effect) %>%
  mutate(hline = unique(corr_spearman[body_system == "All"])) %>%
  ungroup()

# number of body systems (excluding "All")
n_bs = length(unique(corr_ultimate$body_system[corr_ultimate$body_system != "All"]))
vline_color = viridis(n_bs + 1, option = "C", end = 0.9)[1]
vir_cols = viridis(n_bs + 1, option = "C", end = 0.9)[-1]

ggplot(
  corr_ultimate %>% filter(body_system != "All"),
  aes(x = clingraph_threshold, y = corr_spearman, color = body_system)
) +
  geom_ribbon(
    aes(ymin = corr_95ci_lower, ymax = corr_95ci_upper, fill = body_system),
    alpha = 0.4, show.legend = FALSE
  ) +
  geom_line() +
  geom_point(aes(size = nr_pheno_pairs)) +
  scale_size_area(
    max_size = 5,
    name = "# disease pairs"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", color = "grey40") +
  geom_hline(aes(yintercept = hline), linewidth = 0.6, color = vline_color) +
  scale_x_continuous(breaks = seq(-1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
  scale_color_manual(values = vir_cols) +
  scale_fill_manual(values = vir_cols) +
  labs(
    color = "ClinGraph-based\nphenotype similarity",
    x = "Phenotypic similarity threshold\n(ClinGraph cosine similarity)",
    caption = "*Pairs with ClinGraph similarity above the threshold are considered phenotypically similar",
    y = "Genetic similarity ~ Drug sharing\n(spearman correlation)"
  ) +
  ggh4x::facet_grid2(
    drug_effect ~ feature, 
    scales = "free_y", 
    independent = "y"
  ) +
  
  theme_bw(base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 15, margin = margin(t = 20), color = "black"), 
    plot.caption = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    strip.background = element_rect(color = "black", fill = "grey90"),
    legend.title = element_text(size = 15, face = "bold", color = "black"),
    legend.text = element_text(size = 12),
    panel.spacing = unit(1, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, linewidth = 1)))
ggsave("figures/Figure_2.png",
       device = "png", dpi = 600, width = 24, height = 10)