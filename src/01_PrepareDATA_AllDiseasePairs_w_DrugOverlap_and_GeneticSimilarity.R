#### Purpose: to create a data frame of all phenotype pairs with genetic similarity and drug overlap 

library(data.table)
library(dplyr)
library(tidyr)
source("src/00_functions.R")

## ======================== ##
## Load and process data    ##
## ======================== ##

#### GWAS ID - MeSH ID mapper (links GWAS IDs to MeSH disease IDs) ####
gwas_matching_sources = fread("data/gwas_matching_between_sources_manual_final_2.txt", na.strings = "")

#### Drug indications (ChEMBL, RxNORM, SIDER) ####
drug_indications = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/drugs_indications/drug_indications_chembl_rxnorm_sider_combined_extended_chemblid_mesh.txt")
drug_indications = drug_indications %>%
  filter(max_phase_for_ind == 4) %>% # keep only approved drug indications
  distinct() %>%
  filter(mesh_index == "parent") %>% # keep only originally mentioned phenotype
  filter(mesh_id %in% gwas_matching_sources$mesh_id) %>% # keep phenotypes with GWAS data available
  dplyr::select(mesh_id, chembl_id) %>%
  distinct()
# for each phenotype (MeSH ID), create a set of indicated drugs
mesh_chembl_sets = drug_indications %>%
  group_by(mesh_id) %>%
  summarise(chembl_set = list(chembl_id), .groups = "drop")
# create all phenotype pairs
mesh_pairs = expand.grid(
  mesh_id1 = mesh_chembl_sets$mesh_id,
  mesh_id2 = mesh_chembl_sets$mesh_id) %>%
  filter(mesh_id1 != mesh_id2)
# compute overlap coefficient of drug indications for each pair
drug_ind_overlap = mesh_pairs %>%
  rowwise() %>%
  mutate(
    overlap = {
      set1 = mesh_chembl_sets$chembl_set[mesh_chembl_sets$mesh_id == mesh_id1][[1]]
      set2 = mesh_chembl_sets$chembl_set[mesh_chembl_sets$mesh_id == mesh_id2][[1]]
      overlap_coeff(set1, set2)
    }) %>%
  ungroup()

#### Drug side effects (onSIDES) ####
drug_side_effects = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/nSIDES/onSIDES_processed.txt")
drug_side_effects = drug_side_effects %>%
  filter(mesh_id %in% gwas_matching_sources$mesh_id) %>% # keep GWAS phenotypes only
  dplyr::select(mesh_id, chembl_id) %>%
  distinct()
# for each phenotype (MeSH ID), compute set of drugs causing side effects
mesh_chembl_sets = drug_side_effects %>%
  group_by(mesh_id) %>%
  summarise(chembl_set = list(chembl_id), .groups = "drop")
# create all phenotype pairs
mesh_pairs = expand.grid(
  mesh_id1 = mesh_chembl_sets$mesh_id,
  mesh_id2 = mesh_chembl_sets$mesh_id) %>%
  filter(mesh_id1 != mesh_id2)
# compute overlap coefficient of drug side effects for each pair
drug_se_overlap = mesh_pairs %>%
  rowwise() %>%
  mutate(overlap = {
    set1 = mesh_chembl_sets$chembl_set[mesh_chembl_sets$mesh_id == mesh_id1][[1]]
    set2 = mesh_chembl_sets$chembl_set[mesh_chembl_sets$mesh_id == mesh_id2][[1]]
    overlap_coeff(set1, set2)}
  ) %>%
  ungroup()

#### LDSC-based metric ####
ldsc_scores = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/gwasatlas_watanabe/raw/genetic_correlation_v20191115.txt")
ldsc_scores = ldsc_scores %>%
  filter(id1 %in% gwas_matching_sources$ldsc_magma & id2 %in% gwas_matching_sources$ldsc_magma) %>% # keep only GWAS in sample
  dplyr::select(id1, id2, rg) %>% 
  distinct() %>%
  left_join(gwas_matching_sources[,c("ldsc_magma", "mesh_id")], by = c("id1" = "ldsc_magma")) %>%  # add MeSH IDs for both GWAS
  dplyr::rename("id1_mesh_id" = "mesh_id") %>%
  left_join(gwas_matching_sources[,c("ldsc_magma", "mesh_id")], by = c("id2" = "ldsc_magma")) %>%
  dplyr::rename("id2_mesh_id" = "mesh_id")

#### MAGMA-based metric ####
magma_pvalues = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/gwasatlas_watanabe/preprocessed/magma_v20191115_entrezIDs.txt")
magma_pvalues = magma_pvalues %>% 
  dplyr::select(-entrezID) %>%
  distinct()
# keep only GWAS present in our sample
columns_to_keep = c(1, which(colnames(magma_pvalues) %in% gwas_matching_sources$ldsc_magma)) 
magma_pvalues = magma_pvalues[, ..columns_to_keep] ; rm(columns_to_keep)
# remove 10 phenotypes with most missing values
magma_pvalues_nr_na_per_column = data.frame(nr_na = unlist(apply(magma_pvalues, MARGIN = 2, FUN = function(x) {sum(is.na(x))})))
magma_pvalues_nr_na_per_column = magma_pvalues_nr_na_per_column %>% arrange(desc(nr_na))
gwas_to_keep = setdiff(rownames(magma_pvalues_nr_na_per_column), rownames(magma_pvalues_nr_na_per_column)[1:10])
magma_pvalues = magma_pvalues[, ..gwas_to_keep]
# keep only genes with no missing data
magma_pvalues = magma_pvalues[complete.cases(magma_pvalues), ]
magma_pvalues = magma_pvalues %>% dplyr::select(GENE, everything())
rm(magma_pvalues_nr_na_per_column, gwas_to_keep)
# keep only genes significant in ≥1 phenotype (BH<0.05)
magma_pvalues_bh = as.data.frame(apply(magma_pvalues[, -1], MARGIN = 2, FUN = function(x) {p.adjust(x, method = "BH")}))
magma_pvalues_bh = cbind(magma_pvalues[, 1], magma_pvalues_bh)
sig_genes = sort(unique(as.numeric(unlist(apply(magma_pvalues_bh[,-1], MARGIN = 2, FUN = function(x) {which(x < 0.05)})))))
magma_pvalues = magma_pvalues[sig_genes, ]
rm(magma_pvalues_bh, sig_genes)
# compute pairwise spearman correlation across phenotypes
magma_pvalues_spearman_cor_all = cor(magma_pvalues[,-1], method = "spearman")
magma_pvalues_spearman_cor_all = reshape2::melt(magma_pvalues_spearman_cor_all, varnames = c("Disease1", "Disease2"), value.name = "magma_spearman_correlation")
# map disease ID to MeSH
magma_pvalues_spearman_cor_all = magma_pvalues_spearman_cor_all %>%
  left_join(gwas_matching_sources[, c("ldsc_magma", "mesh_id")], by = c("Disease1" = "ldsc_magma")) %>%
  dplyr::rename(Disease1_mesh_id = mesh_id) %>%
  left_join(gwas_matching_sources[, c("ldsc_magma", "mesh_id")], by = c("Disease2" = "ldsc_magma")) %>%
  dplyr::rename(Disease2_mesh_id = mesh_id) %>%
  dplyr::select(Disease1_mesh_id, Disease2_mesh_id, magma_spearman_correlation) %>%
  distinct()

#### S-MultiXcan-based metric ####
smultixcan_values = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/smultixcan/multixcan_all.txt")
# keep only GWAS in sample
columns_to_keep = c(1, which(colnames(smultixcan_values) %in% gwas_matching_sources$smultixcan))
smultixcan_values = smultixcan_values[, ..columns_to_keep] ; rm(columns_to_keep)
smultixcan_values = as.data.frame(smultixcan_values)
# compute pairwise spearman correlation across phenotypes
smultixcan_spearman_cor_all = cor(smultixcan_values[,-1], method = "spearman")
smultixcan_spearman_cor_all = reshape2::melt(smultixcan_spearman_cor_all, varnames = c("Disease1", "Disease2"), value.name = "smultixcan_spearman_cor")
# map disease ID to MeSH
smultixcan_spearman_cor_all = smultixcan_spearman_cor_all %>%
  left_join(gwas_matching_sources[, c("smultixcan", "mesh_id")], by = c("Disease1" = "smultixcan")) %>%
  dplyr::rename(Disease1_mesh_id = mesh_id) %>%
  left_join(gwas_matching_sources[, c("smultixcan", "mesh_id")], by = c("Disease2" = "smultixcan")) %>%
  dplyr::rename(Disease2_mesh_id = mesh_id) %>%
  dplyr::select(Disease1_mesh_id, Disease2_mesh_id, smultixcan_spearman_cor) %>%
  distinct()

#### L2G-based metric ####
ot_gwas_l2g = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/opentargets/processed/l2g.txt")
# keep only GWAS in sample
ot_gwas_l2g = ot_gwas_l2g %>% 
  filter(study_id %in% gwas_matching_sources$opentargets) %>%
  distinct()
# collapse duplicated gene–GWAS pairs by taking max L2G score
ot_gwas_l2g = ot_gwas_l2g %>% 
  dplyr::select(study_id, gene_id, l2g_score = y_proba_full_model) %>%
  arrange(study_id) %>%
  group_by(study_id, gene_id) %>% 
  mutate(l2g_score = max(l2g_score)) %>%
  ungroup() %>%
  distinct()
# all genes with L2G score for at least one phenotype
ot_gwas_l2g_unique_genes = unique(ot_gwas_l2g$gene_id) 
# build gene × disease matrix with L2G scores
temp = data.frame(
  study_id = rep(
    unique(ot_gwas_l2g$study_id), 
    each = length(ot_gwas_l2g_unique_genes)), 
  gene = ot_gwas_l2g_unique_genes
)
temp = left_join(temp, ot_gwas_l2g, by = c("study_id", "gene" = "gene_id"))
ot_gwas_l2g = temp ; rm(temp)
# replace missing scores with 0 (no evidence)
ot_gwas_l2g$l2g_score = ifelse(is.na(ot_gwas_l2g$l2g_score), 0, ot_gwas_l2g$l2g_score)
ot_gwas_l2g = as.data.frame(pivot_wider(ot_gwas_l2g, names_from = study_id, values_from = l2g_score))
# calculate pairwise spearman correlation across phenotypes
ot_l2g_spearman_cor_all = cor(ot_gwas_l2g[,-1], method = "spearman")
ot_l2g_spearman_cor_all = reshape2::melt(ot_l2g_spearman_cor_all, varnames = c("Disease1", "Disease2"), value.name = "ot_l2g_spearman_cor")
# map disease ID to MeSH
ot_l2g_spearman_cor_all = ot_l2g_spearman_cor_all %>%
  left_join(gwas_matching_sources[, c("opentargets", "mesh_id")], by = c("Disease1" = "opentargets")) %>%
  dplyr::rename(Disease1_mesh_id = mesh_id) %>%
  left_join(gwas_matching_sources[, c("opentargets", "mesh_id")], by = c("Disease2" = "opentargets")) %>%
  dplyr::rename(Disease2_mesh_id = mesh_id) %>%
  dplyr::select(Disease1_mesh_id, Disease2_mesh_id, ot_l2g_spearman_cor) %>%
  distinct()

#### COLOC-based metric ####
ot_gwas_coloc = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/opentargets/processed/coloc.txt")
ot_gwas_coloc = ot_gwas_coloc %>% 
  filter(left_study %in% gwas_matching_sources$opentargets) %>% # keep only GWAS in sample
  filter(coloc_h4 >= 0.8) %>% # filter for colocalizations with H4>=0.8
  filter(right_type != "gwas") %>% # filter for colocalization between GWAS and mQTLs
  dplyr::select(left_study, right_gene_id) %>%
  distinct() %>%
  group_by(left_study) %>% # create a list of genes for each phenotype/GWAS
  summarise(genes = list(right_gene_id))
# calculate pairwise overlap coefficient across phenotypes using the disease-genes supported by colocalization with QTLs
ot_gwas_coloc_overlap_coeff = pairwise_overlap(ot_gwas_coloc, overlap_coeff)
# map disease ID to MeSH
ot_gwas_coloc_overlap_coeff = ot_gwas_coloc_overlap_coeff %>%
  left_join(gwas_matching_sources[, c("mesh_id", "opentargets")], by = c("study1" = "opentargets")) %>%
  dplyr::rename(study1_mesh_id = mesh_id) %>%
  left_join(gwas_matching_sources[, c("mesh_id", "opentargets")], by = c("study2" = "opentargets")) %>% 
  dplyr::rename(study2_mesh_id = mesh_id) %>% 
  dplyr::select(study1_mesh_id, study2_mesh_id, ot_coloc_overlap_coefficient = overlap) %>%
  distinct()


## ================================================================================ ##
## Merge all genetic similarity and drug overlap metrics into a single dataframe
## ================================================================================ ##

## create dataframe with all phenotype pairs
pheno_pairs_all_data = t(combn(gwas_matching_sources$mesh_id, 2)) |>
  as.data.frame()
names(pheno_pairs_all_data) = c("phenotype_1", "phenotype_2")

## add genetic similarity metrics
# LDSC-based metric
pheno_pairs_all_data = pheno_pairs_all_data %>%
  left_join(ldsc_scores[, c("id1_mesh_id", "id2_mesh_id", "rg")], by = c("phenotype_1" = "id1_mesh_id", "phenotype_2" = "id2_mesh_id")) %>%
  left_join(ldsc_scores[, c("id1_mesh_id", "id2_mesh_id", "rg")], by = c("phenotype_1" = "id2_mesh_id", "phenotype_2" = "id1_mesh_id")) %>%
  group_by(phenotype_1, phenotype_2) %>%
  mutate(ldsc_score = max(rg.x, rg.y, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(ldsc_score = if_else(ldsc_score == -Inf, NA, ldsc_score)) %>%
  dplyr::select(-rg.x, -rg.y)
# MAGMA-based metric
pheno_pairs_all_data = pheno_pairs_all_data %>%
  left_join(magma_pvalues_spearman_cor_all, by = c("phenotype_1" = "Disease1_mesh_id", "phenotype_2" = "Disease2_mesh_id")) %>%
  left_join(magma_pvalues_spearman_cor_all, by = c("phenotype_1" = "Disease2_mesh_id", "phenotype_2" = "Disease1_mesh_id")) %>%
  group_by(phenotype_1, phenotype_2) %>%
  mutate(magma_cor = max(magma_spearman_correlation.x, magma_spearman_correlation.y, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(magma_cor = if_else(magma_cor == -Inf, NA, magma_cor)) %>%
  dplyr::select(-magma_spearman_correlation.x, -magma_spearman_correlation.y)
# S-MultiXcan-based metric
pheno_pairs_all_data = pheno_pairs_all_data %>%
  left_join(smultixcan_spearman_cor_all, by = c("phenotype_1" = "Disease1_mesh_id", "phenotype_2" = "Disease2_mesh_id")) %>%
  left_join(smultixcan_spearman_cor_all, by = c("phenotype_1" = "Disease2_mesh_id", "phenotype_2" = "Disease1_mesh_id")) %>%
  group_by(phenotype_1, phenotype_2) %>%
  mutate(smultixcan_cor = max(smultixcan_spearman_cor.x, smultixcan_spearman_cor.y, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(smultixcan_cor = if_else(smultixcan_cor == -Inf, NA, smultixcan_cor)) %>%
  dplyr::select(-smultixcan_spearman_cor.x, -smultixcan_spearman_cor.y)
# L2G-based metric
pheno_pairs_all_data = pheno_pairs_all_data %>%
  left_join(ot_l2g_spearman_cor_all, by = c("phenotype_1" = "Disease1_mesh_id", "phenotype_2" = "Disease2_mesh_id")) %>%
  left_join(ot_l2g_spearman_cor_all, by = c("phenotype_1" = "Disease2_mesh_id", "phenotype_2" = "Disease1_mesh_id")) %>%
  group_by(phenotype_1, phenotype_2) %>%
  mutate(l2g_cor = max(ot_l2g_spearman_cor.x, ot_l2g_spearman_cor.y, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(l2g_cor = if_else(l2g_cor == -Inf, NA, l2g_cor)) %>%
  dplyr::select(-ot_l2g_spearman_cor.x, -ot_l2g_spearman_cor.y)
# COLOC-based metric
pheno_pairs_all_data = pheno_pairs_all_data %>%
  left_join(ot_gwas_coloc_overlap_coeff[, c("study1_mesh_id", "study2_mesh_id", "ot_coloc_overlap_coefficient")], by = c("phenotype_1" = "study1_mesh_id", "phenotype_2" = "study2_mesh_id")) %>%
  left_join(ot_gwas_coloc_overlap_coeff[, c("study1_mesh_id", "study2_mesh_id", "ot_coloc_overlap_coefficient")], by = c("phenotype_1" = "study2_mesh_id", "phenotype_2" = "study1_mesh_id")) %>%
  group_by(phenotype_1, phenotype_2) %>%
  mutate(coloc_overcoef = max(ot_coloc_overlap_coefficient.x, ot_coloc_overlap_coefficient.y, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(coloc_overcoef = if_else(coloc_overcoef == -Inf, NA, coloc_overcoef)) %>%
  dplyr::select(-ot_coloc_overlap_coefficient.x, -ot_coloc_overlap_coefficient.y)

## keep pairs with genetic similarity metric from at least one resource
pheno_pairs_all_data = pheno_pairs_all_data %>% 
  filter(!if_all(3:7, is.na))

## add phenotype category (body system)
pheno_pairs_all_data = pheno_pairs_all_data %>%
  left_join(gwas_matching_sources[, c("mesh_id", "trait_category")], by = c("phenotype_1" = "mesh_id")) %>%
  left_join(gwas_matching_sources[, c("mesh_id", "trait_category")], by = c("phenotype_2" = "mesh_id")) %>%
  dplyr::relocate(phenotype_1_category = trait_category.x, .before = ldsc_score) %>%
  dplyr::relocate(phenotype_2_category = trait_category.y, .before = ldsc_score) %>%
  mutate(phenotype_category = if_else(phenotype_1_category == phenotype_2_category, "Same", "Different")) %>%
  dplyr::relocate(phenotype_category, .after = phenotype_2_category)
# annotate whether the phenotypes in the pair affect the same or different body system
pheno_pairs_all_data$phenotype_category = factor(
  pheno_pairs_all_data$phenotype_category, 
  levels = c("Same","Different"),
  labels = c("Same","Different")
)

## add drug overlap metrics
# indications overlap
pheno_pairs_all_data = pheno_pairs_all_data %>%
  left_join(drug_ind_overlap, by = c("phenotype_1" = "mesh_id1", "phenotype_2" = "mesh_id2")) %>%
  left_join(drug_ind_overlap, by = c("phenotype_1" = "mesh_id2", "phenotype_2" = "mesh_id1")) %>%
  group_by(phenotype_1, phenotype_2) %>%
  mutate(drug_ind_overlap_coef = max(overlap.x, overlap.y, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-overlap.x, -overlap.y) %>%
  dplyr::relocate(drug_ind_overlap_coef, .after = phenotype_category) %>%
  mutate(drug_ind_overlap_coef = if_else(drug_ind_overlap_coef == -Inf, NA, drug_ind_overlap_coef))
# side effects overlap
pheno_pairs_all_data = pheno_pairs_all_data %>%
  left_join(drug_se_overlap, by = c("phenotype_1" = "mesh_id1", "phenotype_2" = "mesh_id2")) %>%
  left_join(drug_se_overlap, by = c("phenotype_1" = "mesh_id2", "phenotype_2" = "mesh_id1")) %>%
  group_by(phenotype_1, phenotype_2) %>%
  mutate(drug_se_overlap_coef = max(overlap.x, overlap.y, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-overlap.x, -overlap.y) %>%
  dplyr::relocate(drug_se_overlap_coef, .after = phenotype_category) %>%
  mutate(drug_se_overlap_coef = if_else(drug_se_overlap_coef == -Inf, NA, drug_se_overlap_coef))

## save file
fwrite(
  pheno_pairs_all_data, 
  "data/AllPhenoPairs_DrugOverlap_GeneticSimilarity.txt",
  sep = "\t", 
  row.names = FALSE
)
