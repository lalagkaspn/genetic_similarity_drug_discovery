
#### Purpose: create a dataframe with all drug-known indications and tested disease pairs with the corresponding genetic similarity metrics

library(dplyr)
library(data.table)
library(tidyr)
source("src/00_functions.R")

## ======================== ##
## Load and process data    ##
## ======================== ##

## GWAS ID - MeSH ID mapper (links GWAS IDs to MeSH disease IDs)
gwas_matching_sources = fread("data/gwas_matching_between_sources_manual_final_2.txt", na.strings = "")
gwas_matching_sources = gwas_matching_sources %>% 
  dplyr::select(-mr_coloc_pqtl_2020, -mr_coloc_pqtl_2024, -lava)

## drug indications
drugs_indications = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/drugs_indications/drug_indications_chembl_rxnorm_sider_combined_extended_chemblid_mesh.txt") %>%
  filter(max_phase_for_ind == 4) %>%
  distinct()
chembl_drug_names = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/chemlb_id_drugnames.csv") %>%
  dplyr::select(chembl_id = "ChEMBL ID", name = Name, type = Type) %>%
  distinct()
drugs_indications = drugs_indications %>%
  left_join(chembl_drug_names, by = "chembl_id") %>%
  dplyr::select(chembl_id, name, type, mesh_heading, mesh_id, max_phase_for_ind, source, mesh_index) %>%
  filter(mesh_index == "parent") # keep only parent-level MeSH terms
rm(chembl_drug_names)

### Genetic similarity metrics
## LDSC-based metric
ldsc_scores = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/gwasatlas_watanabe/raw/genetic_correlation_v20191115.txt") %>%
  filter(id1 %in% gwas_matching_sources$ldsc_magma & id2 %in% gwas_matching_sources$ldsc_magma) %>% # keep only GWAS present in our sample
  dplyr::select(id1, id2, rg) %>% 
  distinct() %>%
  left_join(gwas_matching_sources[,c("ldsc_magma", "mesh_id")], by = c("id1" = "ldsc_magma")) %>%
  dplyr::rename("id1_mesh_id" = "mesh_id") %>%
  left_join(gwas_matching_sources[,c("ldsc_magma", "mesh_id")], by = c("id2" = "ldsc_magma")) %>%
  dplyr::rename("id2_mesh_id" = "mesh_id")

## MAGMA-based metric
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

## S-MultiXcan-based metric
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

## L2G-based metric
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

## COLOC-based metric
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

## ================================================================================================================= ##
## Create a dataframe with drug - known indications - novel phenotype triplets and  their genetic similarity metrics ##
## ================================================================================================================= ##

## drug indications in our sample
drugs_treating_gwas = drugs_indications %>% 
  filter(mesh_id %in% gwas_matching_sources$mesh_id) %>% 
  dplyr::select(chembl_id, name, type) %>%
  distinct()
# keep molecules that are drugs and keep only the general drug form
molecules_to_remove = c("ALUMINIUM OXIDE", "ALUMINUM CHLORIDE", "ALUMINUM HYDROXIDE", "AMLODIPINE BESYLATE", "AMMONIA SOLUTION, STRONG", "AMMONIUM LACTATE", 
                        "ARIPIPRAZOLE LAUROXIL", "ATOMOXETINE HYDROCHLORIDE", "AZILSARTAN MEDOXOMIL", "BECLOMETHASONE 17-MONOPROPIONATE", "BECLOMETHASONE DIPROPIONATE", 
                        "BEPOTASTINE BESYLATE", "BETAMETHASONE ACETATE", "BETAMETHASONE BENZOATE", "BETAMETHASONE DIPROPIONATE", "BETAMETHASONE PHOSPHORIC ACID", 
                        "BETAMETHASONE VALERATE", "BOTULINUM TOXIN TYPE A", "BROMCRESOL GREEN", "CALCIUM ACETATE", "CALCIUM CARBONATE", "CALCIUM CITRATE", "CALCIUM GLUBIONATE", 
                        "CALCIUM GLUCEPTATE", "CALCIUM GLUCONATE", "CALCIUM LACTATE", "CALCIUM PHOSPHATE, DIBASIC", "CANDESARTAN CILEXETIL", "CEFUROXIME AXETIL", 
                        "CHLORAMPHENICOL PALMITATE", "CHLORAMPHENICOL SUCCINIC ACID", "CINACALCET HYDROCHLORIDE", "CLINDAMYCIN PALMITATE", "CLINDAMYCIN PHOSPHATE", 
                        "DEXAMETHASONE ACETATE", "DEXAMETHASONE PHOSPHORIC ACID", "DIETHYLSTILBESTROL DIPHOSPHATE", "DULOXETINE HYDROCHLORIDE", "ECHOTHIOPHATE IODIDE", 
                        "ENALAPRILAT ANHYDROUS", "ERYTHROMYCIN ETHYLSUCCINATE", "ESLICARBAZEPINE ACETATE", "ESTRADIOL CYPIONATE", "ESTROGENS, CONJUGATED", 
                        "ESTROGENS, ESTERIFIED", "ETOPOSIDE PHOSPHATE", "FERRIC CARBOXYMALTOSE", "FERRIC CITRATE ANHYDROUS", "FERRIC DERISOMALTOSE", "FERRIC MALTOL", 
                        "FERRIC PYROPHOSPHATE CITRATE", "FERROUS FUMARATE", "FERROUS GLUCONATE ANHYDROUS", "FERROUS GLYCINE SULFATE", "FERROUS SUCCINATE", "FERROUS SULFATE", 
                        "FLUOCINOLONE ACETONIDE", "FLUOROMETHOLONE ACETATE", "FLUPHENAZINE DECANOATE", "FLUTICASONE FUROATE", "FOSPHENYTOIN SODIUM", "GABAPENTIN ENACARBIL", 
                        "HALOBETASOL PROPIONATE", "HYDROCORTISONE ACETATE", "HYDROCORTISONE CYPIONATE", "HYDROCORTISONE HEMISUCCINATE ANHYDROUS", "HYDROCORTISONE PHOSPHORIC ACID",
                        "HYDROCORTISONE PROBUTATE", "HYDROCORTISONE VALERATE", "ISOSORBIDE DINITRATE", "ISOSORBIDE MONONITRATE", "KETOROLAC TROMETHAMINE", "LATANOPROSTENE BUNOD", 
                        "LINACLOTIDE ACETATE", "LISDEXAMFETAMINE DIMESYLATE", "LITHIUM CARBONATE", "LITHIUM CITRATE ANHYDROUS", "LITHIUM ION", "LOTEPREDNOL ETABONATE", 
                        "MAGNESIUM CARBONATE", "MAGNESIUM CITRATE", "MAGNESIUM HYDROXIDE", "MAGNESIUM OXIDE", "MAGNESIUM SULFATE ANHYDROUS", "MELOXICAM SODIUM", 
                        "METHYL AMINOLEVULINATE HYDROCHLORIDE", "METHYLDOPA (RACEMIC)", "METHYLPREDNISOLONE ACETATE", "METHYLPREDNISOLONE HEMISUCCINATE", "MOEXIPRIL HYDROCHLORIDE",
                        "MYCOPHENOLATE MOFETIL", "NORETHINDRONE ACETATE", "OLMESARTAN MEDOXOMIL", "PALIPERIDONE PALMITATE", "PLERIXAFOR OCTAHYDROCHLORIDE", "POTASSIUM BITARTRATE",
                        "POTASSIUM IODIDE", "PREDNISOLONE ACETATE", "PREDNISOLONE PHOSPHORIC ACID", "PREDNISOLONE TEBUTATE", "PROPYLENE GLYCOL DIACETATE", 
                        "PROPYLENE GLYCOL MONOSTEARATE", "QUINIDINE POLYGALACTURONATE", "RADIUM DICHLORIDE", "RIMONABANT HYDROCHLORIDE", "RIZATRIPTAN BENZOATE", 
                        "SALMETEROL XINAFOATE", "SODIUM ASCORBATE", "SODIUM BICARBONATE", "SODIUM CARBONATE", "SODIUM CHLORIDE", "SODIUM CITRATE", "SODIUM FEREDETATE", 
                        "SODIUM FLUORIDE", "SODIUM IODIDE", "SODIUM PHOSPHATE", "SODIUM PHOSPHATE, MONOBASIC", "SODIUM SULFATE ANHYDROUS", "TECHNETIUM SESTAMIBI", 
                        "TECHNETIUM TC 99M SULFUR COLLOID", "TECHNETIUM TC 99M TILMANOCEPT", "TRAMETINIB DIMETHYL SULFOXIDE", "TRASTUZUMAB DERUXTECAN", "TRASTUZUMAB EMTANSINE",
                        "TRIAMCINOLONE ACETONIDE", "TRIAMCINOLONE DIACETATE", "TRIAMCINOLONE HEXACETONIDE", "VERATRUM VIRIDE ROOT", "VERNAKALANT HYDROCHLORIDE", 
                        "ZINC ACETATE ANHYDROUS", "ZINC CHLORIDE", "ZINC OXIDE", "ZINC SULFATE", "ZIPRASIDONE HYDROCHLORIDE")
drugs_treating_gwas = drugs_treating_gwas %>% 
  filter(name != "") %>% 
  filter(!name %in% drugs_treating_gwas[which(duplicated(drugs_treating_gwas$name)), name]) %>%
  filter(!name %in% molecules_to_remove) %>%
  filter(type != "Unknown")
drugs_treating_gwas_unique = unique(drugs_treating_gwas$chembl_id)
rm(molecules_to_remove)

## dataframe with all drug-pheno pairs
drug_gwas_pairs = data.frame(
  drug = rep(drugs_treating_gwas_unique, each = nrow(gwas_matching_sources)),
  phenotype = unique(gwas_matching_sources$mesh_id)
)

## add drug indications
drug_gwas_pairs = drug_gwas_pairs %>%
  left_join(drugs_indications[, c("chembl_id", "mesh_id")], by = c("drug" = "chembl_id")) %>%
  dplyr::relocate(drug_indication = mesh_id, .after = drug) %>%
  filter(drug_indication %in% gwas_matching_sources$mesh_id) %>%
  left_join(gwas_matching_sources[,c(1,3)], by = c("drug_indication" = "mesh_id")) %>%
  dplyr::relocate(drug_indication_category = trait_category, .after = drug_indication) %>%
  left_join(gwas_matching_sources[,c(1,3)], by = c("phenotype" = "mesh_id")) %>%
  dplyr::relocate(phenotype_category = trait_category, .after = phenotype)

## add column to indicate if it is a known indication (1) or not (0)
drugs_indications$drug_indicated_for_pheno = 1
drug_gwas_pairs = left_join(drug_gwas_pairs, unique(drugs_indications[, c("chembl_id", "mesh_id", "drug_indicated_for_pheno")]), by = c("drug" = "chembl_id", "phenotype" = "mesh_id"))
drug_gwas_pairs$drug_indicated_for_pheno = ifelse(is.na(drug_gwas_pairs$drug_indicated_for_pheno), 0, 1)

## add genetic similarity metrics
# LDSC-based metric
drug_gwas_pairs = drug_gwas_pairs %>% 
  left_join(ldsc_scores[, c("id1_mesh_id", "id2_mesh_id", "rg")], by = c("drug_indication" = "id1_mesh_id", "phenotype" = "id2_mesh_id")) %>%
  left_join(ldsc_scores[, c("id1_mesh_id", "id2_mesh_id", "rg")], by = c("drug_indication" = "id2_mesh_id", "phenotype" = "id1_mesh_id")) %>%
  group_by(drug_indication, phenotype) %>%
  mutate(ldsc_score = max(rg.x, rg.y, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(ldsc_score = if_else(ldsc_score == -Inf, NA, ldsc_score)) %>%
  dplyr::select(-rg.x, -rg.y)
# MAGMA-based metric
drug_gwas_pairs = left_join(drug_gwas_pairs, magma_pvalues_spearman_cor_all, by = c("drug_indication" = "Disease1_mesh_id", "phenotype" = "Disease2_mesh_id"))
# S-MultiXcan-based metric
drug_gwas_pairs = left_join(drug_gwas_pairs, smultixcan_spearman_cor_all, by = c("drug_indication" = "Disease1_mesh_id", "phenotype" = "Disease2_mesh_id"))
# L2G-based metric
drug_gwas_pairs = left_join(drug_gwas_pairs, ot_l2g_spearman_cor_all, by = c("drug_indication" = "Disease1_mesh_id", "phenotype" = "Disease2_mesh_id"))
# COLOC-based metric
drug_gwas_pairs = drug_gwas_pairs %>%
  left_join(ot_gwas_coloc_overlap_coeff[, c("study1_mesh_id", "study2_mesh_id", "ot_coloc_overlap_coefficient")], by = c("drug_indication" = "study1_mesh_id", "phenotype" = "study2_mesh_id")) %>%
  left_join(ot_gwas_coloc_overlap_coeff[, c("study1_mesh_id", "study2_mesh_id", "ot_coloc_overlap_coefficient")], by = c("drug_indication" = "study2_mesh_id", "phenotype" = "study1_mesh_id")) %>%
  group_by(drug_indication, phenotype) %>%
  mutate(ot_coloc_overlap_coefficient = max(ot_coloc_overlap_coefficient.x, ot_coloc_overlap_coefficient.y, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(ot_coloc_overlap_coefficient = if_else(ot_coloc_overlap_coefficient == -Inf, NA, ot_coloc_overlap_coefficient)) %>%
  dplyr::select(-ot_coloc_overlap_coefficient.x, -ot_coloc_overlap_coefficient.y)

## remove pairs with no genetic similarity metrics
rows_to_remove = which(apply(drug_gwas_pairs[, c(7:11)], MARGIN = 1, FUN = function(x) {sum(is.na(x))}) == 5)
drug_gwas_pairs = drug_gwas_pairs[-rows_to_remove, ]
rownames(drug_gwas_pairs) = NULL
rm(rows_to_remove)

## add disease annotations and phenotypic similarity data
# ClinGraph embedding cosine similarity
clingraph_cosine_df = fread("data/ClinGraph_cosine_similarity.txt")
drug_gwas_pairs = left_join(drug_gwas_pairs, clingraph_cosine_df, by = c("drug_indication" = "id1_mesh", "phenotype" = "id2_mesh"))
# Open Targets same/different body systems
drug_gwas_pairs$same_body_system_OT = ifelse(drug_gwas_pairs$drug_indication_category == drug_gwas_pairs$phenotype_category, "Same", "Different")
# UK HRCS and ICD10
trait_categories = fread("data/phenotype_categories.txt") %>%
  dplyr::select(mesh_id, trait_area_HRCS, trait_ICD10_top)
drug_gwas_pairs = drug_gwas_pairs %>%
  left_join(trait_categories, by = c("drug_indication" = "mesh_id")) %>%
  left_join(trait_categories, by = c("phenotype" = "mesh_id")) %>%
  dplyr::select(drug, drug_indication, drug_indication_category, drug_indication_category_HRCS = trait_area_HRCS.x, drug_indication_category_ICD10 = trait_ICD10_top.x,
                phenotype, phenotype_category, phenotype_category_HRCS = trait_area_HRCS.y, phenotype_category_ICD10 = trait_ICD10_top.y,
                drug_indicated_for_pheno, ldsc_score, magma_spearman_correlation, smultixcan_spearman_cor, ot_l2g_spearman_cor, ot_coloc_overlap_coefficient,
                cosine_similarity_mean, same_body_system_OT) %>%
  mutate(
    same_body_system_HRCS = if_else(drug_indication_category_HRCS == phenotype_category_HRCS, "Same", "Different"),
    same_body_system_ICD10 = if_else(drug_indication_category_ICD10 == phenotype_category_ICD10, "Same", "Different")
  )

## consensus
# if at least one source says that an indication-phenotype pair is phenotypically similar, then we group it as similar ("Same body system")
drug_gwas_pairs$same_body_system_consensus = ifelse(
  drug_gwas_pairs$same_body_system_OT == "Same" | drug_gwas_pairs$same_body_system_HRCS == "Same" | drug_gwas_pairs$same_body_system_ICD10 == "Same", 
  "Same", 
  "Different")
table(drug_gwas_pairs$same_body_system_consensus)
drug_gwas_pairs = drug_gwas_pairs %>% dplyr::select(
  drug, drug_indication, drug_indication_category, phenotype, phenotype_category, drug_indicated_for_pheno:ot_coloc_overlap_coefficient, same_body_system_consensus, cosine_similarity_mean
)

# remove rows where the drug indication is the same as the tested disease
drug_gwas_pairs = drug_gwas_pairs %>%
  filter(drug_indication != phenotype)

## save data frames
fwrite(drug_gwas_pairs, "data/DrugIND_TestedDisease_GenSim.txt.gz", sep = "\t", row.names = FALSE)
