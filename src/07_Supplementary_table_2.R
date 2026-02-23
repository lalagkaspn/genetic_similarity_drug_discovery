
#### Create supplementary table 2 - side effects

library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)

##  - Drug (chembl_id and name)
##  - Disease (mesh_id, name)
##  - True label
##  - STAN prediction (using all pairs or phenotypically dissimilar pairs - two separate sheets)
##  - Drug's known side effects (mesh_id, name)

wb = createWorkbook()

diseases_names = fread("data/gwas_matching_between_sources_manual_final_2.txt")
diseases_names = diseases_names[,c("mesh_id","mesh_name")]

drug_names = fread("data/chemlb_id_drugnames.csv") %>%
  dplyr::select(chembl_id = "ChEMBL ID", drug_name = Name) %>%
  distinct()

## Phenotypically dissimilar drug side effect - disease pairs 
preds_SE_PhenoDissim = fread("results/SE_STAN_PhenoDissimilar_SD=3.75_OMULT=0.6_TestPred.csv")
data_SE_all = fread("data/DrugSE_TestedDisease_GenSim.txt.gz") %>%
  dplyr::select(drug, drug_side_effect, disease = phenotype, is_side_effect, same_body_system_consensus, cosine_similarity_mean) %>%
  distinct()
# flag drug side effect - disease pairs that are phenotypically similar
data_SE_all$same_body_system_consensus = ifelse(data_SE_all$same_body_system_consensus == "Same" | data_SE_all$cosine_similarity_mean >= -0.1, "Same", "Different")

preds_SE_PhenoDissim = cbind(
  stringr::str_split_fixed(preds_SE_PhenoDissim$V1, pattern = ":", n = 2),
  preds_SE_PhenoDissim[,c("p","label")]
) 
data_SE_all = data_SE_all %>%
  left_join(preds_SE_PhenoDissim, by = c("drug" = "V1", "disease" = "V2", "is_side_effect" = "label")) %>%
  na.omit() %>% # the drug-disease pairs removed are those that all drug side effects are phenotypically similar to the tested disease
  filter(same_body_system_consensus == "Different") %>%
  dplyr::select(-same_body_system_consensus, -cosine_similarity_mean)
data_SE_all = data_SE_all %>% 
  left_join(drug_names, by = c("drug" = "chembl_id")) %>%
  dplyr::relocate(drug_name, .after = drug) %>%
  left_join(diseases_names, by = c("drug_side_effect" = "mesh_id")) %>%
  dplyr::relocate(drug_side_effect_name = mesh_name, .after = p) %>%
  left_join(diseases_names, by = c("disease" = "mesh_id")) %>%
  dplyr::relocate(disease_name = mesh_name, .after = disease)
data_SE_all = data_SE_all[, .(
  drug_name = paste(unique(drug_name), collapse = " ||| "),
  drug_side_effect = paste(unique(drug_side_effect), collapse = " ||| "),
  drug_side_effect_name = paste(unique(drug_side_effect_name), collapse = " ||| "),
  disease_name = first(disease_name),
  is_side_effect = first(is_side_effect),
  p = first(p)),
  by = .(drug, disease)
]
data_SE_all = data_SE_all %>%
  dplyr::select(drug, disease, drug_name, disease_name, is_side_effect, predicted_prob = p, drug_side_effect, drug_side_effect_name) %>%
  distinct()

# add excel sheet
addWorksheet(wb, "pheno_dissimilar_pairs")
writeData(wb, "pheno_dissimilar_pairs", data_SE_all)


## All drug side effect - disease pairs 
preds_SE_all = fread("results/SE_STAN_AllPairs_SD=3.75_OMULT=0.6_TestPred.csv")
data_SE_all = fread("data/DrugSE_TestedDisease_GenSim.txt.gz") %>%
  dplyr::select(drug, drug_side_effect, disease = phenotype, is_side_effect) %>%
  distinct()

preds_SE_all = cbind(
  stringr::str_split_fixed(preds_SE_all$V1, pattern = ":", n = 2),
  preds_SE_all[,c("p","label")]
) 
data_SE_all = data_SE_all %>%
  left_join(preds_SE_all, by = c("drug" = "V1", "disease" = "V2", "is_side_effect" = "label"))
data_SE_all = data_SE_all %>% 
  left_join(drug_names, by = c("drug" = "chembl_id")) %>%
  dplyr::relocate(drug_name, .after = drug) %>%
  left_join(diseases_names, by = c("drug_side_effect" = "mesh_id")) %>%
  dplyr::relocate(drug_side_effect_name = mesh_name, .after = p) %>%
  left_join(diseases_names, by = c("disease" = "mesh_id")) %>%
  dplyr::relocate(disease_name = mesh_name, .after = disease)
data_SE_all = data_SE_all[, .(
  drug_name = paste(unique(drug_name), collapse = " ||| "),
  drug_side_effect = paste(unique(drug_side_effect), collapse = " ||| "),
  drug_side_effect_name = paste(unique(drug_side_effect_name), collapse = " ||| "),
  disease_name = first(disease_name),
  is_side_effect = first(is_side_effect),
  p = first(p)),
  by = .(drug, disease)
]
data_SE_all = data_SE_all %>%
  dplyr::select(drug, disease, drug_name, disease_name, is_side_effect, predicted_prob = p, drug_side_effect, drug_side_effect_name) %>%
  distinct()

# add excel sheet
addWorksheet(wb, "all_pairs")
writeData(wb, "all_pairs", data_SE_all)

# save excel file
saveWorkbook(wb, "tables/Supplemenary_table_2.xlsx", overwrite = TRUE)
