
#### Create supplementary table 1 - indications

library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)

##  - Drug (chembl_id and name)
##  - Disease (mesh_id, name)
##  - True label
##  - Max clinical trial phase reached
##  - STAN prediction (using all pairs or phenotypically dissimilar pairs - two separate sheets)
##  - Drug's known indications (mesh_id, name)

wb = createWorkbook()

diseases_names = fread("data/gwas_matching_between_sources_manual_final_2.txt")
diseases_names = diseases_names[,c("mesh_id","mesh_name")]

drug_names = fread("data/chemlb_id_drugnames.csv") %>%
  dplyr::select(chembl_id = "ChEMBL ID", drug_name = Name) %>%
  distinct()

drug_ct = fread("data/DrugDisease_PerCT.txt")

## Phenotypically dissimilar drug indication - disease pairs 
preds_IND_PhenoDissim = fread("results/IND_STAN_PhenoDissimilar_SD=2.5_OMULT=0.7_TestPred.csv")
data_IND_all = fread("data/DrugIND_TestedDisease_GenSim.txt.gz") %>%
  dplyr::select(drug, drug_indication, disease = phenotype, is_indicated = drug_indicated_for_pheno, same_body_system_consensus, cosine_similarity_mean) %>%
  distinct()
# flag drug indication - disease pairs that are phenotypically similar
data_IND_all$same_body_system_consensus = ifelse(data_IND_all$same_body_system_consensus == "Same" | data_IND_all$cosine_similarity_mean >= -0.1, "Same", "Different")

preds_IND_PhenoDissim = cbind(
  stringr::str_split_fixed(preds_IND_PhenoDissim$V1, pattern = ":", n = 2),
  preds_IND_PhenoDissim[,c("p","label")]
) 
data_IND_all = data_IND_all %>%
  left_join(preds_IND_PhenoDissim, by = c("drug" = "V1", "disease" = "V2", "is_indicated" = "label")) %>%
  na.omit() %>% # the drug-disease pairs removed are those that all drug indications are phenotypically similar to the tested disease
  filter(same_body_system_consensus == "Different") %>%
  dplyr::select(-same_body_system_consensus, -cosine_similarity_mean)
data_IND_all = data_IND_all %>% 
  left_join(drug_names, by = c("drug" = "chembl_id")) %>%
  dplyr::relocate(drug_name, .after = drug) %>%
  left_join(diseases_names, by = c("drug_indication" = "mesh_id")) %>%
  dplyr::relocate(drug_indication_name = mesh_name, .after = p) %>%
  left_join(diseases_names, by = c("disease" = "mesh_id")) %>%
  dplyr::relocate(disease_name = mesh_name, .after = disease)
data_IND_all = data_IND_all[, .(
  drug_name = paste(unique(drug_name), collapse = " ||| "),
  drug_indication = paste(unique(drug_indication), collapse = " ||| "),
  drug_indication_name = paste(unique(drug_indication_name), collapse = " ||| "),
  disease_name = first(disease_name),
  is_indicated = first(is_indicated),
  p = first(p)),
  by = .(drug, disease)
]
data_IND_all = data_IND_all %>%
  dplyr::select(drug, disease, drug_name, disease_name, is_indicated, predicted_prob = p, drug_indication, drug_indication_name) %>%
  distinct() %>%
  left_join(drug_ct, by = c("drug", "disease" = "phenotype")) %>%
  dplyr::relocate(max_CT_phase = max_phase_for_ind, .after = is_indicated)
# for drugs and phenotypes that do not exist in the clinical trials file, we don't know if they have ever been in clinical trials --> set "status unknown"
shared_drugs = intersect(drug_ct$drug, data_IND_all$drug)
shared_phenos = intersect(drug_ct$phenotype, data_IND_all$disease)
data_IND_all$max_CT_phase = ifelse(!data_IND_all$drug %in% shared_drugs | !data_IND_all$disease %in% shared_phenos, "Status unknown", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(is.na(data_IND_all$max_CT_phase), "Never tested", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "Not Applicable", "", data_IND_all$max_CT_phase) # clinical trials of devices or behavioral interventions (according to ClinicalTrials.gov)
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "", "Phase unknown", data_IND_all$max_CT_phase) # has been in clinical trials but phase is unknown
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "1", "Phase I", data_IND_all$max_CT_phase) 
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "2", "Phase II", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "3", "Phase III", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "4", "Approved", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(data_IND_all$is_indicated == 1, "Approved", data_IND_all$max_CT_phase) # some of drug-disease pairs with max_clin_phase_unknown are already approved --> fix that
table(data_IND_all$is_indicated, data_IND_all$max_CT_phase)

# add excel sheet
addWorksheet(wb, "pheno_dissimilar_pairs")
writeData(wb, "pheno_dissimilar_pairs", data_IND_all)

## All drug indication - disease pairs 
preds_IND_all = fread("results/IND_STAN_AllPairs_SD=2.5_OMULT=0.7_TestPred.csv")
data_IND_all = fread("data/DrugIND_TestedDisease_GenSim.txt.gz") %>%
  dplyr::select(drug, drug_indication, disease = phenotype, is_indicated = drug_indicated_for_pheno) %>%
  distinct()

preds_IND_all = cbind(
  stringr::str_split_fixed(preds_IND_all$V1, pattern = ":", n = 2),
  preds_IND_all[,c("p","label")]
) 
data_IND_all = data_IND_all %>%
  left_join(preds_IND_all, by = c("drug" = "V1", "disease" = "V2", "is_indicated" = "label"))
data_IND_all = data_IND_all %>% 
  left_join(drug_names, by = c("drug" = "chembl_id")) %>%
  dplyr::relocate(drug_name, .after = drug) %>%
  left_join(diseases_names, by = c("drug_indication" = "mesh_id")) %>%
  dplyr::relocate(drug_indication_name = mesh_name, .after = p) %>%
  left_join(diseases_names, by = c("disease" = "mesh_id")) %>%
  dplyr::relocate(disease_name = mesh_name, .after = disease)
data_IND_all = data_IND_all[, .(
  drug_name = paste(unique(drug_name), collapse = " ||| "),
  drug_indication = paste(unique(drug_indication), collapse = " ||| "),
  drug_indication_name = paste(unique(drug_indication_name), collapse = " ||| "),
  disease_name = first(disease_name),
  is_indicated = first(is_indicated),
  p = first(p)),
  by = .(drug, disease)
]
data_IND_all = data_IND_all %>%
  dplyr::select(drug, disease, drug_name, disease_name, is_indicated, predicted_prob = p, drug_indication, drug_indication_name) %>%
  distinct() %>%
  left_join(drug_ct, by = c("drug", "disease" = "phenotype")) %>%
  dplyr::relocate(max_CT_phase = max_phase_for_ind, .after = is_indicated)
# for drugs and phenotypes that do not exist in the clinical trials file, we don't know if they have ever been in clinical trials --> set "status unknown"
shared_drugs = intersect(drug_ct$drug, data_IND_all$drug)
shared_phenos = intersect(drug_ct$phenotype, data_IND_all$disease)
data_IND_all$max_CT_phase = ifelse(!data_IND_all$drug %in% shared_drugs | !data_IND_all$disease %in% shared_phenos, "Status unknown", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(is.na(data_IND_all$max_CT_phase), "Never tested", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "Not Applicable", "", data_IND_all$max_CT_phase) # clinical trials of devices or behavioral interventions (according to ClinicalTrials.gov)
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "", "Phase unknown", data_IND_all$max_CT_phase) # has been in clinical trials but phase is unknown
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "1", "Phase I", data_IND_all$max_CT_phase) 
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "2", "Phase II", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "3", "Phase III", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(data_IND_all$max_CT_phase == "4", "Approved", data_IND_all$max_CT_phase)
data_IND_all$max_CT_phase = ifelse(data_IND_all$is_indicated == 1, "Approved", data_IND_all$max_CT_phase) # some of drug-disease pairs with max_clin_phase_unknown are already approved --> fix that
table(data_IND_all$is_indicated, data_IND_all$max_CT_phase)

# add excel sheet
addWorksheet(wb, "all_pairs")
writeData(wb, "all_pairs", data_IND_all)

# save excel file
saveWorkbook(wb, "tables/Supplemenary_table_1.xlsx", overwrite = TRUE)
