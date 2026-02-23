
#### GOAL: train rstan indications model

.libPaths(c("/home/panagiotisnikolaos_lalagkas_student_uml_edu/R/x86_64-pc-linux-gnu-library/4.2",
            "/usr/local/lib/R/site-library",
            "/usr/local/lib/R/library",
            "/modules/spack/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.4.0/r-4.2.0-3kitpfbxevyhxd2adiznenkjqqdbekzs/rlib/R/library"))
library(data.table)
library(dplyr)
library(parallel)
source('src/00_functions.R')
library(rstan)
options(mc.cores=as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
rstan_options(auto_write = TRUE)

## load package test
print("All packages loaded successfully!")

## hyperparameters
sd_man = 2.5
om = 0.7

print(sd_man)
print(om)

## load data
drug_indication_both = fread("data/DrugIND_TestedDisease_GenSim.txt.gz")
drug_indication_both$same_body_system_ClinGraph = ifelse(drug_indication_both$cosine_similarity_mean >= -0.1, "Same", "Different")
# consensus
# if at least one source says that an indication-phenotype pair is phenotypically similar, then we group it as similar ("Same body system")
drug_indication_both$same_body_system_consensus = ifelse(
  drug_indication_both$same_body_system_ClinGraph == "Same" | drug_indication_both$same_body_system_consensus == "Same", 
  "Same", 
  "Different")
table(drug_indication_both$same_body_system_consensus)
drug_indication_both = drug_indication_both %>% dplyr::select(
  drug, drug_indication, drug_indication_category, phenotype, phenotype_category, drug_indicated_for_pheno:ot_coloc_overlap_coefficient, same_body_system_consensus
)

diff = drug_indication_both %>% 
  filter(same_body_system_consensus == "Different") %>% 
  arrange(drug, phenotype) %>%
  dplyr::select(-same_body_system_consensus) %>%
  as.data.frame()

same = drug_indication_both %>% 
  filter(same_body_system_consensus == "Same") %>% 
  arrange(drug, phenotype) %>%
  dplyr::select(-same_body_system_consensus) %>%
  as.data.frame()

## run stan
pair_names = paste(same$drug, same$phenotype, sep=":")
di_pairs = unique(pair_names)
nsplit = 4
results = list()

splits = rep(1:nsplit, length(di_pairs)/nsplit)
for(i in 1:nsplit){
  cat("RUNNING", i, " with ", sum(splits==i),"\n")
  deb = mult_split_eval_IND(same, 0, c(sd_man), c(om), di_pairs[splits==i])
  results[[i]]= deb[[4]]
}

bound = bind_rows(results, .id = "column_label")
roc(bound$label, bound$p)
save_path = paste0("results/IND_STAN_PhenoSimilar_SD=", sd_man, "_OMULT=", om, "_TestPred.csv", collapse = "")
write.csv(bound,file=save_path)
