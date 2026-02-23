
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
sd_man = 3.75
om = 0.6

print(sd_man)
print(om)

## load data
drug_se_both = fread("data/DrugSE_TestedDisease_GenSim.txt.gz")
drug_se_both =  drug_se_both %>%
  dplyr::rename(drug_indication = drug_side_effect, drug_indication_category = drug_side_effect_category, drug_indicated_for_pheno = is_side_effect) %>%
  arrange(drug, phenotype) %>%
  dplyr::select(-same_body_system_consensus, -cosine_similarity_mean) %>%
  as.data.frame()

## run stan
pair_names = paste(drug_se_both$drug, drug_se_both$phenotype, sep=":")
di_pairs = unique(pair_names)
nsplit = 4
results = list()

splits = rep(1:nsplit, length(di_pairs)/nsplit)
for(i in 1:nsplit){
  cat("RUNNING", i, " with ", sum(splits==i),"\n")
  deb = mult_split_eval_SE(drug_se_both, 0, c(3), c(.6), di_pairs[splits==i])
  results[[i]]= deb[[4]]
}

bound = bind_rows(results, .id = "column_label")
roc(bound$label, bound$p)
save_path = paste0("results/SE_STAN_AllPairs_SD=", sd_man, "_OMULT=", om, "_TestPred.csv", collapse = "")
write.csv(bound,file=save_path)
