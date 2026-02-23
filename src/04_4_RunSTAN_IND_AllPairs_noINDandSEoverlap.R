
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
# indications
drug_indication_both = fread("data/DrugIND_TestedDisease_GenSim.txt.gz")
# side effects
drug_se_both = fread("data/DrugSE_TestedDisease_GenSim.txt.gz")

## remove drug-phenotype pairs that are both indications and side effects
indications = rbind(
  drug_indication_both %>% dplyr::select(drug, drug_indication),
  drug_indication_both %>% filter(drug_indicated_for_pheno == 1) %>% dplyr::select(drug, drug_indication = phenotype)
)
indications = unique(indications)
side_effects = rbind(
  drug_se_both %>% dplyr::select(drug, drug_side_effect),
  drug_se_both %>% filter(is_side_effect == 1) %>% dplyr::select(drug, drug_side_effect = phenotype)
)
side_effects = unique(side_effects)
indications$drug_pheno = paste(indications$drug, indications$drug_indication, sep = "_")
side_effects$drug_pheno = paste(side_effects$drug, side_effects$drug_side_effect, sep = "_")
both_ind_se = intersect(indications$drug_pheno, side_effects$drug_pheno)
drug_indication_both$drug_known = paste(drug_indication_both$drug, drug_indication_both$drug_indication, sep = "_")
drug_se_both$drug_known = paste(drug_se_both$drug, drug_se_both$drug_side_effect, sep = "_")
# remove drug-pheno pairs being both ind and se so not to influence the model training
drug_indication_both$drug_pheno = paste(drug_indication_both$drug, drug_indication_both$phenotype, sep = "_")
drug_se_both$drug_pheno = paste(drug_se_both$drug, drug_se_both$phenotype, sep = "_")
# remove drug-pheno pairs being both ind and se so not to influence the estimation of maximum genetic similarity
drug_indication_both = drug_indication_both %>% filter(!drug_pheno %in% both_ind_se, !drug_known %in% both_ind_se) %>% dplyr::select(-drug_pheno)
drug_se_both = drug_se_both %>% filter(!drug_pheno %in% both_ind_se, !drug_known %in% both_ind_se) %>% dplyr::select(-drug_pheno)

## keep only drugs and phenotypes for which we have information if they are indication and/or side effect
# keep drugs for which we have information about indications and side effects 
common_drug = intersect(drug_indication_both$drug, drug_se_both$drug)
drug_indication_both = drug_indication_both %>%
  filter(drug %in% common_drug) %>%
  distinct()
drug_se_both = drug_se_both %>%
  filter(drug %in% common_drug) %>%
  distinct()
# keep phenotypes for which we have information about being indication or side effect of kept drugs
common_pheno = intersect(drug_indication_both$phenotype, drug_se_both$phenotype)
drug_indication_both = drug_indication_both %>%
  filter(
    drug_indication %in% common_pheno,
    phenotype %in% common_pheno
  ) %>%
  distinct()
drug_se_both = drug_se_both %>%
  filter(
    drug_side_effect %in% common_pheno,
    phenotype %in% common_pheno
  ) %>%
  distinct()

drug_indication_both = drug_indication_both %>% 
  arrange(drug, phenotype) %>%
  dplyr::select(drug, drug_indication, drug_indication_category, phenotype, phenotype_category, drug_indicated_for_pheno, ldsc_score, magma_spearman_correlation, smultixcan_spearman_cor, ot_l2g_spearman_cor, ot_coloc_overlap_coefficient) %>%
  as.data.frame()

## run stan
pair_names = paste(drug_indication_both$drug, drug_indication_both$phenotype, sep=":")
di_pairs = unique(pair_names)
nsplit = 6
results = list()

splits = rep(1:nsplit, length(di_pairs)/nsplit)
for(i in 1:nsplit){
  cat("RUNNING", i, " with ", sum(splits==i),"\n")
  deb = mult_split_eval_IND(drug_indication_both, 0, c(sd_man), c(om), di_pairs[splits==i])
  results[[i]]= deb[[4]]
}

bound = bind_rows(results, .id = "column_label")
roc(bound$label, bound$p)
save_path = paste0("results/IND_STAN_AllPairs_noINDandSEoverlap_SD=", sd_man, "_OMULT=", om, "_TestPred.csv", collapse = "")
write.csv(bound,file=save_path)
