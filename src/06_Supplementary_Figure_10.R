
#### Purpose: calculate the percentage of predicted drug-disease pairs where the drug targets genes with direct genetic evidence (GWAS variants below a p-value threhsold)

library(data.table)
library(dplyr)
library(biomaRt)
library(GenomicRanges)
library(ggplot2)
library(viridis)
source("src/00_functions.R")

## In this script, we use the predictions from the indications model that used phenotypically dissimilar drug indication-disease pairs.
## We obtained GWAS summary statistics from the GWAS Catalog for each disease whenever available (99 of 178 diseases).
## For each phenotype, we performed the following steps:
## - identify all unique drugs with predicted probabilities for that phenotype.
## - retrieve the gene targets of each drug and their genomic coordinates.
## - extract GWAS p-values for SNPs located within the genomic regions of those gene targets.
## - for each drug, find the minimum GWAS p-value across all SNPs mapped to its gene targets.
## - create binary indicator variables denoting whether the minimum p-value falls below predefined significance thresholds:
##     - 5e-08 (genome-wide significance)
##     - 5e-06
##     - 5e-04
##     - 5e-02
## To run this script, you must first download the disease-specific GWAS summary statistics from the GWAS Catalog.
## As redistribution of these data is not permitted, we instead provide the processed results required to reproduce the supplementary figure 10.
## The code below loads these results and regenerates Supplementary Figure 10.
## For transparency and reproducibility, we also include the full analysis code so users can review the complete workflow.

## load data with direct genetic support information
gensup = fread("results/DirectGeneticEvidence_SupFig10.txt", sep="\t")
## baseline probability of a drug-disease pair being indication - we get it from the test set data
IND_preds = fread("results/IND_STAN_PhenoDissimilar_SD=2.5_OMULT=0.7_TestPred.csv")
baseline_prob = prop.table(table(IND_preds$label))[2]
## split the data into drug-disease pairs with predicted probability below or above the baseline
gensup_above_baseline = gensup %>% filter(pred_prob > baseline_prob)
gensup_below_baseline = gensup %>% filter(pred_prob <= baseline_prob)

## calculate percentage of drug-disease pairs with direct genetic evidence across different threhsolds of GWAS p-value and predicted probability
# pairs with probability above the baseline
x = gensup_above_baseline
x = reshape2::melt(x, c("drug","phenotype","pred_prob"), colnames(x)[5:ncol(x)])
x$variable = gsub("pval_less_5Eneg08", 5e-08, x$variable)
x$variable = gsub("pval_less_5Eneg06", 5e-06, x$variable)
x$variable = gsub("pval_less_5Eneg04", 5e-04, x$variable)
x$variable = gsub("pval_less_5Eneg02", 5e-02, x$variable)
x$variable = gsub("pval_less_1Eneg01", 1e-01, x$variable)
x$variable = as.numeric(x$variable)
colnames(x) = c("drug", "phenotype", "pred_prob", "gwas_threshold", "has_SNP")
gwas_support_above_baseline = expand.grid(
  prob_threshold = c(0.01, 0.1, 0.2),
  gwas_threshold = c(5e-08, 5e-06, 5e-04, 5e-02),
  stringsAsFactors = FALSE
)
gwas_support_above_baseline$predicted = NA
gwas_support_above_baseline$supported_count = NA
gwas_support_above_baseline$supported_percentage = NA
for (i in 1:nrow(gwas_support_above_baseline)) {
  prob_thresh = gwas_support_above_baseline[i,"prob_threshold"]
  gwas_thresh = gwas_support_above_baseline[i,"gwas_threshold"]
  temp = x %>%
    filter(pred_prob >= prob_thresh) %>%
    filter(gwas_threshold == gwas_thresh)
  gwas_support_above_baseline[i, "predicted"] = nrow(temp)
  gwas_support_above_baseline[i, "supported_count"] = sum(temp$has_SNP == 1)
  gwas_support_above_baseline[i, "supported_percentage"] = sum(temp$has_SNP == 1) / nrow(temp) * 100
}
# pairs with probability below the baseline
x = gensup_below_baseline
x = reshape2::melt(x, c("drug","phenotype","pred_prob"), colnames(x)[5:ncol(x)])
x$variable = gsub("pval_less_5Eneg08", 5e-08, x$variable)
x$variable = gsub("pval_less_5Eneg06", 5e-06, x$variable)
x$variable = gsub("pval_less_5Eneg04", 5e-04, x$variable)
x$variable = gsub("pval_less_5Eneg02", 5e-02, x$variable)
x$variable = gsub("pval_less_1Eneg01", 1e-01, x$variable)
x$variable = as.numeric(x$variable)
colnames(x) = c("drug", "phenotype", "pred_prob", "gwas_threshold", "has_SNP")
gwas_support_below_baseline = expand.grid(
  prob_threshold = c(0),
  gwas_threshold = c(5e-08, 5e-06, 5e-04, 5e-02),
  stringsAsFactors = FALSE
)
gwas_support_below_baseline$predicted = NA
gwas_support_below_baseline$supported_count = NA
gwas_support_below_baseline$supported_percentage = NA
for (i in 1:nrow(gwas_support_below_baseline)) {
  gwas_thresh = gwas_support_below_baseline[i,"gwas_threshold"]
  temp = x %>%
    filter(pred_prob <= baseline_prob) %>%
    filter(gwas_threshold == gwas_thresh)
  gwas_support_below_baseline[i, "predicted"] = nrow(temp)
  gwas_support_below_baseline[i, "supported_count"] = sum(temp$has_SNP == 1)
  gwas_support_below_baseline[i, "supported_percentage"] = sum(temp$has_SNP == 1) / nrow(temp) * 100
}

# ggplot(gwas_support_above_baseline, aes(x = factor(gwas_threshold), y = factor(prob_threshold),
#                                           fill = supported_percentage)) +
#   geom_tile(color = "white", linewidth = 0.5) +
#   geom_text(aes(label = sprintf("%.01f%%", supported_percentage)),
#             color = "white", size = 5, fontface = "bold") +
#   scale_fill_viridis_c(
#     option = "plasma",
#     name = "% Supported",
#     labels = function(x) paste0(x, "%"),
#     guide = guide_colorbar(barwidth = 1, barheight = 15),
#     begin = 0.15,  # Skip the darkest purple/black
#     direction = -1,
#     end = 0.85     # Stop before the lightest yellow
#   ) +
#   # scale_x_continuous(breaks = unique(x$gwas_threshold), expand = c(0, 0)) +
#   scale_y_discrete(
#     name = "Prediction Probability Threshold",
#     labels = function(x) sprintf("%.3f", as.numeric(x))
#   ) +
#   labs(
#     title = "GWAS Support Across Prediction Thresholds",
#     subtitle = "Darker colors indicate higher percentage of GWAS support",
#     x = "GWAS Threshold"
#   ) +
#   theme_minimal(base_size = 20, base_family = "Arial") +
#   theme(
#     plot.title = element_text(face = "bold", size = 24, hjust = 0.5, color = "black"),
#     plot.subtitle = element_text(size = 18, hjust = 0.5, color = "black"),
#     axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
#     panel.grid = element_blank(),
#     legend.position = "right"
#   )

ultimate = expand.grid(
  gwas_threshold = c(5e-08, 5e-06, 5e-04, 5e-02),
  top_pred_threshold = c(0.01, 0.1, 0.2),
  pvalue_onesided = NA,
  nr_pred_pairs = NA
)
for (i in 1:nrow(ultimate)) {
  gwas_threshold = ultimate[i, "gwas_threshold"]
  pred_threshold = ultimate[i, "top_pred_threshold"]
  x = gensup_above_baseline %>% filter(pred_prob >= pred_threshold)
  
  if (gwas_threshold == 5e-08) {
    temp = prop.test(
      c(sum(x$pval_less_5Eneg08), sum(gensup_below_baseline$pval_less_5Eneg08)), 
      c(nrow(x), nrow(gensup_below_baseline)),
      alternative = "two.sided", correct = FALSE)
    ultimate[i, "pvalue_onesided"] = temp$p.value
    ultimate[i, "nr_pred_pairs"] = nrow(x)
  }
  
  if (gwas_threshold == 5e-06) {
    temp = prop.test(
      c(sum(x$pval_less_5Eneg06), sum(gensup_below_baseline$pval_less_5Eneg06)), 
      c(nrow(x), nrow(gensup_below_baseline)),
      alternative = "two.sided", correct = FALSE)
    ultimate[i, "pvalue_onesided"] = temp$p.value
    ultimate[i, "nr_pred_pairs"] = nrow(x)
  }
  
  if (gwas_threshold == 5e-04) {
    temp = prop.test(
      c(sum(x$pval_less_5Eneg04), sum(gensup_below_baseline$pval_less_5Eneg04)), 
      c(nrow(x), nrow(gensup_below_baseline)),
      alternative = "two.sided", correct = FALSE)
    ultimate[i, "pvalue_onesided"] = temp$p.value
    ultimate[i, "nr_pred_pairs"] = nrow(x)
  }
  
  if (gwas_threshold == 5e-02) {
    temp = prop.test(
      c(sum(x$pval_less_5Eneg02), sum(gensup_below_baseline$pval_less_5Eneg02)), 
      c(nrow(x), nrow(gensup_below_baseline)),
      alternative = "two.sided", correct = FALSE)
    ultimate[i, "pvalue_onesided"] = temp$p.value
    ultimate[i, "nr_pred_pairs"] = nrow(x)
  }
}

ultimate_plot = rbind(
  gwas_support_below_baseline %>% 
    dplyr::select(prob_threshold, gwas_threshold, supported_count, supported_percentage, predicted) %>%
    mutate(type = "below_baseline") %>%
    filter(prob_threshold == 0), 
  gwas_support_above_baseline %>%
    dplyr::select(prob_threshold, gwas_threshold, supported_count, supported_percentage, predicted) %>% 
    mutate(type = "above_baseline") %>%
    filter(prob_threshold >= 0.01, prob_threshold <= 0.2)
)
ultimate$pvalue_stars = ifelse(ultimate$pvalue_onesided < 0.05, "*", "")
ultimate$pvalue_stars = ifelse(ultimate$pvalue_onesided < 0.01, "**", ultimate$pvalue_stars)
ultimate$pvalue_stars = ifelse(ultimate$pvalue_onesided < 0.001, "***", ultimate$pvalue_stars)
ultimate_plot = ultimate_plot %>%
  left_join(ultimate[,c("gwas_threshold","top_pred_threshold","pvalue_stars")], by = c("gwas_threshold", "prob_threshold" = "top_pred_threshold")) %>%
  mutate(pvalue_stars = if_else(is.na(pvalue_stars), "", pvalue_stars)) %>%
  group_by(gwas_threshold) %>%
  mutate(y_val = max(supported_percentage) + 3)
ultimate_plot$prob_threshold = gsub("0.01", "≥0.01", ultimate_plot$prob_threshold)
ultimate_plot$prob_threshold = gsub("0.1", "≥0.1", ultimate_plot$prob_threshold)
ultimate_plot$prob_threshold = gsub("0.2", "≥0.2", ultimate_plot$prob_threshold)
ultimate_plot$prob_threshold[1:4] = "≤0.004 (baseline prob)"
ultimate_plot$gwas_threshold = paste("<", ultimate_plot$gwas_threshold, sep = "")
ultimate_plot$gwas_threshold = factor(ultimate_plot$gwas_threshold, levels = ultimate_plot$gwas_threshold, labels = ultimate_plot$gwas_threshold)
ggplot(ultimate_plot, aes(x = factor(gwas_threshold), y = supported_percentage, fill = factor(prob_threshold))) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.2) +
  scale_fill_manual(values = c(
    "#B0B0B0",
    "#9ECAE1",
    "#4292C6",
    "#08306B"
  )) +
  geom_text(position = position_dodge(width = 0.8),
    aes(
      x = factor(gwas_threshold),
      y = y_val,
      label = pvalue_stars
    ),
    size = 5, fontface = "bold") +
  labs(
    x = "GWAS threshold",
    y = "% of genetically supported drug-disease pairs",
    fill = "Drug-disease pairs\npredicted probability"
  ) +
  scale_y_continuous(breaks = seq(0,100, 10)) +
  theme_bw(base_size = 20, base_family = "Arial") +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
ggsave("figures/Supplementary_figure_10.png",
       device = "png", dpi = 600, width = 12, height = 7)


#### code we used to generate the "DirectGeneticEvidence_SupFig10.txt" file ####

# ## predictions from the indication model (phenotypically dissimilar pairs)
# IND_preds = fread("results/IND_STAN_PhenoDissimilar_SD=2.5_OMULT=0.7_TestPred.csv")
# # split drug-disease pairs into two columns
# IND_preds = cbind(
#   stringr::str_split_fixed(IND_preds$V1, pattern = ":", n = 2),
#   IND_preds[, c(3,4)]
# )
# colnames(IND_preds) = c("drug", "phenotype", "PredProb", "is_indicated")
# # baseline probability of being an indication
# baseline_prob = prop.table(table(IND_preds$is_indicated))[2]
# 
# ## drug gene-targets
# # from Open Targets Platform v22.09
# drug_targets = fread("data/drug_targets.txt")
# # filter for drugs in our sample
# drug_targets = drug_targets %>% filter(drug %in% IND_preds$drug)
# 
# ## for each gene-target, find position in the genome (both GRCh37 and GRCh38)
# ensg = unique(drug_targets$target) # all unique gene-targets
# # hg38
# mart_hg38 = useEnsembl(
#   biomart = "genes",
#   dataset = "hsapiens_gene_ensembl",
#   GRCh = 38 # it is the default, ignore the warning
# )
# hg38_coords = getBM(
#   attributes = c("ensembl_gene_id", "chromosome_name",
#                  "start_position", "end_position"),
#   filters = "ensembl_gene_id",
#   values = ensg,
#   mart = mart_hg38
# )
# colnames(hg38_coords) = c("target", "chromosome_hg38", "start_hg38", "end_hg38")
# # expand gene regions by +-5kb to include variants upstream and downstrea
# hg38_coords$start_hg38 = hg38_coords$start_hg38 - 5000
# hg38_coords$end_hg38 = hg38_coords$end_hg38 + 5000
# # hg19
# mart_hg19 = useEnsembl(
#   biomart = "genes",
#   dataset = "hsapiens_gene_ensembl",
#   GRCh = 37
# )
# hg19_coords = getBM(
#   attributes = c("ensembl_gene_id", "chromosome_name",
#                  "start_position", "end_position"),
#   filters = "ensembl_gene_id",
#   values = ensg,
#   mart = mart_hg19
# )
# colnames(hg19_coords) = c("target", "chromosome_hg19", "start_hg19", "end_hg19")
# # expand gene regions by +-5kb to include variants upstream and downstrea
# hg19_coords$start_hg19 = hg19_coords$start_hg19 - 5000
# hg19_coords$end_hg19 = hg19_coords$end_hg19 + 5000
# # add genomic positions to drug gene-target dataframe
# drug_targets = drug_targets %>%
#   left_join(hg38_coords, by = "target") %>%
#   left_join(hg19_coords, by = "target") %>%
#   na.omit()
# # remove gene-targets in MT/X/Y chromosomes
# drug_targets = drug_targets %>%
#   filter(!chromosome_hg38 %in% c("MT", "X", "Y")) %>%
#   filter(!chromosome_hg19 %in% c("MT", "X", "Y"))
# drug_targets$chromosome_hg38 = paste("chr", drug_targets$chromosome_hg38, sep = "")
# drug_targets$chromosome_hg19 = paste("chr", drug_targets$chromosome_hg19, sep = "")
# 
# # filter predictions for drugs with known gene-targets
# IND_preds = IND_preds %>% filter(drug %in% drug_targets$drug)
# 
# ## empty dataframe to store the calculated percentage of predicted drug-disease pairs with genetic support for disease GWAS variants below a certain p-value threhsold
# gensup = data.frame()
# 
# ## we obtained GWAS summary statistics from the GWAS Catalog for each disease whenever available (99 of 178 diseases).
# ## for each phenotype, we performed the following steps:
# ## - identify all unique drugs with predicted probabilities for that phenotype.
# ## - retrieve the gene targets of each drug and their genomic coordinates.
# ## - extract GWAS p-values for SNPs located within the genomic regions of those gene targets.
# ## - for each drug, find the minimum GWAS p-value across all SNPs mapped to its gene targets.
# ## - create binary indicator variables denoting whether the minimum p-value falls below predefined significance thresholds:
# ##     - 5e-08 (genome-wide significance)
# ##     - 5e-06
# ##     - 5e-04
# ##     - 5e-02
# ##     - 0.1
# 
# #### D000152 ####
# temp_pheno = "D000152"
# temp_gwas = fread("D000152_hg19.tsv")
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D000544 ####
# temp_pheno = "D000544"
# temp_gwas = fread("D000544_hg19.txt.gz")
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location = PosGRCh37, p_value = p) %>%
#   distinct()
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D000690 ####
# temp_pheno = "D000690"
# temp_gwas = fread("D000690_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D000740 ####
# temp_pheno = "D000740"
# temp_gwas = fread("D000740_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# ### D000787 ####
# temp_pheno = "D000787"
# temp_gwas = fread("D000787_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D001008 ####
# temp_pheno = "D001008"
# temp_gwas = fread("D001008_hg19")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P.value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D001145 ####
# temp_pheno = "D001145"
# temp_gwas = fread("D001145_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D001172 ####
# temp_pheno = "D001172"
# temp_gwas = fread("D001172_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# #### D001201 ####
# temp_pheno = "D001201"
# temp_gwas = fread("D001201_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D001249 ####
# temp_pheno = "D001249"
# temp_gwas = fread("D001249_hg19.txt.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome= CHR, base_pair_location = BP, p_value = P) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# #### D001281 ####
# temp_pheno = "D001281"
# temp_gwas = fread("D001281_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# #### D001289 ####
# temp_pheno = "D001289"
# temp_gwas = fread("D001289_hg19")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# #### D001714 ####
# temp_pheno = "D001714"
# temp_gwas = fread("D001714_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = "#CHROM", base_pair_location = POS, p_value = PVAL) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# #### D001851 ####
# temp_pheno = "D001851"
# temp_gwas = fread("D001851_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# 
# #### D001943 ####
# temp_pheno = "D001943"
# temp_gwas = fread("D001943_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# #### D001991 ####
# temp_pheno = "D001991"
# temp_gwas = fread("D001991_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D002294 ####
# temp_pheno = "D002294"
# temp_gwas = fread("D002294_hg19.txt")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = 'Chr ID', base_pair_location = 'Chr Position', p_value = 'P-value') %>%
#   distinct() %>%
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# #### D002349 ####
# temp_pheno = "D002349"
# temp_gwas = fread("D002349_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>%
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D002386 ####
# temp_pheno = "D002386"
# temp_gwas = fread("D002386_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>%
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D002481 ####
# temp_pheno = "D002481"
# temp_gwas = fread("D002481_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>%
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# 
# #### D002543 ####
# temp_pheno = "D002543"
# temp_gwas = fread("D002543_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D002544 ####
# temp_pheno = "D002544"
# temp_gwas = fread("D002544_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D002764 ####
# temp_pheno = "D002764"
# temp_gwas = fread("D002764_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D002769 ####
# temp_pheno = "D002769"
# temp_gwas = fread("D002769_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D003093 ####
# temp_pheno = "D003093"
# temp_gwas = fread("D003093_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D003248 ####
# temp_pheno = "D003248"
# temp_gwas = fread("D003248_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D003424 ####
# temp_pheno = "D003424"
# temp_gwas = fread("D003424_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D003556 ####
# temp_pheno = "D003556"
# temp_gwas = fread("D003556_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D003680 ####
# temp_pheno = "D003680"
# temp_gwas = fread("D003680_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# #### D003922 ####
# temp_pheno = "D003922"
# temp_gwas = fread("D003922_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# 
# #### D004415 ####
# temp_pheno = "D004415"
# temp_gwas = fread("D004415_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D004485 ####
# temp_pheno = "D004485"
# temp_gwas = fread("D004485_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D004941 ####
# temp_pheno = "D004941"
# temp_gwas = fread("D004941_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D005076 ####
# temp_pheno = "D005076"
# temp_gwas = fread("D005076_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D005242 ####
# temp_pheno = "D005242"
# temp_gwas = fread("D005242_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D005348 ####
# temp_pheno = "D005348"
# temp_gwas = fread("D005348_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D005705 ####
# temp_pheno = "D005705"
# temp_gwas = fread("D005705_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D005764 ####
# temp_pheno = "D005764"
# temp_gwas = fread("D005764_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D005901 ####
# temp_pheno = "D005901"
# temp_gwas = fread("D005901_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D006042 ####
# temp_pheno = "D006042"
# temp_gwas = fread("D006042_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D006073 ####
# temp_pheno = "D006073"
# temp_gwas = fread("D006073_hg19.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D006333 ####
# temp_pheno = "D006333"
# temp_gwas = fread("D006333_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D006552 ####
# temp_pheno = "D006552"
# temp_gwas = fread("D006552_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D006937 ####
# temp_pheno = "D006937"
# temp_gwas = fread("D006937_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D006973 ####
# temp_pheno = "D006973"
# temp_gwas = fread("D006973_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D007022 ####
# temp_pheno = "D007022"
# temp_gwas = fread("D007022_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D007037 ####
# temp_pheno = "D007037"
# temp_gwas = fread("D007037_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D007676 ####
# temp_pheno = "D007676"
# temp_gwas = fread("D007676_hg19.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D007762 ####
# temp_pheno = "D007762"
# temp_gwas = fread("D007762_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D007766 ####
# temp_pheno = "D007766"
# temp_gwas = fread("D007766_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D008105 ####
# temp_pheno = "D008105"
# temp_gwas = fread("D008105_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D008175 ####
# temp_pheno = "D008175"
# temp_gwas = fread("D008175_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D008228 ####
# temp_pheno = "D008228"
# temp_gwas = fread("D008228_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D008286 ####
# temp_pheno = "D008286"
# temp_gwas = fread("D008286_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D009203 ####
# temp_pheno = "D009203"
# temp_gwas = fread("D009203_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D009260 ####
# temp_pheno = "D009260"
# temp_gwas = fread("D009260_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D009771 ####
# temp_pheno = "D009771"
# temp_gwas = fread("D009771_hg19.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D010024 ####
# temp_pheno = "D010024"
# temp_gwas = fread("D010024_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D010195 ####
# temp_pheno = "D010195"
# temp_gwas = fread("D010195_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D010689 ####
# temp_pheno = "D010689"
# temp_gwas = fread("D010689_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D010996 ####
# temp_pheno = "D010996"
# temp_gwas = fread("D010996_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D010998 ####
# temp_pheno = "D010998"
# temp_gwas = fread("D010998_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D011030 ####
# temp_pheno = "D011030"
# temp_gwas = fread("D011030_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D011111 ####
# temp_pheno = "D011111"
# temp_gwas = fread("D011111_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D011141 ####
# temp_pheno = "D011141"
# temp_gwas = fread("D011141_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D011471 ####
# temp_pheno = "D011471"
# temp_gwas = fread("D011471_hg19.txt")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = `Chr ID`, base_pair_location = `Chr Position`, p_value = `P-value`) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D011472 ####
# temp_pheno = "D011472"
# temp_gwas = fread("D011472_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D011537 ####
# temp_pheno = "D011537"
# temp_gwas = fread("D011537_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D011565 ####
# temp_pheno = "D011565"
# temp_gwas = fread("D011565_hg19.txt")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = `Chr ID`, base_pair_location = `Chr Position`, p_value = `P-value`) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D011704 ####
# temp_pheno = "D011704"
# temp_gwas = fread("D011704_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D012163 ####
# temp_pheno = "D012163"
# temp_gwas = fread("D012163_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D012559 ####
# temp_pheno = "D012559"
# temp_gwas = fread("D012559_hg19.txt.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D012852 ####
# temp_pheno = "D012852"
# temp_gwas = fread("D012852_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D012878 ####
# temp_pheno = "D012878"
# temp_gwas = fread("D012878_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D012891 ####
# temp_pheno = "D012891"
# temp_gwas = fread("D012891_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D013167 ####
# temp_pheno = "D013167"
# temp_gwas = fread("D013167_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D013276 ####
# temp_pheno = "D013276"
# temp_gwas = fread("D013276_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D013345 ####
# temp_pheno = "D013345"
# temp_gwas = fread("D013345_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D013575 ####
# temp_pheno = "D013575"
# temp_gwas = fread("D013575_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D014012 ####
# temp_pheno = "D014012"
# temp_gwas = fread("D014012_hg19.txt.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D014060 ####
# temp_pheno = "D014060"
# temp_gwas = fread("D014060_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D014581 ####
# temp_pheno = "D014581"
# temp_gwas = fread("D014581_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D014648 ####
# temp_pheno = "D014648"
# temp_gwas = fread("D014648_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D015179 ####
# temp_pheno = "D015179"
# temp_gwas = fread("D015179_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D015212 ####
# temp_pheno = "D015212"
# temp_gwas = fread("D015212_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D015228 ####
# temp_pheno = "D015228"
# temp_gwas = fread("D015228_hg19.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D016055 ####
# temp_pheno = "D016055"
# temp_gwas = fread("D016055_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D016491 ####
# temp_pheno = "D016491"
# temp_gwas = fread("D016491_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D017098 ####
# temp_pheno = "D017098"
# temp_gwas = fread("D017098_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D018798 ####
# temp_pheno = "D018798"
# temp_gwas = fread("D018798_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D019226 ####
# temp_pheno = "D019226"
# temp_gwas = fread("D019226_hg19.txt.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D019966 ####
# temp_pheno = "D019966"
# temp_gwas = fread("D019966_hg38.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 38
#   )
# )
# 
# #### D020177 ####
# temp_pheno = "D020177"
# temp_gwas = fread("D020177_hg19.txt.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D043183 ####
# temp_pheno = "D043183"
# temp_gwas = fread("D043183_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D051436 ####
# temp_pheno = "D051436"
# temp_gwas = fread("D051436_hg19.txt.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = Chr, base_pair_location = Pos_b37, p_value = `P-value`) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D058186 ####
# temp_pheno = "D058186"
# temp_gwas = fread("D058186_hg19.tsv.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome, base_pair_location, p_value) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D065631 ####
# temp_pheno = "D065631"
# temp_gwas = fread("D065631_hg19.txt.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D000067877 ####
# temp_pheno = "D000067877"
# temp_gwas = fread("D000067877_hg19.txt.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D000076385 ####
# temp_pheno = "D000076385"
# temp_gwas = fread("D000076385_hg19.txt.gz")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = CHR, base_pair_location = BP, p_value = P) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
# 
# #### D000083242 ####
# temp_pheno = "D000083242"
# temp_gwas = fread("D000083242_hg19.txt")
# temp_gwas = temp_gwas %>%
#   dplyr::select(chromosome = `Chr ID`, base_pair_location = `Chr Position`, p_value = `P-value`) %>%
#   distinct() %>% 
#   na.omit()
# temp_gwas$chromosome = paste("chr", temp_gwas$chromosome, sep = "")
# gensup = rbind(
#   gensup, 
#   gensup_subthreshold(
#     preds = IND_preds,
#     drug_targets = drug_targets, 
#     gwas = temp_gwas, 
#     pheno = temp_pheno,
#     genome_v = 19
#   )
# )
#
# # save
# fwrite(gensup, "results/DirectGeneticEvidence_SupFig10.txt", sep="\t", row.names = FALSE)
