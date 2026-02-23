
#### Purpose: to check the relationship between genetic similarity and true drug indications or side effects across different thresholds of phenotypic similarity

library(dplyr)
library(data.table)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(grid)

### ------- INDICATIONS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

# ====================== #
# Load and process data
# ====================== #

drug_IND = fread("data/DrugIND_TestedDisease_GenSim.txt.gz")

## across different threhsolds of clingraph
results_same = data.frame(
  clingraph_threshold = seq(-0.4, 0.9, 0.1),
  ldsc_same_diff_p = NA,
  magma_same_diff_p = NA,
  smultixcan_same_diff_p = NA,
  l2g_same_diff_p = NA,
  coloc_same_diff_p = NA
)
results_different = data.frame(
  clingraph_threshold = seq(-0.4, 0.9, 0.1),
  ldsc_same_diff_p = NA,
  magma_same_diff_p = NA,
  smultixcan_same_diff_p = NA,
  l2g_same_diff_p = NA,
  coloc_same_diff_p = NA
)

for (i in 1:nrow(results_same)) {
  threshold = results_same[i, "clingraph_threshold"]
  drug_IND = drug_IND %>% 
    mutate(clingraph = if_else(cosine_similarity_mean >= threshold, "Same", "Different")) %>%
    mutate(consensus = if_else(clingraph == "Same", "Same", "Different"))
  drug_IND$consensus = factor(drug_IND$consensus, levels = c("Same", "Different"), labels = c("Same", "Different"))
  
  # split data based on ClinGraph similarity and disease annotations
  drug_IND_Same = drug_IND %>% filter(consensus == "Same")
  drug_IND_Different = drug_IND %>% filter(consensus == "Different")
  
  ## get maximum gen_sim features for each drug-phenotype pair
  # same body system (drug indication - novel phenotype)
  drug_IND_Same = setDT(drug_IND_Same)
  drug_IND_Same = drug_IND_Same[, lapply(.SD, max, na.rm = TRUE), 
                                by = .(drug, phenotype, drug_indicated_for_pheno), 
                                .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
  drug_IND_Same$drug_indicated_for_pheno = factor(drug_IND_Same$drug_indicated_for_pheno, levels = c(1,0), labels = c(1,0))
  # different body system (drug indication - novel phenotype)
  drug_IND_Different = setDT(drug_IND_Different)
  drug_IND_Different = drug_IND_Different[, lapply(.SD, max, na.rm = TRUE), 
                                          by = .(drug, phenotype, drug_indicated_for_pheno), 
                                          .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
  drug_IND_Different$drug_indicated_for_pheno = factor(drug_IND_Different$drug_indicated_for_pheno, levels = c(1,0), labels = c(1,0))
  
  ## SAME
  results_same[i, "ldsc_same_diff_p"] = wilcox.test(ldsc_score ~ drug_indicated_for_pheno, data = drug_IND_Same, alternative = "greater")$p.value
  results_same[i, "magma_same_diff_p"] = wilcox.test(magma_spearman_correlation ~ drug_indicated_for_pheno, data = drug_IND_Same, alternative = "greater")$p.value
  results_same[i, "smultixcan_same_diff_p"] = wilcox.test(smultixcan_spearman_cor ~ drug_indicated_for_pheno, data = drug_IND_Same, alternative = "greater")$p.value
  results_same[i, "l2g_same_diff_p"] = wilcox.test(ot_l2g_spearman_cor ~ drug_indicated_for_pheno, data = drug_IND_Same, alternative = "greater")$p.value
  results_same[i, "coloc_same_diff_p"] = wilcox.test(ot_coloc_overlap_coefficient ~ drug_indicated_for_pheno, data = drug_IND_Same, alternative = "greater")$p.value
  
  ## DIFFERENT
  results_different[i, "ldsc_same_diff_p"] = wilcox.test(ldsc_score ~ drug_indicated_for_pheno, data = drug_IND_Different, alternative = "greater")$p.value
  results_different[i, "magma_same_diff_p"] = wilcox.test(magma_spearman_correlation ~ drug_indicated_for_pheno, data = drug_IND_Different, alternative = "greater")$p.value
  results_different[i, "smultixcan_same_diff_p"] = wilcox.test(smultixcan_spearman_cor ~ drug_indicated_for_pheno, data = drug_IND_Different, alternative = "greater")$p.value
  results_different[i, "l2g_same_diff_p"] = wilcox.test(ot_l2g_spearman_cor ~ drug_indicated_for_pheno, data = drug_IND_Different, alternative = "greater")$p.value
  results_different[i, "coloc_same_diff_p"] = wilcox.test(ot_coloc_overlap_coefficient ~ drug_indicated_for_pheno, data = drug_IND_Different, alternative = "greater")$p.value
}

# Combine results
results <- bind_rows(
  results_same %>% mutate(type = "≥ threshold"),
  results_different %>% mutate(type = "< threshold")
)

# Convert clingraph_threshold to factor (preserving order)
results$clingraph_threshold <- factor(results$clingraph_threshold, 
                                      levels = sort(unique(results$clingraph_threshold)))

# Pivot to long format for easier plotting
results_long <- results %>%
  tidyr::pivot_longer(
    cols = c(ldsc_same_diff_p, magma_same_diff_p, smultixcan_same_diff_p, l2g_same_diff_p, coloc_same_diff_p),
    names_to = "metric",
    values_to = "pvalue"
  ) %>%
  mutate(metric = recode(metric,
                         ldsc_same_diff_p = "LDSC",
                         magma_same_diff_p = "MAGMA",
                         smultixcan_same_diff_p = "S-MultiXcan",
                         l2g_same_diff_p = "L2G",
                         coloc_same_diff_p = "COLOC"))

# Define colors for metrics
metric_colors <- c(
  "LDSC" = "blue",
  "MAGMA" = "red",
  "S-MultiXcan" = "green",
  "L2G" = "purple",
  "COLOC" = "orange"
)

results_long$type = factor(results_long$type, levels = c("≥ threshold", "< threshold"), labels = c("≥ threshold", "< threshold"))
results_long$metric = factor(results_long$metric, levels = results_long$metric, labels = results_long$metric)
# Plot
ggplot(results_long, aes(x = clingraph_threshold, y = -log10(pvalue), color = metric, linetype = type, group = interaction(metric, type))) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(values = metric_colors) +
  labs(
    x = "Phenotypic similarity threshold\n(ClinGraph cosine similarity)",
    y = expression(-log[10](p)),
    title = "",
    linetype = "Drug indication-disease\nphenotypic similarity",
    color = ""
  ) +
  theme_classic(base_size = 16, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  guides(color = "none") +
  facet_wrap(~metric, scales = "free", nrow = 2)
ggsave("figures/Supplementary_figure_5.png",
       device = "png", dpi = 600, width = 20, height = 10)

rm(list = ls())

### ------- SIDE EFFECTS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ####

## load data
drug_SE = fread("data/DrugSE_TestedDisease_GenSim.txt.gz")

## across different threhsolds of clingraph
results_same = data.frame(
  clingraph_threshold = seq(-0.4, 0.9, 0.1),
  ldsc_same_diff_p = NA,
  magma_same_diff_p = NA,
  smultixcan_same_diff_p = NA,
  l2g_same_diff_p = NA,
  coloc_same_diff_p = NA
)
results_different = data.frame(
  clingraph_threshold = seq(-0.4, 0.9, 0.1),
  ldsc_same_diff_p = NA,
  magma_same_diff_p = NA,
  smultixcan_same_diff_p = NA,
  l2g_same_diff_p = NA,
  coloc_same_diff_p = NA
)

for (i in 1:nrow(results_same)) {
  threshold = results_same[i, "clingraph_threshold"]
  drug_SE = drug_SE %>% 
    mutate(clingraph = if_else(cosine_similarity_mean >= threshold, "Same", "Different")) %>%
    mutate(consensus = if_else(clingraph == "Same", "Same", "Different"))
  drug_SE$consensus = factor(drug_SE$consensus, levels = c("Same", "Different"), labels = c("Same", "Different"))
  
  # split data based on ClinGraph similarity and disease annotations
  drug_SE_Same = drug_SE %>% filter(consensus == "Same")
  drug_SE_Different = drug_SE %>% filter(consensus == "Different")
  
  ## get maximum gen_sim features for each drug-phenotype pair
  # same body system (drug side effect - novel phenotype)
  drug_SE_Same = setDT(drug_SE_Same)
  drug_SE_Same = drug_SE_Same[, lapply(.SD, max, na.rm = TRUE), 
                              by = .(drug, phenotype, is_side_effect), 
                              .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
  drug_SE_Same$is_side_effect = factor(drug_SE_Same$is_side_effect, levels = c(1,0), labels = c(1,0))
  # different body system (drug side effect - novel phenotype)
  drug_SE_Different = setDT(drug_SE_Different)
  drug_SE_Different = drug_SE_Different[, lapply(.SD, max, na.rm = TRUE), 
                                        by = .(drug, phenotype, is_side_effect), 
                                        .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
  drug_SE_Different$is_side_effect = factor(drug_SE_Different$is_side_effect, levels = c(1,0), labels = c(1,0))
  
  ## SAME
  results_same[i, "ldsc_same_diff_p"] = wilcox.test(ldsc_score ~ is_side_effect, data = drug_SE_Same, alternative = "greater")$p.value
  results_same[i, "magma_same_diff_p"] = wilcox.test(magma_spearman_correlation ~ is_side_effect, data = drug_SE_Same, alternative = "greater")$p.value
  results_same[i, "smultixcan_same_diff_p"] = wilcox.test(smultixcan_spearman_cor ~ is_side_effect, data = drug_SE_Same, alternative = "greater")$p.value
  results_same[i, "l2g_same_diff_p"] = wilcox.test(ot_l2g_spearman_cor ~ is_side_effect, data = drug_SE_Same, alternative = "greater")$p.value
  results_same[i, "coloc_same_diff_p"] = wilcox.test(ot_coloc_overlap_coefficient ~ is_side_effect, data = drug_SE_Same, alternative = "greater")$p.value
  
  ## DIFFERENT
  results_different[i, "ldsc_same_diff_p"] = wilcox.test(ldsc_score ~ is_side_effect, data = drug_SE_Different, alternative = "greater")$p.value
  results_different[i, "magma_same_diff_p"] = wilcox.test(magma_spearman_correlation ~ is_side_effect, data = drug_SE_Different, alternative = "greater")$p.value
  results_different[i, "smultixcan_same_diff_p"] = wilcox.test(smultixcan_spearman_cor ~ is_side_effect, data = drug_SE_Different, alternative = "greater")$p.value
  results_different[i, "l2g_same_diff_p"] = wilcox.test(ot_l2g_spearman_cor ~ is_side_effect, data = drug_SE_Different, alternative = "greater")$p.value
  results_different[i, "coloc_same_diff_p"] = wilcox.test(ot_coloc_overlap_coefficient ~ is_side_effect, data = drug_SE_Different, alternative = "greater")$p.value
}

# Combine results
results <- bind_rows(
  results_same %>% mutate(type = "≥ threshold"),
  results_different %>% mutate(type = "< threshold")
)

# Convert clingraph_threshold to factor (preserving order)
results$clingraph_threshold <- factor(
  results$clingraph_threshold, 
  levels = sort(unique(results$clingraph_threshold)))

# Pivot to long format for easier plotting
results_long <- results %>%
  tidyr::pivot_longer(
    cols = c(ldsc_same_diff_p, magma_same_diff_p, smultixcan_same_diff_p, l2g_same_diff_p, coloc_same_diff_p),
    names_to = "metric",
    values_to = "pvalue"
  ) %>%
  mutate(metric = recode(metric,
                         ldsc_same_diff_p = "LDSC",
                         magma_same_diff_p = "MAGMA",
                         smultixcan_same_diff_p = "S-MultiXcan",
                         l2g_same_diff_p = "L2G",
                         coloc_same_diff_p = "COLOC"))

# Define colors for metrics
metric_colors <- c(
  "LDSC" = "blue",
  "MAGMA" = "red",
  "S-MultiXcan" = "green",
  "L2G" = "purple",
  "COLOC" = "orange"
)

results_long$type = factor(results_long$type, levels = c("≥ threshold", "< threshold"), labels = c("≥ threshold", "< threshold"))
results_long$metric = factor(results_long$metric, levels = results_long$metric, labels = results_long$metric)
# Plot
ggplot(results_long, aes(x = clingraph_threshold, y = -log10(pvalue), color = metric, linetype = type, group = interaction(metric, type))) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(values = metric_colors) +
  labs(
    x = "Phenotypic similarity threshold\n(Clingraph cosine similarity)",
    y = expression(-log[10](p)),
    title = "",
    linetype = "Drug side effect-disease\nphenotypic similarity",
    color = ""
  ) +
  theme_classic(base_size = 16, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  guides(color = "none") +
  facet_wrap(~metric, scales = "free", nrow = 2)
ggsave("figures/Supplementary_figure_7.png",
       device = "png", dpi = 600, width = 20, height = 10)
