
#### Purpose: to check the relationship between genetic similarity and true drug indications and side effects stratified by body system (OT, UK HRCS, ICD10)

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

### add disease annotations (UK HRCS & ICD10)
trait_categories = fread("data/phenotype_categories.txt") %>%
  dplyr::select(mesh_id, trait_area_HRCS, trait_ICD10_top)
drug_IND = drug_IND %>%
  left_join(trait_categories, by = c("drug_indication" = "mesh_id")) %>%
  left_join(trait_categories, by = c("phenotype" = "mesh_id")) %>%
  dplyr::select(drug, drug_indication, phenotype, drug_indicated_for_pheno, ldsc_score:ot_coloc_overlap_coefficient,
                drug_indication_category_OT = drug_indication_category, phenotype_category_OT = phenotype_category,
                drug_indication_category_HRCS = trait_area_HRCS.x, phenotype_category_HRCS = trait_area_HRCS.y,
                drug_indication_category_ICD10 = trait_ICD10_top.x, phenotype_category_ICD10 = trait_ICD10_top.y) %>%
  mutate(
    category = if_else(
      drug_indication_category_OT == phenotype_category_OT | drug_indication_category_HRCS == phenotype_category_HRCS | drug_indication_category_ICD10 == phenotype_category_ICD10,
      "Same", "Different")) %>%
  dplyr::select(drug, drug_indication, phenotype, category, drug_indicated_for_pheno, ldsc_score:ot_coloc_overlap_coefficient)

# split data based on disease annotations
drug_IND_Same = drug_IND %>% filter(category == "Same")
drug_IND_Different = drug_IND %>% filter(category == "Different")
rm(drug_IND)

## get maximum gen_sim features for each drug-phenotype pair
# same body system (drug indication - novel phenotype)
drug_IND_Same = setDT(drug_IND_Same)
drug_IND_Same = drug_IND_Same[, lapply(.SD, max, na.rm = TRUE), 
                              by = .(drug, phenotype, drug_indicated_for_pheno), 
                              .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
drug_IND_Same$drug_indicated_for_pheno = factor(drug_IND_Same$drug_indicated_for_pheno, levels = c(0,1), labels = c(0,1))
# different body system (drug indication - novel phenotype)
drug_IND_Different = setDT(drug_IND_Different)
drug_IND_Different = drug_IND_Different[, lapply(.SD, max, na.rm = TRUE), 
                                        by = .(drug, phenotype, drug_indicated_for_pheno), 
                                        .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
drug_IND_Different$drug_indicated_for_pheno = factor(drug_IND_Different$drug_indicated_for_pheno, levels = c(0,1), labels = c(0,1))

# ================= #
# Rank-sum tests
# ================= #

## same body system (drug indication - novel phenotype)
# LDSC score
data_ldsc_ind_same = drug_IND_Same %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, ldsc_score) %>% 
  na.omit()
data_ldsc_ind_same$drug_indicated_for_pheno = factor(data_ldsc_ind_same$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_ldsc_ind_same_noOutliers = data_ldsc_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(ldsc_score, 0.25),
    Q3 = quantile(ldsc_score, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(ldsc_score >= (Q1 - 1.5 * IQR) & ldsc_score <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_ldsc_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_ldsc_ind_same = ggplot(data_ldsc_ind_same, aes(x = drug_indicated_for_pheno, y = ldsc_score, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "LDSC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_ldsc_ind_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_ldsc_ind_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_ldsc_ind_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_ldsc_ind_same_noOutliers$ldsc_score) + 0.15) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_ldsc_ind_same_noOutliers$ldsc_score) + 0.1*max(data_ldsc_ind_same_noOutliers$ldsc_score), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_ldsc_ind_same_noOutliers$ldsc_score), max(data_ldsc_ind_same_noOutliers$ldsc_score) + 0.2*max(data_ldsc_ind_same_noOutliers$ldsc_score))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# MAGMA spearman correlation
data_magma_ind_same = drug_IND_Same %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, magma_cor = magma_spearman_correlation) %>% 
  na.omit()
data_magma_ind_same$drug_indicated_for_pheno = factor(data_magma_ind_same$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_magma_ind_same_noOutliers = data_magma_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(magma_cor, 0.25),
    Q3 = quantile(magma_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(magma_cor >= (Q1 - 1.5 * IQR) & magma_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_magma_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_magma_ind_same = ggplot(data_magma_ind_same, aes(x = drug_indicated_for_pheno, y = magma_cor, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "MAGMA",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_magma_ind_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_magma_ind_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_magma_ind_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_magma_ind_same_noOutliers$magma_cor) + 0.01*max(data_magma_ind_same_noOutliers$magma_cor)) +
  scale_y_continuous(breaks = seq(-0.03, 0.25, 0.05)) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_magma_ind_same_noOutliers$magma_cor) + 0.1*max(data_magma_ind_same_noOutliers$magma_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_magma_ind_same_noOutliers$magma_cor), max(data_magma_ind_same_noOutliers$magma_cor) + 0.2*max(data_magma_ind_same_noOutliers$magma_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# S-MultiXcan spearman correlaton
data_smultixcan_ind_same = drug_IND_Same %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, smultixcan_cor = smultixcan_spearman_cor) %>% 
  na.omit()
data_smultixcan_ind_same$drug_indicated_for_pheno = factor(data_smultixcan_ind_same$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_smultixcan_ind_same_noOutliers = data_smultixcan_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(smultixcan_cor, 0.25),
    Q3 = quantile(smultixcan_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(smultixcan_cor >= (Q1 - 1.5 * IQR) & smultixcan_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_smultixcan_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_smultixcan_ind_same = ggplot(data_smultixcan_ind_same, aes(x = drug_indicated_for_pheno, y = smultixcan_cor, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "S-MultiXcan",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_smultixcan_ind_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_smultixcan_ind_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_smultixcan_ind_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_smultixcan_ind_same_noOutliers$smultixcan_cor) - 0.015) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_smultixcan_ind_same_noOutliers$smultixcan_cor) + 0.1*max(data_smultixcan_ind_same_noOutliers$smultixcan_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_smultixcan_ind_same_noOutliers$smultixcan_cor), max(data_smultixcan_ind_same_noOutliers$smultixcan_cor) + 0.2*max(data_smultixcan_ind_same_noOutliers$smultixcan_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# Open Targets L2G spearman correlation
data_l2g_ind_same = drug_IND_Same %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, l2g_cor = ot_l2g_spearman_cor) %>% 
  na.omit()
data_l2g_ind_same$drug_indicated_for_pheno = factor(data_l2g_ind_same$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_l2g_ind_same_noOutliers = data_l2g_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(l2g_cor, 0.25),
    Q3 = quantile(l2g_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(l2g_cor >= (Q1 - 1.5 * IQR) & l2g_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_l2g_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_l2g_ind_same = ggplot(data_l2g_ind_same, aes(x = drug_indicated_for_pheno, y = l2g_cor, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "Open Targets L2G",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_l2g_ind_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_l2g_ind_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_l2g_ind_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_l2g_ind_same_noOutliers$l2g_cor) + 0.1*max(data_l2g_ind_same_noOutliers$l2g_cor)) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_l2g_ind_same_noOutliers$l2g_cor) + 0.1*max(data_l2g_ind_same_noOutliers$l2g_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_l2g_ind_same_noOutliers$l2g_cor), max(data_l2g_ind_same_noOutliers$l2g_cor) + 0.2*max(data_l2g_ind_same_noOutliers$l2g_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# Open Targets COLOC overlap coefficient
data_coloc_ind_same = drug_IND_Same %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, coloc_overcoef = ot_coloc_overlap_coefficient) %>% 
  na.omit()
data_coloc_ind_same$drug_indicated_for_pheno = factor(data_coloc_ind_same$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_coloc_ind_same_noOutliers = data_coloc_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(coloc_overcoef, 0.25),
    Q3 = quantile(coloc_overcoef, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(coloc_overcoef >= (Q1 - 1.5 * IQR) & coloc_overcoef <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_coloc_ind_same %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_coloc_ind_same = ggplot(data_coloc_ind_same, aes(x = drug_indicated_for_pheno, y = coloc_overcoef, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "COLOC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_coloc_ind_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_coloc_ind_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_coloc_ind_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_coloc_ind_same_noOutliers$coloc_overcoef) + 0.09) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_coloc_ind_same_noOutliers$coloc_overcoef) + 0.1*max(data_coloc_ind_same_noOutliers$coloc_overcoef), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_coloc_ind_same_noOutliers$coloc_overcoef), max(data_coloc_ind_same_noOutliers$coloc_overcoef) + 0.2*max(data_coloc_ind_same_noOutliers$coloc_overcoef))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")

IND_same = annotate_figure(
  ggarrange(plot_ldsc_ind_same, plot_magma_ind_same, plot_smultixcan_ind_same, plot_l2g_ind_same, plot_coloc_ind_same, 
            ncol = 5, nrow = 1, align = "hv"), 
  left = text_grob(
    "Max genetic similarity between\nknown indications of a drug and a test disease",
    size = 20, face = "bold", family = "Arial", rot = 90
  ),
  bottom = text_grob(
    "Is the test disease a known indication of the drug?",
    size = 20, face = "bold", family = "Arial"
  ),       
) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

## different body system (drug indication - novel phenotype)
# LDSC score
data_ldsc_ind_different = drug_IND_Different %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, ldsc_score) %>% 
  na.omit()
data_ldsc_ind_different$drug_indicated_for_pheno = factor(data_ldsc_ind_different$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_ldsc_ind_different_noOutliers = data_ldsc_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(ldsc_score, 0.25),
    Q3 = quantile(ldsc_score, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(ldsc_score >= (Q1 - 1.5 * IQR) & ldsc_score <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_ldsc_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_ldsc_ind_different = ggplot(data_ldsc_ind_different, aes(x = drug_indicated_for_pheno, y = ldsc_score, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "LDSC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_ldsc_ind_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_ldsc_ind_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_ldsc_ind_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_ldsc_ind_different_noOutliers$ldsc_score) + 0.06) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_ldsc_ind_different_noOutliers$ldsc_score) + 0.05, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_ldsc_ind_different_noOutliers$ldsc_score), max(data_ldsc_ind_different_noOutliers$ldsc_score) + 0.15)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# MAGMA spearman correlation
data_magma_ind_different = drug_IND_Different %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, magma_cor = magma_spearman_correlation) %>% 
  na.omit()
data_magma_ind_different$drug_indicated_for_pheno = factor(data_magma_ind_different$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_magma_ind_different_noOutliers = data_magma_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(magma_cor, 0.25),
    Q3 = quantile(magma_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(magma_cor >= (Q1 - 1.5 * IQR) & magma_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_magma_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_magma_ind_different = ggplot(data_magma_ind_different, aes(x = drug_indicated_for_pheno, y = magma_cor, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "MAGMA",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_magma_ind_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_magma_ind_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_magma_ind_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(-0.03, 0.15, 0.05)) +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_magma_ind_different_noOutliers$magma_cor) + 0.031) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_magma_ind_different_noOutliers$magma_cor) + 0.025, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_magma_ind_different_noOutliers$magma_cor), max(data_magma_ind_different_noOutliers$magma_cor) + 0.05)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# S-MultiXcan spearman correlaton
data_smultixcan_ind_different = drug_IND_Different %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, smultixcan_cor = smultixcan_spearman_cor) %>% 
  na.omit()
data_smultixcan_ind_different$drug_indicated_for_pheno = factor(data_smultixcan_ind_different$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_smultixcan_ind_different_noOutliers = data_smultixcan_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(smultixcan_cor, 0.25),
    Q3 = quantile(smultixcan_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(smultixcan_cor >= (Q1 - 1.5 * IQR) & smultixcan_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_smultixcan_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_smultixcan_ind_different = ggplot(data_smultixcan_ind_different, aes(x = drug_indicated_for_pheno, y = smultixcan_cor, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "S-MultiXcan",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_smultixcan_ind_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_smultixcan_ind_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_smultixcan_ind_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(-0.03, 0.07, by = 0.02)) + 
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_smultixcan_ind_different_noOutliers$smultixcan_cor) + 0.005) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_smultixcan_ind_different_noOutliers$smultixcan_cor) + 0.0125, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_smultixcan_ind_different_noOutliers$smultixcan_cor), max(data_smultixcan_ind_different_noOutliers$smultixcan_cor) + 0.025)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# Open Targets L2G spearman correlation
data_l2g_ind_different = drug_IND_Different %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, l2g_cor = ot_l2g_spearman_cor) %>% 
  na.omit()
data_l2g_ind_different$drug_indicated_for_pheno = factor(data_l2g_ind_different$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_l2g_ind_different_noOutliers = data_l2g_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(l2g_cor, 0.25),
    Q3 = quantile(l2g_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(l2g_cor >= (Q1 - 1.5 * IQR) & l2g_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_l2g_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_l2g_ind_different = ggplot(data_l2g_ind_different, aes(x = drug_indicated_for_pheno, y = l2g_cor, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "Open Targets L2G",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_l2g_ind_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_l2g_ind_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_l2g_ind_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_l2g_ind_different_noOutliers$l2g_cor) + 0.008) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_l2g_ind_different_noOutliers$l2g_cor) + 0.018, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_l2g_ind_different_noOutliers$l2g_cor), max(data_l2g_ind_different_noOutliers$l2g_cor) + 0.05)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# Open Targets COLOC overlap coefficient
data_coloc_ind_different = drug_IND_Different %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, coloc_overcoef = ot_coloc_overlap_coefficient) %>% 
  na.omit()
data_coloc_ind_different$drug_indicated_for_pheno = factor(data_coloc_ind_different$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_coloc_ind_different_noOutliers = data_coloc_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(coloc_overcoef, 0.25),
    Q3 = quantile(coloc_overcoef, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(coloc_overcoef >= (Q1 - 1.5 * IQR) & coloc_overcoef <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_coloc_ind_different %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_coloc_ind_different = ggplot(data_coloc_ind_different, aes(x = drug_indicated_for_pheno, y = coloc_overcoef, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "COLOC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_coloc_ind_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_coloc_ind_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_coloc_ind_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_coloc_ind_different_noOutliers$coloc_overcoef)) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_coloc_ind_different_noOutliers$coloc_overcoef) + 0.03, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_coloc_ind_different_noOutliers$coloc_overcoef), max(data_coloc_ind_different_noOutliers$coloc_overcoef) + 0.075)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")

IND_different = annotate_figure(
  ggarrange(plot_ldsc_ind_different, plot_magma_ind_different, plot_smultixcan_ind_different, plot_l2g_ind_different, plot_coloc_ind_different, 
            ncol = 5, nrow = 1, align = "hv"), 
  left = text_grob(
    "Max genetic similarity between\nknown indications of a drug and a test disease",
    size = 20, face = "bold", family = "Arial", rot = 90
  ),
  bottom = text_grob(
    "Is the test disease a known indication of the drug?",
    size = 20, face = "bold", family = "Arial"
  ),   
) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# plot all together
ggarrange(
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  grobTree(rectGrob(gp = gpar(fill = "#BF3984", col = NA)),
           textGrob("Same body system", gp = gpar(col = "white", fontface = "bold", fontsize = 26))),
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  IND_same,
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  grobTree(rectGrob(gp = gpar(fill = "#FCCE25", col = NA)),
           textGrob("Different body system", gp = gpar(col = "black", fontface = "bold", fontsize = 26))),
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  IND_different,
  ncol = 1,
  heights = c(0.3, 0.55, 0.3, 7, 0.3, 0.55, 0.3, 7)
)
ggsave("figures/Supplementary_figure_4.png",
       device = "png", dpi = 600, width = 26, height = 15)

rm(list = ls())

### ------- SIDE EFFECTS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ####

# ====================== #
# Load and process data
# ====================== #

drug_SE = fread("data/DrugSE_TestedDisease_GenSim.txt.gz")

### add disease annotations (UK HRCS & ICD10)
trait_categories = fread("data/phenotype_categories.txt") %>%
  dplyr::select(mesh_id, trait_area_HRCS, trait_ICD10_top)
drug_SE = drug_SE %>%
  left_join(trait_categories, by = c("drug_side_effect" = "mesh_id")) %>%
  left_join(trait_categories, by = c("phenotype" = "mesh_id")) %>%
  dplyr::select(drug, drug_side_effect, phenotype, is_side_effect, ldsc_score:ot_coloc_overlap_coefficient,
                drug_side_effect_category_OT = drug_side_effect_category, phenotype_category_OT = phenotype_category,
                drug_side_effect_category_HRCS = trait_area_HRCS.x, phenotype_category_HRCS = trait_area_HRCS.y,
                drug_side_effect_category_ICD10 = trait_ICD10_top.x, phenotype_category_ICD10 = trait_ICD10_top.y) %>%
  mutate(
    category = if_else(
      drug_side_effect_category_OT == phenotype_category_OT | drug_side_effect_category_HRCS == phenotype_category_HRCS | drug_side_effect_category_ICD10 == phenotype_category_ICD10,
      "Same", "Different")) %>%
  dplyr::select(drug, drug_side_effect, phenotype, category, is_side_effect, ldsc_score:ot_coloc_overlap_coefficient)

# split data based on ClinGraph similarity and disease annotations
drug_SE_Same = drug_SE %>% filter(category == "Same")
drug_SE_Different = drug_SE %>% filter(category == "Different")

## for each drug-novel_phenotype pair, keep  maximum gen_sim among all known drug side effects
# same body system (drug side effect - novel phenotype)
drug_SE_Same = setDT(drug_SE_Same)
drug_SE_Same = drug_SE_Same[, lapply(.SD, max, na.rm = TRUE), 
                            by = .(drug, phenotype, is_side_effect), 
                            .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
drug_SE_Same$is_side_effect = factor(drug_SE_Same$is_side_effect, levels = c(0,1), labels = c(0,1))
# different body system (drug side effect - novel phenotype)
drug_SE_Different = setDT(drug_SE_Different)
drug_SE_Different = drug_SE_Different[, lapply(.SD, max, na.rm = TRUE), 
                                      by = .(drug, phenotype, is_side_effect), 
                                      .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
drug_SE_Different$is_side_effect = factor(drug_SE_Different$is_side_effect, levels = c(0,1), labels = c(0,1))

# ================= #
# Rank-sum tests
# ================= #

## same body system (drug indication - novel phenotype)
# LDSC score
data_ldsc_se_same = drug_SE_Same %>%
  dplyr::select(drug, phenotype, is_side_effect, ldsc_score) %>% 
  na.omit()
data_ldsc_se_same$is_side_effect = factor(data_ldsc_se_same$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_ldsc_se_same_noOutliers = data_ldsc_se_same %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(ldsc_score, 0.25),
    Q3 = quantile(ldsc_score, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(ldsc_score >= (Q1 - 1.5 * IQR) & ldsc_score <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_ldsc_se_same %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_ldsc_se_same = ggplot(data_ldsc_se_same, aes(x = is_side_effect, y = ldsc_score, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "LDSC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_ldsc_se_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_ldsc_se_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_ldsc_se_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_ldsc_se_same_noOutliers$ldsc_score) + 0.06) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_ldsc_se_same_noOutliers$ldsc_score) + 0.1*max(data_ldsc_se_same_noOutliers$ldsc_score), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_ldsc_se_same_noOutliers$ldsc_score), max(data_ldsc_se_same_noOutliers$ldsc_score) + 0.2*max(data_ldsc_se_same_noOutliers$ldsc_score))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# MAGMA spearman correlation
data_magma_se_same = drug_SE_Same %>%
  dplyr::select(drug, phenotype, is_side_effect, magma_cor = magma_spearman_correlation) %>% 
  na.omit()
data_magma_se_same$is_side_effect = factor(data_magma_se_same$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_magma_se_same_noOutliers = data_magma_se_same %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(magma_cor, 0.25),
    Q3 = quantile(magma_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(magma_cor >= (Q1 - 1.5 * IQR) & magma_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_magma_se_same %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_magma_se_same = ggplot(data_magma_se_same, aes(x = is_side_effect, y = magma_cor, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "MAGMA",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_magma_se_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_magma_se_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_magma_se_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_magma_se_same_noOutliers$magma_cor) - 0.0125) +
  scale_y_continuous(breaks = seq(-0.03, 0.12, 0.05)) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_magma_se_same_noOutliers$magma_cor) + 0.1*max(data_magma_se_same_noOutliers$magma_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_magma_se_same_noOutliers$magma_cor), max(data_magma_se_same_noOutliers$magma_cor) + 0.2*max(data_magma_se_same_noOutliers$magma_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# S-MultiXcan spearman correlaton
data_smultixcan_se_same = drug_SE_Same %>%
  dplyr::select(drug, phenotype, is_side_effect, smultixcan_cor = smultixcan_spearman_cor) %>% 
  na.omit()
data_smultixcan_se_same$is_side_effect = factor(data_smultixcan_se_same$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_smultixcan_se_same_noOutliers = data_smultixcan_se_same %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(smultixcan_cor, 0.25),
    Q3 = quantile(smultixcan_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(smultixcan_cor >= (Q1 - 1.5 * IQR) & smultixcan_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_smultixcan_se_same %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_smultixcan_se_same = ggplot(data_smultixcan_se_same, aes(x = is_side_effect, y = smultixcan_cor, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "S-MultiXcan",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_smultixcan_se_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_smultixcan_se_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_smultixcan_se_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(-0.5, 0.1, 0.02)) +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_smultixcan_se_same_noOutliers$smultixcan_cor) - 0.025) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_smultixcan_se_same_noOutliers$smultixcan_cor) + 0.1*max(data_smultixcan_se_same_noOutliers$smultixcan_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_smultixcan_se_same_noOutliers$smultixcan_cor), max(data_smultixcan_se_same_noOutliers$smultixcan_cor) + 0.2*max(data_smultixcan_se_same_noOutliers$smultixcan_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# Open Targets L2G spearman correlation
data_l2g_se_same = drug_SE_Same %>%
  dplyr::select(drug, phenotype, is_side_effect, l2g_cor = ot_l2g_spearman_cor) %>% 
  na.omit()
data_l2g_se_same$is_side_effect = factor(data_l2g_se_same$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_l2g_se_same_noOutliers = data_l2g_se_same %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(l2g_cor, 0.25),
    Q3 = quantile(l2g_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(l2g_cor >= (Q1 - 1.5 * IQR) & l2g_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_l2g_se_same %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_l2g_se_same = ggplot(data_l2g_se_same, aes(x = is_side_effect, y = l2g_cor, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "Open Targets L2G",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_l2g_se_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_l2g_se_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_l2g_se_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_l2g_se_same_noOutliers$l2g_cor) - 0.023) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_l2g_se_same_noOutliers$l2g_cor) + 0.1*max(data_l2g_se_same_noOutliers$l2g_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_l2g_se_same_noOutliers$l2g_cor), max(data_l2g_se_same_noOutliers$l2g_cor) + 0.2*max(data_l2g_se_same_noOutliers$l2g_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# Open Targets COLOC overlap coefficient
data_coloc_se_same = drug_SE_Same %>%
  dplyr::select(drug, phenotype, is_side_effect, coloc_overcoef = ot_coloc_overlap_coefficient) %>% 
  na.omit()
data_coloc_se_same$is_side_effect = factor(data_coloc_se_same$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_coloc_se_same_noOutliers = data_coloc_se_same %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(coloc_overcoef, 0.25),
    Q3 = quantile(coloc_overcoef, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(coloc_overcoef >= (Q1 - 1.5 * IQR) & coloc_overcoef <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_coloc_se_same %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_coloc_se_same = ggplot(data_coloc_se_same, aes(x = is_side_effect, y = coloc_overcoef, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "COLOC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_coloc_se_same$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_coloc_se_same$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_coloc_se_same), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_coloc_se_same_noOutliers$coloc_overcoef) + 0.1*max(data_coloc_se_same_noOutliers$coloc_overcoef)) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_coloc_se_same_noOutliers$coloc_overcoef) + 0.07, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_coloc_se_same_noOutliers$coloc_overcoef), max(data_coloc_se_same_noOutliers$coloc_overcoef) + 0.2*max(data_coloc_se_same_noOutliers$coloc_overcoef))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")

SE_same = annotate_figure(
  ggarrange(plot_ldsc_se_same, plot_magma_se_same, plot_smultixcan_se_same, plot_l2g_se_same, plot_coloc_se_same, 
            ncol = 5, nrow = 1, align = "hv"), 
  left = text_grob(
    "Max genetic similarity between\nside effects of a drug and a test disease",
    size = 20, face = "bold", family = "Arial", rot = 90
  ),
  bottom = text_grob(
    "Is the test disease a known side effect of the drug?",
    size = 20, face = "bold", family = "Arial"
  ), 
) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

## different body system (drug indication - novel phenotype)
# LDSC score
data_ldsc_se_different = drug_SE_Different %>%
  dplyr::select(drug, phenotype, is_side_effect, ldsc_score) %>% 
  na.omit()
data_ldsc_se_different$is_side_effect = factor(data_ldsc_se_different$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_ldsc_se_different_noOutliers = data_ldsc_se_different %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(ldsc_score, 0.25),
    Q3 = quantile(ldsc_score, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(ldsc_score >= (Q1 - 1.5 * IQR) & ldsc_score <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_ldsc_se_different %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_ldsc_se_different = ggplot(data_ldsc_se_different, aes(x = is_side_effect, y = ldsc_score, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "LDSC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_ldsc_se_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_ldsc_se_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_ldsc_se_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_ldsc_se_different_noOutliers$ldsc_score) + 0.08) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_ldsc_se_different_noOutliers$ldsc_score) + 0.08, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_ldsc_se_different_noOutliers$ldsc_score), max(data_ldsc_se_different_noOutliers$ldsc_score) + 0.15)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# MAGMA spearman correlation
data_magma_se_different = drug_SE_Different %>%
  dplyr::select(drug, phenotype, is_side_effect, magma_cor = magma_spearman_correlation) %>% 
  na.omit()
data_magma_se_different$is_side_effect = factor(data_magma_se_different$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_magma_se_different_noOutliers = data_magma_se_different %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(magma_cor, 0.25),
    Q3 = quantile(magma_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(magma_cor >= (Q1 - 1.5 * IQR) & magma_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_magma_se_different %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_magma_se_different = ggplot(data_magma_se_different, aes(x = is_side_effect, y = magma_cor, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "MAGMA",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_magma_se_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_magma_se_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_magma_se_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(-0.03, 0.15, 0.05)) +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_magma_se_different_noOutliers$magma_cor) + 0.035) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_magma_se_different_noOutliers$magma_cor) + 0.033, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_magma_se_different_noOutliers$magma_cor), max(data_magma_se_different_noOutliers$magma_cor) + 0.05)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# S-MultiXcan spearman correlaton
data_smultixcan_se_different = drug_SE_Different %>%
  dplyr::select(drug, phenotype, is_side_effect, smultixcan_cor = smultixcan_spearman_cor) %>% 
  na.omit()
data_smultixcan_se_different$is_side_effect = factor(data_smultixcan_se_different$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_smultixcan_se_different_noOutliers = data_smultixcan_se_different %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(smultixcan_cor, 0.25),
    Q3 = quantile(smultixcan_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(smultixcan_cor >= (Q1 - 1.5 * IQR) & smultixcan_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_smultixcan_se_different %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_smultixcan_se_different = ggplot(data_smultixcan_se_different, aes(x = is_side_effect, y = smultixcan_cor, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "S-MultiXcan",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_smultixcan_se_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_smultixcan_se_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_smultixcan_se_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(-0.03, 0.07, by = 0.02)) + 
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_smultixcan_se_different_noOutliers$smultixcan_cor) + 0.008) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_smultixcan_se_different_noOutliers$smultixcan_cor) + 0.0155, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_smultixcan_se_different_noOutliers$smultixcan_cor), max(data_smultixcan_se_different_noOutliers$smultixcan_cor) + 0.025)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# Open Targets L2G spearman correlation
data_l2g_se_different = drug_SE_Different %>%
  dplyr::select(drug, phenotype, is_side_effect, l2g_cor = ot_l2g_spearman_cor) %>% 
  na.omit()
data_l2g_se_different$is_side_effect = factor(data_l2g_se_different$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_l2g_se_different_noOutliers = data_l2g_se_different %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(l2g_cor, 0.25),
    Q3 = quantile(l2g_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(l2g_cor >= (Q1 - 1.5 * IQR) & l2g_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_l2g_se_different %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_l2g_se_different = ggplot(data_l2g_se_different, aes(x = is_side_effect, y = l2g_cor, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "Open Targets L2G",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_l2g_se_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_l2g_se_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_l2g_se_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_l2g_se_different_noOutliers$l2g_cor) + 0.01) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_l2g_se_different_noOutliers$l2g_cor) + 0.025, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_l2g_se_different_noOutliers$l2g_cor), max(data_l2g_se_different_noOutliers$l2g_cor) + 0.05)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# Open Targets COLOC overlap coefficient
data_coloc_se_different = drug_SE_Different %>%
  dplyr::select(drug, phenotype, is_side_effect, coloc_overcoef = ot_coloc_overlap_coefficient) %>% 
  na.omit()
data_coloc_se_different$is_side_effect = factor(data_coloc_se_different$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_coloc_se_different_noOutliers = data_coloc_se_different %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(coloc_overcoef, 0.25),
    Q3 = quantile(coloc_overcoef, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(coloc_overcoef >= (Q1 - 1.5 * IQR) & coloc_overcoef <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_coloc_se_different %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_coloc_se_different = ggplot(data_coloc_se_different, aes(x = is_side_effect, y = coloc_overcoef, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "COLOC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_coloc_se_different$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_coloc_se_different$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_coloc_se_different), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_coloc_se_different_noOutliers$coloc_overcoef) + 0.014) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_coloc_se_different_noOutliers$coloc_overcoef) + 0.026, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_coloc_se_different_noOutliers$coloc_overcoef), max(data_coloc_se_different_noOutliers$coloc_overcoef) + 0.075)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")

SE_different = annotate_figure(
  ggarrange(plot_ldsc_se_different, plot_magma_se_different, plot_smultixcan_se_different, plot_l2g_se_different, plot_coloc_se_different, 
            ncol = 5, nrow = 1, align = "hv"), 
  left = text_grob(
    "Max genetic similarity between\nside effects of a drug and a test disease",
    size = 20, face = "bold", family = "Arial", rot = 90
  ),
  bottom = text_grob(
    "Is the test disease a known side effect of the drug?",
    size = 20, face = "bold", family = "Arial"
  ), 
) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# plot all together
ggarrange(
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  grobTree(rectGrob(gp = gpar(fill = "#BF3984", col = NA)),
           textGrob("Same body system", gp = gpar(col = "white", fontface = "bold", fontsize = 26))),
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  SE_same,
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  grobTree(rectGrob(gp = gpar(fill = "#FCCE25", col = NA)),
           textGrob("Different body system", gp = gpar(col = "black", fontface = "bold", fontsize = 26))),
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  SE_different,
  ncol = 1,
  heights = c(0.3, 0.55, 0.3, 7, 0.3, 0.55, 0.3, 7)
)
ggsave("figures/Supplementary_figure_6.png",
       device = "png", dpi = 600, width = 26, height = 15)
