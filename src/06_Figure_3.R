
#### Purpose: check the relationship between genetic similarity and true drug indications and side effects

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

## for each drug-phenotype pair, take the maximum gen_sim score between the tested phenotype and all known indications of the drug
drug_IND = setDT(drug_IND)
drug_IND = drug_IND[, lapply(.SD, max, na.rm = TRUE), 
                              by = .(drug, phenotype, phenotype_category, drug_indicated_for_pheno), 
                              .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
drug_IND$drug_indicated_for_pheno = factor(drug_IND$drug_indicated_for_pheno, levels = c(0,1), labels = c(0,1))

# ================= #
# Rank-sum tests
# ================= #

# LDSC-based metric
data_ldsc_ind = drug_IND %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, ldsc_score, phenotype_category) %>% 
  na.omit()
data_ldsc_ind$drug_indicated_for_pheno = factor(data_ldsc_ind$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_ldsc_ind_noOutliers = data_ldsc_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(ldsc_score, 0.25),
    Q3 = quantile(ldsc_score, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(ldsc_score >= (Q1 - 1.5 * IQR) & ldsc_score <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_ldsc_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_ldsc_ind = ggplot(data_ldsc_ind, aes(x = drug_indicated_for_pheno, y = ldsc_score, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "LDSC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_ldsc_ind$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_ldsc_ind$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_ldsc_ind), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_ldsc_ind_noOutliers$ldsc_score) + 0.08) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_ldsc_ind_noOutliers$ldsc_score) + 0.1*max(data_ldsc_ind_noOutliers$ldsc_score), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_ldsc_ind_noOutliers$ldsc_score), max(data_ldsc_ind_noOutliers$ldsc_score) + 0.2*max(data_ldsc_ind_noOutliers$ldsc_score))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# MAGMA-based metric
data_magma_ind = drug_IND %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, magma_cor = magma_spearman_correlation, phenotype_category) %>% 
  na.omit()
data_magma_ind$drug_indicated_for_pheno = factor(data_magma_ind$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_magma_ind_noOutliers = data_magma_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(magma_cor, 0.25),
    Q3 = quantile(magma_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(magma_cor >= (Q1 - 1.5 * IQR) & magma_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_magma_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_magma_ind = ggplot(data_magma_ind, aes(x = drug_indicated_for_pheno, y = magma_cor, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "MAGMA",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_magma_ind$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_magma_ind$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_magma_ind), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_magma_ind_noOutliers$magma_cor) + 0.01*max(data_magma_ind_noOutliers$magma_cor)) +
  scale_y_continuous(breaks = seq(-0.03, 0.25, 0.05)) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_magma_ind_noOutliers$magma_cor) + 0.1*max(data_magma_ind_noOutliers$magma_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_magma_ind_noOutliers$magma_cor), max(data_magma_ind_noOutliers$magma_cor) + 0.2*max(data_magma_ind_noOutliers$magma_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# S-MultiXcan-based metric
data_smultixcan_ind = drug_IND %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, smultixcan_cor = smultixcan_spearman_cor, phenotype_category) %>% 
  na.omit()
data_smultixcan_ind$drug_indicated_for_pheno = factor(data_smultixcan_ind$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_smultixcan_ind_noOutliers = data_smultixcan_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(smultixcan_cor, 0.25),
    Q3 = quantile(smultixcan_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(smultixcan_cor >= (Q1 - 1.5 * IQR) & smultixcan_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_smultixcan_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_smultixcan_ind = ggplot(data_smultixcan_ind, aes(x = drug_indicated_for_pheno, y = smultixcan_cor, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "S-MultiXcan",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_smultixcan_ind$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_smultixcan_ind$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_smultixcan_ind), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_smultixcan_ind_noOutliers$smultixcan_cor) - 0.02) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_smultixcan_ind_noOutliers$smultixcan_cor) + 0.1*max(data_smultixcan_ind_noOutliers$smultixcan_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_smultixcan_ind_noOutliers$smultixcan_cor), max(data_smultixcan_ind_noOutliers$smultixcan_cor) + 0.2*max(data_smultixcan_ind_noOutliers$smultixcan_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# L2G-based metric
data_l2g_ind = drug_IND %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, l2g_cor = ot_l2g_spearman_cor, phenotype_category) %>% 
  na.omit()
data_l2g_ind$drug_indicated_for_pheno = factor(data_l2g_ind$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_l2g_ind_noOutliers = data_l2g_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(l2g_cor, 0.25),
    Q3 = quantile(l2g_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(l2g_cor >= (Q1 - 1.5 * IQR) & l2g_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_l2g_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_l2g_ind = ggplot(data_l2g_ind, aes(x = drug_indicated_for_pheno, y = l2g_cor, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "Open Targets L2G",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_l2g_ind$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_l2g_ind$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_l2g_ind), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_l2g_ind_noOutliers$l2g_cor) + 0.1*max(data_l2g_ind_noOutliers$l2g_cor)) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_l2g_ind_noOutliers$l2g_cor) + 0.1*max(data_l2g_ind_noOutliers$l2g_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_l2g_ind_noOutliers$l2g_cor), max(data_l2g_ind_noOutliers$l2g_cor) + 0.2*max(data_l2g_ind_noOutliers$l2g_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# COLOC-based metric
data_coloc_ind = drug_IND %>%
  dplyr::select(drug, phenotype, drug_indicated_for_pheno, coloc_overcoef = ot_coloc_overlap_coefficient, phenotype_category) %>% 
  na.omit()
data_coloc_ind$drug_indicated_for_pheno = factor(data_coloc_ind$drug_indicated_for_pheno, levels = c(0,1), labels = c("No", "Yes"))
data_coloc_ind_noOutliers = data_coloc_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  mutate(
    Q1 = quantile(coloc_overcoef, 0.25),
    Q3 = quantile(coloc_overcoef, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(coloc_overcoef >= (Q1 - 1.5 * IQR) & coloc_overcoef <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_coloc_ind %>%
  group_by(drug_indicated_for_pheno) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_coloc_ind = ggplot(data_coloc_ind, aes(x = drug_indicated_for_pheno, y = coloc_overcoef, fill = drug_indicated_for_pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "COLOC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_coloc_ind$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_coloc_ind$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_coloc_ind), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_coloc_ind_noOutliers$coloc_overcoef) + 0.12*max(data_coloc_ind_noOutliers$coloc_overcoef)) +
  geom_label(data = data_label, aes(x = drug_indicated_for_pheno, y = max(data_coloc_ind_noOutliers$coloc_overcoef) + 0.1*max(data_coloc_ind_noOutliers$coloc_overcoef), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_coloc_ind_noOutliers$coloc_overcoef), max(data_coloc_ind_noOutliers$coloc_overcoef) + 0.2*max(data_coloc_ind_noOutliers$coloc_overcoef))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")

IND = annotate_figure(
  ggarrange(plot_ldsc_ind, plot_magma_ind, plot_smultixcan_ind, plot_l2g_ind, plot_coloc_ind, 
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

### ------- SIDE EFFECTS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ####

# ====================== #
# Load and process data
# ====================== #

drug_SE = fread("data/DrugSE_TestedDisease_GenSim.txt.gz")

## for each drug-phenotype pair, take the maximum gen_sim score between the tested phenotype and all known side effects of the drug
drug_SE = setDT(drug_SE)
drug_SE = drug_SE[, lapply(.SD, max, na.rm = TRUE), 
                            by = .(drug, phenotype, phenotype_category, is_side_effect), 
                            .SDcols = ldsc_score:ot_coloc_overlap_coefficient]
drug_SE$is_side_effect = factor(drug_SE$is_side_effect, levels = c(0,1), labels = c(0,1))

# ================= #
# Rank-sum tests
# ================= #

# LDSC-based metric
data_ldsc_se = drug_SE %>%
  dplyr::select(drug, phenotype, is_side_effect, ldsc_score, phenotype_category) %>% 
  na.omit()
data_ldsc_se$is_side_effect = factor(data_ldsc_se$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_ldsc_se_noOutliers = data_ldsc_se %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(ldsc_score, 0.25),
    Q3 = quantile(ldsc_score, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(ldsc_score >= (Q1 - 1.5 * IQR) & ldsc_score <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_ldsc_se %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_ldsc_se = ggplot(data_ldsc_se, aes(x = is_side_effect, y = ldsc_score, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "LDSC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_ldsc_se$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_ldsc_se$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_ldsc_se), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_ldsc_se_noOutliers$ldsc_score) + 0.08) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_ldsc_se_noOutliers$ldsc_score) + 0.1*max(data_ldsc_se_noOutliers$ldsc_score), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_ldsc_se_noOutliers$ldsc_score), max(data_ldsc_se_noOutliers$ldsc_score) + 0.2*max(data_ldsc_se_noOutliers$ldsc_score))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# MAGMA-based metric
data_magma_se = drug_SE %>%
  dplyr::select(drug, phenotype, is_side_effect, magma_cor = magma_spearman_correlation, phenotype_category) %>% 
  na.omit()
data_magma_se$is_side_effect = factor(data_magma_se$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_magma_se_noOutliers = data_magma_se %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(magma_cor, 0.25),
    Q3 = quantile(magma_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(magma_cor >= (Q1 - 1.5 * IQR) & magma_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_magma_se %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_magma_se = ggplot(data_magma_se, aes(x = is_side_effect, y = magma_cor, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "MAGMA",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_magma_se$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_magma_se$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_magma_se), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_magma_se_noOutliers$magma_cor) - 0.0125) +
  scale_y_continuous(breaks = seq(-0.03, 0.12, 0.05)) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_magma_se_noOutliers$magma_cor) + 0.1*max(data_magma_se_noOutliers$magma_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_magma_se_noOutliers$magma_cor), max(data_magma_se_noOutliers$magma_cor) + 0.2*max(data_magma_se_noOutliers$magma_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# S-MultiXcan-based metric
data_smultixcan_se = drug_SE %>%
  dplyr::select(drug, phenotype, is_side_effect, smultixcan_cor = smultixcan_spearman_cor, phenotype_category) %>% 
  na.omit()
data_smultixcan_se$is_side_effect = factor(data_smultixcan_se$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_smultixcan_se_noOutliers = data_smultixcan_se %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(smultixcan_cor, 0.25),
    Q3 = quantile(smultixcan_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(smultixcan_cor >= (Q1 - 1.5 * IQR) & smultixcan_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_smultixcan_se %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_smultixcan_se = ggplot(data_smultixcan_se, aes(x = is_side_effect, y = smultixcan_cor, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "S-MultiXcan",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_smultixcan_se$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_smultixcan_se$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_smultixcan_se), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(-0.5, 0.1, 0.02)) +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_smultixcan_se_noOutliers$smultixcan_cor) - 0.027) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_smultixcan_se_noOutliers$smultixcan_cor) + 0.1*max(data_smultixcan_se_noOutliers$smultixcan_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_smultixcan_se_noOutliers$smultixcan_cor), max(data_smultixcan_se_noOutliers$smultixcan_cor) + 0.2*max(data_smultixcan_se_noOutliers$smultixcan_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# L2G-based metric
data_l2g_se = drug_SE %>%
  dplyr::select(drug, phenotype, is_side_effect, l2g_cor = ot_l2g_spearman_cor, phenotype_category) %>% 
  na.omit()
data_l2g_se$is_side_effect = factor(data_l2g_se$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_l2g_se_noOutliers = data_l2g_se %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(l2g_cor, 0.25),
    Q3 = quantile(l2g_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(l2g_cor >= (Q1 - 1.5 * IQR) & l2g_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_l2g_se %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_l2g_se = ggplot(data_l2g_se, aes(x = is_side_effect, y = l2g_cor, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "Open Targets L2G",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_l2g_se$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_l2g_se$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_l2g_se), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_l2g_se_noOutliers$l2g_cor) - 0.005) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_l2g_se_noOutliers$l2g_cor) + 0.1*max(data_l2g_se_noOutliers$l2g_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_l2g_se_noOutliers$l2g_cor), max(data_l2g_se_noOutliers$l2g_cor) + 0.2*max(data_l2g_se_noOutliers$l2g_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")
# COLOC-based metric
data_coloc_se = drug_SE %>%
  dplyr::select(drug, phenotype, is_side_effect, coloc_overcoef = ot_coloc_overlap_coefficient, phenotype_category) %>% 
  na.omit()
data_coloc_se$is_side_effect = factor(data_coloc_se$is_side_effect, levels = c(0,1), labels = c("No", "Yes"))
data_coloc_se_noOutliers = data_coloc_se %>%
  group_by(is_side_effect) %>%
  mutate(
    Q1 = quantile(coloc_overcoef, 0.25),
    Q3 = quantile(coloc_overcoef, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(coloc_overcoef >= (Q1 - 1.5 * IQR) & coloc_overcoef <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label = data_coloc_se %>%
  group_by(is_side_effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
plot_coloc_se = ggplot(data_coloc_se, aes(x = is_side_effect, y = coloc_overcoef, fill = is_side_effect)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "COLOC",
       subtitle = paste0(
         "Drugs: ", format(length(unique(data_coloc_se$drug)), big.mark = ","), " | ",
         "Diseases: ", format(length(unique(data_coloc_se$phenotype)), big.mark = ","), " | ",
         "Total pairs: ", format(nrow(data_coloc_se), big.mark = ",")
       )
  ) +
  xlab("") +
  ylab("") +
  stat_signif(comparisons=list(c("Yes", "No")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 5, family = "Arial", y_position = max(data_coloc_se_noOutliers$coloc_overcoef) + 0.12*max(data_coloc_se_noOutliers$coloc_overcoef)) +
  geom_label(data = data_label, aes(x = is_side_effect, y = max(data_coloc_se_noOutliers$coloc_overcoef) + 0.07, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#1B9E77")) +
  coord_cartesian(ylim = c(min(data_coloc_se_noOutliers$coloc_overcoef), max(data_coloc_se_noOutliers$coloc_overcoef) + 0.2*max(data_coloc_se_noOutliers$coloc_overcoef))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")

SE = annotate_figure(
  ggarrange(plot_ldsc_se, plot_magma_se, plot_smultixcan_se, plot_l2g_se, plot_coloc_se, 
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

## Plot IND and SE plots together
ggarrange(
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  grobTree(rectGrob(gp = gpar(fill = "#BF3984", col = NA)),
           textGrob("Drug indications", gp = gpar(col = "white", fontface = "bold", fontsize = 26))),
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  IND,
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  grobTree(rectGrob(gp = gpar(fill = "#FCCE25", col = NA)),
           textGrob("Drug side effects", gp = gpar(col = "black", fontface = "bold", fontsize = 26))),
  grobTree(rectGrob(gp = gpar(fill = "white", col = NA)),
           textGrob("")),
  SE,
  ncol = 1,
  heights = c(0.3, 0.55, 0.3, 7, 0.3, 0.55, 0.3, 7)
)
ggsave("figures/Figure_3.png",
       device = "png", dpi = 600, width = 26, height = 15)
