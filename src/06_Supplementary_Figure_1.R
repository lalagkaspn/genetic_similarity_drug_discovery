
#### Purpose: to check whether the developed genetic similarity and drug overlap metrics capture known biology
## Group categories based on Open Targets, ICD10 and UK HRCS disease labels (no clingraph similarity here!)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)

## ============ ##
## Load data    ##
## ============ ##

data_AllPhenoPairs = fread("data/AllPhenoPairs_DrugOverlap_GeneticSimilarity.txt")

### add phenotype annotations
# ClinGraph
clingraph_cosine_df = fread("data/ClinGraph_cosine_similarity.txt")
# UK HRCS & ICD10
trait_categories = fread("data/phenotype_categories.txt") %>%
  dplyr::select(mesh_id, trait_area_HRCS, trait_ICD10_top)
data_AllPhenoPairs = data_AllPhenoPairs %>%
  left_join(trait_categories, by = c("phenotype_1" = "mesh_id")) %>%
  left_join(trait_categories, by = c("phenotype_2" = "mesh_id")) %>%
  left_join(clingraph_cosine_df, by = c("phenotype_1" = "id1_mesh", "phenotype_2" = "id2_mesh")) %>%
  dplyr::select(phenotype_1, phenotype_2, drug_se_overlap_coef, drug_ind_overlap_coef, ldsc_score, magma_cor, smultixcan_cor, l2g_cor, coloc_overcoef,
                phenotype_1_category_OT = phenotype_1_category, phenotype_2_category_OT = phenotype_2_category,
                phenotype_1_category_HRCS = trait_area_HRCS.x, phenotype_2_category_HRCS = trait_area_HRCS.y,
                phenotype_1_category_ICD10 = trait_ICD10_top.x, phenotype_2_category_ICD10 = trait_ICD10_top.y,
                clingraph_sim = cosine_similarity_mean) %>%
  mutate(
    same_body_system_consensus_OT_HRCS_ICD10 = if_else(
      phenotype_1_category_OT == phenotype_2_category_OT | phenotype_1_category_HRCS == phenotype_2_category_HRCS | phenotype_1_category_ICD10 == phenotype_2_category_ICD10,
      "Same", "Different")) %>%
  dplyr::select(phenotype_1, phenotype_2, same_body_system_consensus_OT_HRCS_ICD10, clingraph_sim, drug_se_overlap_coef:coloc_overcoef)

## Use only disease annotations (OT, UK HRCS, ICD-10)
data_AllPhenoPairs$phenotype_category = factor(
  data_AllPhenoPairs$same_body_system_consensus_OT_HRCS_ICD10,
  levels = c("Same","Different"),
  labels = c("Same","Different")
)

# ================================================================ #
# Genetic similarity of phenotype pairs stratified by body system  #
# ================================================================ #

### For visualization reasons only, outliers, defined as as values below Q1 − 1.5×IQR or above Q3 + 1.5×IQR, are removed (but included in all statistical analyses)

## LDSC-based metric
data_ldsc = data_AllPhenoPairs %>% 
  dplyr::select(phenotype_1, phenotype_2, phenotype_category, ldsc_score) %>% 
  na.omit() %>%
  distinct()
data_ldsc_noOutliers = data_ldsc %>%
  group_by(phenotype_category) %>%
  mutate(
    Q1 = quantile(ldsc_score, 0.25),
    Q3 = quantile(ldsc_score, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(ldsc_score >= (Q1 - 1.5 * IQR) & ldsc_score <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_ldsc_label = data_ldsc %>%
  group_by(phenotype_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
gg_ldsc_BodySystem = ggplot(data_ldsc, aes(x = phenotype_category, y = ldsc_score, fill = phenotype_category)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "LDSC", 
       subtitle = paste0(
         "Phenotypes: ", length(unique(c(data_ldsc$phenotype_1, data_ldsc$phenotype_2))), " (", format(nrow(data_ldsc), big.mark = ","), " pairs)"
       ),
       fill = "Body system"
  ) +
  xlab("Body system") +
  ylab("Genome-wide correlation (Rg)") +
  scale_y_continuous(breaks = seq(-0.5, 0.6, 0.1)) +
  scale_fill_manual(values = c("Same" = "#BF3984", "Different" = "#FCCE25")) +
  stat_signif(comparisons=list(c("Same", "Different")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 6, family = "Arial", y_position = max(data_ldsc_noOutliers$ldsc_score) + 0.02) +
  geom_label(data = data_ldsc_label, aes(x = phenotype_category, y = max(data_ldsc_noOutliers$ldsc_score) + 0.10*max(data_ldsc_noOutliers$ldsc_score), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  coord_cartesian(ylim = c(min(data_ldsc_noOutliers$ldsc_score), max(data_ldsc_noOutliers$ldsc_score) + 0.25*max(data_ldsc_noOutliers$ldsc_score))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold", margin = margin(t = 1)),
        axis.title.y = element_text(size = 16, color = "black", face = "bold", margin = margin(r = 1)),
        legend.position = "none")

## MAGMA-based metric
data_magma = data_AllPhenoPairs %>% 
  dplyr::select(phenotype_1, phenotype_2, phenotype_category, magma_cor) %>% 
  na.omit() %>%
  distinct()
data_magma_noOutliers = data_magma %>%
  group_by(phenotype_category) %>%
  mutate(
    Q1 = quantile(magma_cor, 0.25),
    Q3 = quantile(magma_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(magma_cor >= (Q1 - 1.5 * IQR) & magma_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_magma_label = data_magma %>%
  group_by(phenotype_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
gg_magma_BodySystem = ggplot(data_magma, aes(x = phenotype_category, y = magma_cor, fill = phenotype_category)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "MAGMA", 
       subtitle = paste0(
         "Phenotypes: ", length(unique(c(data_magma$phenotype_1, data_magma$phenotype_2))), " (", format(nrow(data_magma), big.mark = ","), " pairs)"
       ),
       fill = "Body system"
  ) +
  xlab("Body system") +
  ylab("Gene-based similarity\n(spearman correlation)") +
  scale_y_continuous(breaks = seq(-0.04, 0.1, 0.02)) +
  scale_fill_manual(values = c("Same" = "#BF3984", "Different" = "#FCCE25")) +
  stat_signif(comparisons=list(c("Same", "Different")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 6, family = "Arial", y_position = max(data_magma_noOutliers$magma_cor) - 0.017) +
  geom_label(data = data_magma_label, aes(x = phenotype_category, y = max(data_magma_noOutliers$magma_cor) + 0.12*max(data_magma_noOutliers$magma_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  coord_cartesian(ylim = c(min(data_magma_noOutliers$magma_cor), max(data_magma_noOutliers$magma_cor) + 0.25*max(data_magma_noOutliers$magma_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold", margin = margin(t = 1)),
        axis.title.y = element_text(size = 16, color = "black", face = "bold", margin = margin(r = 1)),
        legend.position = "none")

## S-MultiXcan-based metric
data_smultixcan = data_AllPhenoPairs %>% 
  dplyr::select(phenotype_1, phenotype_2, phenotype_category, smultixcan_cor) %>% 
  na.omit() %>%
  distinct()
data_smultixcan_noOutliers = data_smultixcan  %>% 
  group_by(phenotype_category) %>%
  mutate(
    Q1 = quantile(smultixcan_cor, 0.25),
    Q3 = quantile(smultixcan_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(smultixcan_cor >= (Q1 - 1.5 * IQR) & smultixcan_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_smultixcan_label = data_smultixcan %>%
  group_by(phenotype_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
gg_smultixcan_BodySystem = ggplot(data_smultixcan, aes(x = phenotype_category, y = smultixcan_cor, fill = phenotype_category)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "S-MultiXcan", 
       subtitle = paste0(
         "Phenotypes: ", length(unique(c(data_smultixcan$phenotype_1, data_smultixcan$phenotype_2))), " (", format(nrow(data_smultixcan), big.mark = ","), " pairs)"
       ),
       fill = "Body system"
  ) +
  xlab("Body system") +
  ylab("Genetic regulation similarity\n(spearman correlation)") +
  scale_y_continuous(breaks = seq(-0.04, 0.05, 0.01)) +
  scale_fill_manual(values = c("Same" = "#BF3984", "Different" = "#FCCE25")) +
  stat_signif(comparisons=list(c("Same", "Different")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 6, family = "Arial", y_position = max(data_smultixcan_noOutliers$smultixcan_cor) - 0.033) +
  geom_label(data = data_smultixcan_label, aes(x = phenotype_category, y = max(data_smultixcan_noOutliers$smultixcan_cor) + 0.11*max(data_smultixcan_noOutliers$smultixcan_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  coord_cartesian(ylim = c(min(data_smultixcan_noOutliers$smultixcan_cor), max(data_smultixcan_noOutliers$smultixcan_cor) + 0.25*max(data_smultixcan_noOutliers$smultixcan_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold", margin = margin(t = 1)),
        axis.title.y = element_text(size = 16, color = "black", face = "bold", margin = margin(r = 1)),
        legend.position = "none")

## L2G-based metric
data_l2g = data_AllPhenoPairs %>% 
  dplyr::select(phenotype_1, phenotype_2, phenotype_category, l2g_cor) %>%
  na.omit() %>%
  distinct()
data_l2g_noOutliers = data_l2g %>% 
  group_by(phenotype_category) %>%
  mutate(
    Q1 = quantile(l2g_cor, 0.25),
    Q3 = quantile(l2g_cor, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(l2g_cor >= (Q1 - 1.5 * IQR) & l2g_cor <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_l2g_label = data_l2g %>%
  group_by(phenotype_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
gg_l2g_BodySystem = ggplot(data_l2g, aes(x = phenotype_category, y = l2g_cor, fill = phenotype_category)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "Open Targets L2G", 
       subtitle = paste0(
         "Phenotypes: ", length(unique(c(data_l2g$phenotype_1, data_l2g$phenotype_2))), " (", format(nrow(data_l2g), big.mark = ","), " pairs)"
       ),
       fill = "Body system"
  ) +
  xlab("Body system") +
  ylab("L2G score similarity\n(spearman correlation)") +
  scale_y_continuous(breaks = seq(-0.04, 0.05, 0.01)) +
  scale_fill_manual(values = c("Same" = "#BF3984", "Different" = "#FCCE25")) +
  stat_signif(comparisons=list(c("Same", "Different")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 6, family = "Arial", y_position = max(data_l2g_noOutliers$l2g_cor) - 0.04) +
  geom_label(data = data_l2g_label, aes(x = phenotype_category, y = max(data_l2g_noOutliers$l2g_cor) + 0.10*max(data_l2g_noOutliers$l2g_cor), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  coord_cartesian(ylim = c(min(data_l2g_noOutliers$l2g_cor), max(data_l2g_noOutliers$l2g_cor) + 0.25*max(data_l2g_noOutliers$l2g_cor))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold", margin = margin(t = 1)),
        axis.title.y = element_text(size = 16, color = "black", face = "bold", margin = margin(r = 1)),
        legend.position = "none")

## COLOC-based metric
data_coloc = data_AllPhenoPairs %>% 
  dplyr::select(phenotype_1, phenotype_2, phenotype_category, coloc_overcoef) %>% 
  na.omit() %>%
  distinct()
data_coloc_noOutliers = data_coloc %>% 
  group_by(phenotype_category) %>%
  mutate(
    Q1 = quantile(coloc_overcoef, 0.25),
    Q3 = quantile(coloc_overcoef, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(coloc_overcoef >= (Q1 - 1.5 * IQR) & coloc_overcoef <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label_coloc = data_coloc %>%
  group_by(phenotype_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
# boxplot
data_label_coloc[1, "label"] = sprintf("N[\"pairs same\"] == \"%s\"", format(as.numeric(data_label_coloc[1,"n"]), big.mark = ","))
data_label_coloc[2, "label"] = sprintf("N[\"pairs different\"] == \"%s\"", format(as.numeric(data_label_coloc[2,"n"]), big.mark = ","))
gg_coloc_BodySystem = ggplot(data_coloc, aes(x = phenotype_category, y = coloc_overcoef, fill = phenotype_category)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "COLOC", 
       subtitle = paste0(
         "Phenotypes: ", length(unique(c(data_coloc$phenotype_1, data_coloc$phenotype_2))), " (", format(nrow(data_coloc), big.mark = ","), " pairs)"
       ),
       fill = "Body system"
  ) +
  xlab("Body system") +
  ylab("Similarity of colocalized GWAS-molQTL signals\n(overlap coefficient)") +
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  scale_fill_manual(values = c("Same" = "#BF3984", "Different" = "#FCCE25")) +
  stat_signif(comparisons=list(c("Same", "Different")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 6, family = "Arial", y_position = max(data_coloc_noOutliers$coloc_overcoef) - 0.026) +
  geom_label(data = data_label_coloc, aes(x = phenotype_category, y = max(data_coloc_noOutliers$coloc_overcoef) + 0.1*max(data_coloc_noOutliers$coloc_overcoef), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  coord_cartesian(ylim = c(min(data_coloc_noOutliers$coloc_overcoef), max(data_coloc_noOutliers$coloc_overcoef) + 0.23*max(data_coloc_noOutliers$coloc_overcoef))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold", margin = margin(t = 1)),
        axis.title.y = element_text(size = 16, color = "black", face = "bold", margin = margin(r = 1)),
        legend.position = "none")

# ========================================================== #
# Drug overlap of phenotype pairs stratified by body system
# ========================================================== #

## Drug indications
data_drugIND = data_AllPhenoPairs %>%
  dplyr::select(phenotype_1, phenotype_2, phenotype_category, drug_ind_overlap_coef) %>%
  distinct() %>%
  na.omit()
data_drugIND_noOutliers = data_drugIND %>%
  group_by(phenotype_category) %>%
  mutate(
    Q1 = quantile(drug_ind_overlap_coef, 0.25),
    Q3 = quantile(drug_ind_overlap_coef, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(drug_ind_overlap_coef >= (Q1 - 1.5 * IQR) & drug_ind_overlap_coef <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label_drugIND = data_drugIND %>%
  group_by(phenotype_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
# boxplot
data_label_drugIND[1, "label"] = sprintf("N[\"pairs same\"] == \"%s\"", format(as.numeric(data_label_drugIND[1,"n"]), big.mark = ","))
data_label_drugIND[2, "label"] = sprintf("N[\"pairs different\"] == \"%s\"", format(as.numeric(data_label_drugIND[2,"n"]), big.mark = ","))
gg_drugIND = ggplot(data_drugIND, aes(x = phenotype_category, y = drug_ind_overlap_coef, fill = phenotype_category)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "Drug Indications", 
       subtitle = paste0(
         "Phenotypes: ", length(unique(c(data_drugIND$phenotype_1, data_drugIND$phenotype_2))), " (", format(nrow(data_drugIND), big.mark = ","), " pairs)"
       ),
       fill = "Body system"
  ) +
  xlab("Body system") +
  ylab("Similarity of indicated drugs\n(overlap coefficient)") +
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  scale_fill_manual(values = c("Same" = "#BF3984", "Different" = "#FCCE25")) +
  stat_signif(comparisons=list(c("Same", "Different")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 6, family = "Arial", y_position = max(data_drugIND_noOutliers$drug_ind_overlap_coef)) +
  geom_label(data = data_label_drugIND, aes(x = phenotype_category, y = max(data_drugIND_noOutliers$drug_ind_overlap_coef) + 0.12*max(data_drugIND_noOutliers$drug_ind_overlap_coef), label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  coord_cartesian(ylim = c(min(data_drugIND_noOutliers$drug_ind_overlap_coef), max(data_drugIND_noOutliers$drug_ind_overlap_coef) + 0.23*max(data_drugIND_noOutliers$drug_ind_overlap_coef))) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold", margin = margin(t = 1)),
        axis.title.y = element_text(size = 16, color = "black", face = "bold", margin = margin(r = 1)),
        legend.position = "none")

## Drug side effects
data_drugSE = data_AllPhenoPairs %>% 
  dplyr::select(phenotype_1, phenotype_2, phenotype_category, drug_se_overlap_coef) %>%
  distinct() %>% 
  na.omit()
data_drugSE_noOutliers = data_drugSE %>%
  group_by(phenotype_category) %>%
  mutate(
    Q1 = quantile(drug_se_overlap_coef, 0.25),
    Q3 = quantile(drug_se_overlap_coef, 0.75),
    IQR = Q3 - Q1
  ) %>%
  filter(drug_se_overlap_coef >= (Q1 - 1.5 * IQR) & drug_se_overlap_coef <= (Q3 + 1.5 * IQR)) %>%
  ungroup()
data_label_drugSE = data_drugSE %>%
  group_by(phenotype_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = sprintf("N[\"pairs\"] == \"%s\"", format(n, big.mark = ",")))
gg_drugSE = ggplot(data_drugSE, aes(x = phenotype_category, y = drug_se_overlap_coef, fill = phenotype_category)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(title = "Drug Side Effects", 
       subtitle = paste0(
         "Phenotypes: ", length(unique(c(data_drugSE$phenotype_1, data_drugSE$phenotype_2))), " (", format(nrow(data_drugSE), big.mark = ","), " pairs)"
       ),
       fill = "Body system"
  ) +
  xlab("Body system") +
  ylab("Similarity of side effect drugs\n(overlap coefficient)") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = c("Same" = "#BF3984", "Different" = "#FCCE25")) +
  stat_signif(comparisons=list(c("Same", "Different")), test="wilcox.test", test.args=list(alternative = "two.sided"),
              map_signif_level = function(p) sprintf("p = %.2g", p),
              tip_length = 0,  textsize = 6, family = "Arial", y_position = max(data_drugSE_noOutliers$drug_se_overlap_coef) + 0.09) +
  geom_label(data = data_label_drugSE, aes(x = phenotype_category, y = max(data_drugSE_noOutliers$drug_se_overlap_coef) + 0.07, label = label),
             fill = "white", size = 5, label.size = 0.5, parse = TRUE) +
  coord_cartesian(ylim = c(min(data_drugSE_noOutliers$drug_se_overlap_coef), max(data_drugSE_noOutliers$drug_se_overlap_coef) + 0.16)) +
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold", margin = margin(t = 1)),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        legend.position = "none")

# ============== #
# Combined plot
# ============== #

legend = ggplot(data_drugSE, aes(x = phenotype_category, y = drug_se_overlap_coef, color = phenotype_category)) +
  geom_point(size = 6) +
  theme_void() +
  labs(fill = "Body system", color = "Body system") +
  scale_color_manual(values = c("Same" = "#BF3984", "Different" = "#FCCE25")) +
  scale_x_discrete(expand = c(-0.1,-0.1)) +
  scale_y_continuous(limits = c(-0.1,-0.1)) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = c(0.5,0.5),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size =  23),
    legend.title = element_text(size = 25, face = "bold")
  )
ggarrange(
  gg_ldsc_BodySystem, gg_magma_BodySystem, gg_smultixcan_BodySystem, gg_l2g_BodySystem,
  gg_coloc_BodySystem, gg_drugIND, gg_drugSE, legend, 
  ncol = 4, nrow = 2, align = "hv", labels = list("A","B","C","D","E","F","G",""),
  font.label = list(size = 25, face = "bold", family = "Arial", color = "black")
)
ggsave(filename = "figures/Supplementary_figure_1.png",
       device = "png", dpi = 600, width = 24, height = 12)
