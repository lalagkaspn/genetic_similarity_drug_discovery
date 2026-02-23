
#### Purpose: to test whether one model can predict the other

library(data.table)
library(dplyr)
library(pROC)
library(ggplot2)
library(ggpubr)

# ======================= #
# Load and process data
# ======================= #

### NOTE: we use predictions from models trained after removing drug-phenotype pairs being both indications and side effects to exclude the possibility of inflated signal due to overlapping labels
## INDICATIONS models
data_I = fread("results/IND_STAN_AllPairs_noINDandSEoverlap_SD=2.5_OMULT=0.7_TestPred.csv")
data_I = data_I %>% dplyr::select(drug_pheno = V1, is_indicated = label, PredProb_Indicated = p)
## SIDE EFFECTS models
data_SE = fread("results/SE_STAN_AllPairs_noINDandSEoverlap_SD=3.75_OMULT=0.6_TestPred.csv")
data_SE = data_SE %>% dplyr::select(drug_pheno = V1, is_side_effect = label, PredProb_SideEffect = p)

## find drug-phenotype pairs with predicted probabilities available in both indications and side effect models
common = intersect(data_SE$drug_pheno, data_I$drug_pheno)
data_I_common = data_I %>% filter(drug_pheno %in% common)
data_SE_common = data_SE %>% filter(drug_pheno %in% common)
common_ultimate = left_join(data_I_common, data_SE_common, by = "drug_pheno") %>%
  dplyr::select(drug_pheno, is_indicated, is_side_effect, PredProb_Indicated, PredProb_SideEffect)
common_ultimate$is_indicated = factor(common_ultimate$is_indicated, levels = common_ultimate$is_indicated, labels = common_ultimate$is_indicated)
common_ultimate$is_side_effect = factor(common_ultimate$is_side_effect, levels = common_ultimate$is_side_effect, labels = common_ultimate$is_side_effect)
common_ultimate$is_ind_se = ifelse(common_ultimate$is_indicated == 1, 1, 0)
common_ultimate$is_ind_se = ifelse(common_ultimate$is_side_effect == 1, 2, common_ultimate$is_ind_se)
common_ultimate$is_ind_se = ifelse(common_ultimate$is_indicated == 1 & common_ultimate$is_side_effect == 1, 3, common_ultimate$is_ind_se)
common_ultimate$is_ind_se = factor(common_ultimate$is_ind_se, levels = c(0,3,1,2), labels = c("No Indication or Side Effect", "Indication & Side Effect", "Indication", "Side Effect"))

# ============================================================ #
# Indications model can predict side effects and vice versa
# ============================================================ #

# Indications model predicting indications
roc(common_ultimate$is_indicated ~ common_ultimate$PredProb_Indicated, levels = c(0,1), direction = "<") # 0.7403
# Indications model predicting side effects
roc(common_ultimate$is_side_effect ~ common_ultimate$PredProb_Indicated, levels = c(0,1), direction = "<") # 0.5249
# Side effects model predicting side effects
roc(common_ultimate$is_side_effect ~ common_ultimate$PredProb_SideEffect, levels = c(0,1), direction = "<") # 0.5689
# Side effects model predicting indications
roc(common_ultimate$is_indicated ~ common_ultimate$PredProb_SideEffect, levels = c(0,1), direction = "<") # 0.5424

## Prepare data frame with AUROCs
roc_list = list(
  "IND model, Pred IND" = roc(common_ultimate$is_indicated ~ common_ultimate$PredProb_Indicated, levels = c(0,1), direction = "<"),
  "IND model, Pred SE"  = roc(common_ultimate$is_side_effect ~ common_ultimate$PredProb_Indicated, levels = c(0,1), direction = "<"),
  
  "SE model, Pred IND"   = roc(common_ultimate$is_indicated ~ common_ultimate$PredProb_SideEffect, levels = c(0,1), direction = "<"),
  "SE model, Pred SE"    = roc(common_ultimate$is_side_effect ~ common_ultimate$PredProb_SideEffect, levels = c(0,1), direction = "<")
)

## Convert each to a data.frame and tag it
roc_dfs = lapply(names(roc_list), function(name) {
  roc_obj = roc_list[[name]]
  data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    label = name,
    AUC = as.numeric(auc(roc_obj))
  )
})
roc_dfs = bind_rows(roc_dfs)

## Add meta data
roc_dfs$model_general = ifelse(grepl("SE model", roc_dfs$label), "Side effects", "NA")
roc_dfs$model_general = ifelse(grepl("IND model", roc_dfs$label), "Indications", roc_dfs$model_general)
roc_dfs$model_general = factor(roc_dfs$model_general, levels = c("Indications", "Side effects"), labels = c("Indications model","Side effects model"))
roc_dfs$predicting = ifelse(grepl("Pred IND", roc_dfs$label), "Indications", "Side effects")
roc_dfs$predicting = factor(roc_dfs$predicting, levels = c("Indications", "Side effects"), labels = c("Indications", "Side effects"))
roc_dfs$auc_label = ifelse(roc_dfs$predicting == "Indications", "Predicting indications", "Predicting side effects")
roc_dfs$auc_label = factor(roc_dfs$auc_label, levels = c("Predicting indications", "Predicting side effects"), labels = c("predicting indications", "predicting side effects"))

## Plot
# Indications model predicting side effects
auc_labels_df = roc_dfs %>% 
  filter(
    model_general == "Indications model",
    predicting %in% c("Indications", "Side effects")
  ) %>%
  dplyr::select(model_general, predicting, AUC, auc_label) %>% 
  distinct() %>%
  group_by(model_general) %>%
  arrange(desc(AUC), .by_group = TRUE) %>%
  mutate(y = c(0.16, 0.1)) %>%
  ungroup()
ggplot(roc_dfs %>% filter(
  model_general == "Indications model",
  predicting %in% c("Indications", "Side effects")), aes(x = FPR, y = TPR, color = predicting, group = label)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "grey", show.legend = FALSE) +
  ggh4x::facet_grid2(
    scales = "free",
    strip = ggh4x::strip_themed(
      background_x = element_rect(fill = "#BF3984", color = "#BF3984"),
      text_x = element_text(face = "bold", color = "white", size = 18),
      background_y = element_blank(),
      text_y = element_blank()
    )) +
  scale_color_manual(values = c("#1B9E77","#D95F02")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), expand = expansion(mult = c(0, 0))) +
  geom_text(
    data = auc_labels_df,
    aes(x = 0.7, y = y, label = paste0("AUC = ", round(AUC, 3)), color = predicting),
    inherit.aes = FALSE,
    size = 5, fontface = "bold", show.legend = FALSE, hjust = 0, vjust = 0
  ) +
  labs(
    x = "1 - Specificity (FPR)",
    y = "Sensitivity (TPR)",
    color = "Prediction Target"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16, color = "black", face = "bold"),
    legend.text = element_text(size = 14, color = "black"),
    panel.spacing = unit(1, "lines"),
    axis.title.x = element_text(size = 16, color = "black", margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, color = "black", margin = margin(r = 15)),
    axis.text = element_text(size = 14, color = "black")
  )

# Side effects model predicting indications
auc_labels_df = roc_dfs %>% 
  filter(
    model_general == "Side effects model",
    predicting %in% c("Indications", "Side effects")
  ) %>%
  dplyr::select(model_general, predicting, AUC, auc_label) %>% 
  distinct() %>%
  group_by(model_general) %>%
  arrange(desc(AUC), .by_group = TRUE) %>%
  mutate(y = c(0.16, 0.1)) %>%
  ungroup()
ggplot(roc_dfs %>% filter(
    model_general == "Side effects model",
    predicting %in% c("Indications", "Side effects")), aes(x = FPR, y = TPR, color = predicting, group = label)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "grey", show.legend = FALSE) +
  ggh4x::facet_grid2(
    scales = "free",
    strip = ggh4x::strip_themed(
      background_x = element_rect(fill = "#BF3984", color = "#BF3984"),
      text_x = element_text(face = "bold", color = "white", size = 18),
      background_y = element_blank(),
      text_y = element_blank()
    )) +
  scale_color_manual(values = c("#1B9E77","#D95F02")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), expand = expansion(mult = c(0, 0))) +
  geom_text(
    data = auc_labels_df,
    aes(x = 0.7, y = y, label = paste0("AUC = ", round(AUC, 3)), color = predicting),
    inherit.aes = FALSE,
    size = 5, fontface = "bold", show.legend = FALSE, hjust = 0, vjust = 0
  ) +
  labs(
    x = "1 - Specificity (FPR)",
    y = "Sensitivity (TPR)",
    color = "Prediction Target"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16, color = "black", face = "bold"),
    legend.text = element_text(size = 14, color = "black"),
    panel.spacing = unit(1, "lines"),
    axis.title.x = element_text(size = 16, color = "black", margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, color = "black", margin = margin(r = 15)),
    axis.text = element_text(size = 14, color = "black")
  )

# ============================================================================================ #
# Permutations to rule out the possibility that the above results are observed just by chance
# ============================================================================================ #

## Indications model predicting side effects
IND_PredSE_PermAUROC = c()
for (i in 1:1000) {
  data_temp = common_ultimate
  data_temp$is_side_effect = sample(data_temp$is_side_effect, replace = FALSE)
  IND_PredSE_PermAUROC = c(
    IND_PredSE_PermAUROC,
    roc(data_temp$is_side_effect ~ data_temp$PredProb_Indicated, levels = c(0, 1), direction = "<")$auc
  )
  cat(i, "- 1000", "\r")
}
## Side effects model predicting indications
SE_PredIND_PermAUROC = c()
for (i in 1:1000) {
  data_temp = common_ultimate
  data_temp$is_indicated = sample(data_temp$is_indicated, replace = FALSE)
  SE_PredIND_PermAUROC = c(
    SE_PredIND_PermAUROC,
    roc(data_temp$is_indicated ~ data_temp$PredProb_SideEffect, levels = c(0, 1), direction = "<")$auc
  )
  cat(i, "- 1000", "\r")
}

## Merge results
cross_model_permAUROC = data.frame(
  permutation = paste0("Permutation_", 1:1000),
  IND_PredSE = IND_PredSE_PermAUROC,
  SE_PredIND = SE_PredIND_PermAUROC
)
cross_model_permAUROC = reshape2::melt(cross_model_permAUROC, "permutation", colnames(cross_model_permAUROC[2:ncol(cross_model_permAUROC)]))
cross_model_permAUROC$model_name = rep(c("INDICATIONS model predicting SIDE EFFECTS", "SIDE EFFECTS model predicting INDICATIONS"), each = 1000)
cross_model_permAUROC$model_name = factor(cross_model_permAUROC$model_name, levels = cross_model_permAUROC$model_name, labels = cross_model_permAUROC$model_name)
cross_model_permAUROC$value_observed_AUROC = rep(
  c(
    roc_dfs %>% filter(label == "IND model, Pred SE") %>% pull(AUC) %>% unique(), # IND model - predicting SE
    roc_dfs %>% filter(label == "SE model, Pred IND") %>% pull(AUC) %>% unique()  # SE model - predicting IND
  ),
  each = 1000)
cross_model_permAUROC$label_observedAUROC = "Observed"
cross_model_permAUROC$label_permutedAUROC = "Permuted"
cross_model_permAUROC$permutation_pvalue = c(
  rep(sum(roc_dfs %>% filter(label == "IND model, Pred SE") %>% pull(AUC) %>% unique() <= IND_PredSE_PermAUROC) / length(IND_PredSE_PermAUROC), 1000),  # IND model - SAME - predicting SE
  rep(sum(roc_dfs %>% filter(label == "SE model, Pred IND") %>% pull(AUC) %>% unique() <= SE_PredIND_PermAUROC) / length(SE_PredIND_PermAUROC), 1000)   # SE model - DIFFERENT - predicting IND
)
cross_model_permAUROC$permutation_pvalue = ifelse(
  cross_model_permAUROC$permutation_pvalue == 0,
  "p[perm]*'<0.001'",
  paste0("p[perm]*'=",cross_model_permAUROC$permutation_pvalue, "'")
)

## Plot
## NOTE: Due to the stochastic nature of the permutations, the code below will generate a figure slightly different than the one
## in the paper. Howeover, the conclusions should be the same.

# Indications model predicting side effects
ind_pred_se_permutations =
  ggplot(cross_model_permAUROC %>% filter(variable == "IND_PredSE"),
         aes(x = value)) +
  geom_histogram(fill = "#DDEAF2", color = "#3A506B", linewidth = 0.4) +
  geom_vline(aes(xintercept = value_observed_AUROC), color = "#D7263D", linetype = "dashed", linewidth = 1) +
  geom_text(
    aes(x = rep(0.524, 1000), 
        y = rep(140, 1000), 
        label = label_observedAUROC),
    size = 5, color = "#D7263D", hjust = 1) +
  geom_text(
    aes(x = rep(0.5, 1000), 
        y = rep(140, 1000),
        label = label_permutedAUROC),
    size = 5, color = "#3A506B") +
  geom_text(
    aes(x = rep(0.524, 1000),
        y = rep(5, 1000),
        label = permutation_pvalue),
    parse = TRUE, size = 5, color = "black", hjust = 1
  ) +
  labs(
    title = "INDICATIONS model predicts SIDE EFFECTS",
    x = "AUROC",
    y = "Count"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.01), labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(breaks = seq(0, 150, 25)) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.title.x = element_text(size = 16, color = "black", margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, color = "black", margin = margin(r = 15)),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )

## Side effects model predicting indications
se_pred_ind_permutations =
  ggplot(cross_model_permAUROC %>% filter(variable == "SE_PredIND"),
         aes(x = value)) +
  geom_histogram(fill = "#DDEAF2", color = "#3A506B", linewidth = 0.4) +
  geom_vline(aes(xintercept = value_observed_AUROC), color = "#D7263D", linetype = "dashed", linewidth = 1) +
  geom_text(
    aes(x = rep(0.541, 1000), 
        y = rep(125, 1000), 
        label = label_observedAUROC),
    size = 5, color = "#D7263D", hjust = 1) +
  geom_text(
    aes(x = rep(0.5, 1000), 
        y = rep(125, 1000),
        label = label_permutedAUROC),
    size = 5, color = "#3A506B") +
  geom_text(
    aes(x = rep(0.541, 1000),
        y = rep(5, 1000),
        label = permutation_pvalue),
    parse = TRUE, size = 5, color = "black", hjust = 1
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.02), labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(breaks = seq(0, 150, 25)) +
  labs(
    title = "SIDE EFFECTS model predicts INDICATIONS",
    x = "AUROC",
    y = "Count"
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.title.x = element_text(size = 16, color = "black", margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, color = "black", margin = margin(r = 15)),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )

ggarrange(
  ind_pred_se_permutations,
  se_pred_ind_permutations,
  ncol = 2, nrow = 1,
  labels = "AUTO",
  font.label = list(size = 30)
)
ggsave("figures/Figure_6.png",
       device = "png", dpi = 600, width = 20, height = 6)
