
### Purpose: Distribution of phenotypic similarity values among disease pairs affecting same or different body systems (Open Targets, UK HRCS, ICD10)

library(data.table)
library(dplyr)
library(DescTools)
library(ggplot2)
library(viridisLite)

## ============ ##
## Load data    ##
## ============ ##

data_AllPhenoPairs = fread("data/AllPhenoPairs_DrugOverlap_GeneticSimilarity.txt", data.table = FALSE)

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
      "Similar", "Dissimilar")) %>%
  dplyr::select(phenotype_1, phenotype_2, same_body_system_consensus_OT_HRCS_ICD10, clingraph_sim, drug_se_overlap_coef:coloc_overcoef)

data_AllPhenoPairs$same_body_system_consensus_OT_HRCS_ICD10 = factor(
  data_AllPhenoPairs$same_body_system_consensus_OT_HRCS_ICD10,
  level = c("Similar", "Dissimilar"),
  labels = c("Same", "Different")
)
ggplot(data_AllPhenoPairs, aes(x = clingraph_sim, fill = same_body_system_consensus_OT_HRCS_ICD10)) +
  # geom_histogram(bins = 200, alpha = 0.5) +
  geom_density(alpha = 0.5, color = NA) +
  scale_x_continuous(breaks = seq(-1, 1, 0.1)) +
  labs(
    title = "Distribution of ClinGraph cosine similarity values among disease pairs\naffecting same or different body systems according to Open Targets, ICD10 and UK HRCS disease labels",
    x = "ClinGraph cosine similarity",
    y = "Density",
    fill = "Body system"
  ) +
  # geom_vline(xintercept = -0.1, color = "red") +
  theme_bw(base_size = 16, base_family = "Arial") +
  # arrow for increasing similarity
  annotate(
    "segment",
    x = 0.7, xend = 1,
    y = -0.02, yend = -0.02,
    arrow = arrow(length = unit(0.25, "cm")), 
    color = "black",
    size = 1.5
  ) +
  annotate(
    "text",
    x = 0.8, y = -0.05,
    label = "Increasing phenotypic similarity",
    hjust = 0.5,
    size = 4
  ) +
  # arrow for decreasing similarity
  annotate(
    "segment",
    x = -0.7, xend = -1,
    y = -0.02, yend = -0.02,
    arrow = arrow(length = unit(0.25, "cm")), 
    color = "black",
    size = 1.5
  ) +
  annotate(
    "text",
    x = -0.8, y = -0.05,
    label = "Decreasing phenotypic similarity",
    hjust = 0.5,
    size = 4
  )
ggsave("figures/Supplementary_figure_3.png",
       device = "png", dpi = 600, width = 15, height = 8)
  
  