#### Purpose: get ClinGraph embeddings for our phenotypes and use them to calculate phenotype similarity

library(dplyr)
library(data.table)
library(lsa)
library(tidyr)

# ====================== #
# Load and process data
# ====================== #

#### GWAS ID - MeSH ID mapper (links GWAS IDs to MeSH disease IDs) ####
gwas_matching_sources = fread("data/gwas_matching_between_sources_manual_final_2.txt", na.strings = "")

#### UMLS RRF file that maps UMLS CUIs to MeSH IDs
umls_rrf = data.table::fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/MRCONSO.RRF", sep = "|", quote = "")
umls_rrf = umls_rrf[, -19]
colnames(umls_rrf) = c("CUI", "LAT", "TS", "LUI", "STT", "SUI", "ISPREF", "AUI", "SAUI", "SCUI", "SDUI", "SAB", "TTY", "CODE" , "STR", "SRL", "SUPPRESS", "CVF")
umls_rrf = umls_rrf %>% filter(LAT == "ENG" & SAB == "MSH")

umls_rrf_filt = umls_rrf %>% 
  # filter for phenotypes in our sample
  filter(CODE %in% gwas_matching_sources$mesh_id) %>%
  dplyr::select(CUI, CODE) %>%
  distinct()
# manually add two missing UMLS CUIs to the closest cui
umls_rrf_filt = rbind(
  umls_rrf_filt,
  data.frame(
    "CUI" = c("C0036505", "C0007117"),
    "CODE" = c("D012626", "D002280") # matched to D012627 
  )
)

#### ClinGraph embeddings
clingraph_embeddings = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/clin_graph/ClinVec_umls.csv", header = TRUE)
clingraph_nodes = fread("/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/data/clin_graph/ClinGraph_nodes.csv")

clingraph_nodes = clingraph_nodes %>%
  filter(ntype == "UMLS_CUI") %>%
  mutate(node_id = stringr::str_remove(node_id, ":umls_cui")) %>%
  filter(node_id %in% umls_rrf_filt$CUI) %>%
  left_join(umls_rrf_filt, by = c("node_id" = "CUI")) %>%
  na.omit()
setdiff(gwas_matching_sources$mesh_id, clingraph_nodes$CODE) # all phenotypes have at least one UMLS CUI
clingraph_embeddings = clingraph_embeddings %>%
  filter(V1 %in% clingraph_nodes$node_index) %>%
  as.data.frame()
rownames(clingraph_embeddings) = clingraph_embeddings$V1
clingraph_embeddings = clingraph_embeddings[,-1]
clingraph_embeddings = as.matrix(t(clingraph_embeddings))
# cosine similarity
clingraph_cosine = cosine(clingraph_embeddings)
clingraph_cosine_df = as.data.frame(as.table(clingraph_cosine))
colnames(clingraph_cosine_df) = c("id1", "id2", "value")
# add mesh ID
clingraph_nodes$node_index = as.character(clingraph_nodes$node_index)
clingraph_cosine_df = clingraph_cosine_df %>%
  left_join(clingraph_nodes[,c("node_index","CODE")], by = c("id1" = "node_index")) %>%
  left_join(clingraph_nodes[,c("node_index","CODE")], by = c("id2" = "node_index")) %>%
  dplyr::rename(id1_mesh = CODE.x, id2_mesh = CODE.y) %>%
  # remove self-pairs (same mesh_id)
  filter(id1_mesh != id2_mesh)
# take mean similarity for mesh_ids with multiple UMLS CUI codes mapped to one MeSH code
clingraph_cosine_df = clingraph_cosine_df %>%
  group_by(id1_mesh, id2_mesh) %>%
  mutate(value = mean(value)) %>%
  ungroup() %>%
  dplyr::select(id1_mesh, id2_mesh, cosine_similarity_mean = value) %>%
  distinct()
# save
fwrite(clingraph_cosine_df, "data/ClinGraph_cosine_similarity.txt", sep = "\t")
