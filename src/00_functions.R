#### Overlap coefficient ####
overlap_coeff = function(a, b) {
  intersection = length(intersect(a, b))
  min_vector_size = min(length(a), length(b))
  return (intersection / min_vector_size)
}

#### Apply a function pairwise ####
pairwise_overlap = function(data, overlap_fn) {
  studies = data$left_study
  results = data.frame()
  combn(1:nrow(data), 2, function(idx) {
    study1 = studies[idx[1]]
    study2 = studies[idx[2]]
    genes1 = unlist(data$genes[idx[1]])
    genes2 = unlist(data$genes[idx[2]])
    overlap_score = overlap_fn(genes1, genes2)
    # Store the results
    results <<- rbind(results, data.frame(
      study1 = study1,
      study2 = study2,
      overlap = overlap_score
    ))
  })
  return(results)
}

#### Get spearman correlation p-value ####
spearman_cor_pvalue = function(data, metric, drug_effect, body_system) {
  if (drug_effect == "indication") {
    if (body_system == "Same") {
      data = data %>% filter(phenotype_category == "Same")
      return(cor.test(data[,grepl(metric, colnames(data))], data$drug_ind_overlap_coef,
                      use = "complete.obs", method = "spearman")$p.value)
    }
    if (body_system == "Different") {
      data = data %>% filter(phenotype_category == "Different")
      return(cor.test(data[,grepl(metric, colnames(data))], data$drug_ind_overlap_coef,
                      use = "complete.obs", method = "spearman")$p.value)
    }
    if (body_system == "All") {
      return(cor.test(data[,grepl(metric, colnames(data))], data$drug_ind_overlap_coef,
                      use = "complete.obs", method = "spearman")$p.value)
    }
  }
  if (drug_effect == "side effect") {
    if (body_system == "Same") {
      data = data %>% filter(phenotype_category == "Same")
      return(cor.test(data[,grepl(metric, colnames(data))], data$drug_se_overlap_coef,
                      use = "complete.obs", method = "spearman")$p.value)
    }
    if (body_system == "Different") {
      data = data %>% filter(phenotype_category == "Different")
      return(cor.test(data[,grepl(metric, colnames(data))], data$drug_se_overlap_coef,
                      use = "complete.obs", method = "spearman")$p.value)
    }
    if (body_system == "All") {
      return(cor.test(data[,grepl(metric, colnames(data))], data$drug_se_overlap_coef,
                      use = "complete.obs", method = "spearman")$p.value)
    }
  }
}

#### Get speraman correlation confidence intervals ####
spearman_cor_ci = function(data, metric, drug_effect, body_system, conf_int) {
  if (drug_effect == "indication") {
    if (body_system == "Same") {
      data = data %>% filter(phenotype_category == "Same")
      return(SpearmanRho(data[,grepl(metric, colnames(data))], data$drug_ind_overlap_coef,
                         use = "complete.obs", conf.level = conf_int))
    }
    if (body_system == "Different") {
      data = data %>% filter(phenotype_category == "Different")
      return(SpearmanRho(data[,grepl(metric, colnames(data))], data$drug_ind_overlap_coef,
                         use = "complete.obs", conf.level = conf_int))
    }
    if (body_system == "All") {
      return(SpearmanRho(data[,grepl(metric, colnames(data))], data$drug_ind_overlap_coef,
                         use = "complete.obs", conf.level = conf_int))
    }
  }
  if (drug_effect == "side effect") {
    if (body_system == "Same") {
      data = data %>% filter(phenotype_category == "Same")
      return(SpearmanRho(data[,grepl(metric, colnames(data))], data$drug_se_overlap_coef,
                         use = "complete.obs", conf.level = conf_int))
    }
    if (body_system == "Different") {
      data = data %>% filter(phenotype_category == "Different")
      return(SpearmanRho(data[,grepl(metric, colnames(data))], data$drug_se_overlap_coef,
                         use = "complete.obs", conf.level = conf_int))
    }
    if (body_system == "All") {
      return(SpearmanRho(data[,grepl(metric, colnames(data))], data$drug_se_overlap_coef,
                         use = "complete.obs", conf.level = conf_int))
    }
  }
}

#### Add number of disease pairs tested in each spearman correlation analysis ####
nr_pheno_pairs_tested = function(data, body_system, metric, drug_effect) {
  if (drug_effect == "indication") {
    if (body_system == "All") {
      return(
        data %>% 
          dplyr::select(phenotype_1, phenotype_2, drug_ind_overlap_coef, grep(metric, colnames(data))) %>%
          na.omit() %>%
          nrow()
      )
    }
    if (body_system == "Same") {
      return(
        data %>% filter(phenotype_category == "Same") %>% 
          dplyr::select(phenotype_1, phenotype_2, drug_ind_overlap_coef, grep(metric, colnames(data))) %>%
          na.omit() %>%
          nrow()
      )
    }
    if (body_system == "Different") {
      return(
        data %>% 
          filter(phenotype_category == "Different") %>% 
          dplyr::select(phenotype_1, phenotype_2, drug_ind_overlap_coef, grep(metric, colnames(data))) %>% 
          na.omit() %>%
          nrow()
      )
    }
  }
  if (drug_effect == "side effect") {
    if (body_system == "All") {
      return(
        data %>% 
          dplyr::select(phenotype_1, phenotype_2, drug_se_overlap_coef, grep(metric, colnames(data))) %>%
          na.omit() %>%
          nrow()
      )
    }
    if (body_system == "Same") {
      return(
        data %>% filter(phenotype_category == "Same") %>% 
          dplyr::select(phenotype_1, phenotype_2, drug_se_overlap_coef, grep(metric, colnames(data))) %>%
          na.omit() %>%
          nrow()
      )
    }
    if (body_system == "Different") {
      return(
        data %>% 
          filter(phenotype_category == "Different") %>% 
          dplyr::select(phenotype_1, phenotype_2, drug_se_overlap_coef, grep(metric, colnames(data))) %>% 
          na.omit() %>%
          nrow()
      )
    }
  }
}

#### estimate 95% Katz CIs for the relative success analysis ####
RS_ci95 = function(a,b,c,d) {
  # a: high similarity, success
  # b: high similarity, failure (total - success)
  # c: low similarity, success
  # d: low similarity, failure (total - success)
  RS <- (a / (a + b)) / (c / (c + d))
  SE_logRS = sqrt(1/a - 1/(a + b) + 1/c - 1/(c + d))
  CI = exp(log(RS) + c(-1.96, 1.96) * SE_logRS)
  return(c(RS, CI))
}

#### stan indications ####
library(pROC)
mult_split_eval_IND = function(di, subsel, sd_manuals, omults,di_pairs_in=c()){
  evidences = colnames(di)[7:ncol(di)]
  di$pair_names = paste(di$drug, di$phenotype, sep=":")
  mod_mu = stan_model('/project/pi_rachel_melamed_uml_edu/Panos/CD_CD_shared_genetics_therapeutics/genetic_similarity_drug_discovery/src/06.12_multi_evidence.stan',save_dso=T)
  
  di_pairs = unique(di$pair_names) 
  nsplit = 2
  if(length(di_pairs_in) > 0){
    di_pairs = di_pairs_in
    cat("using saved di_pairs",di_pairs_in[1:5],"\n")
  }else if(subsel > 0){
    set.seed(11042024) # to always get the same pairs
    di_pairs = sample(di_pairs, size=subsel)
  }
  di_subsel = di[di$pair_names %in% di_pairs,]
  N_di = length(di_pairs) 
  pair_id = 1:N_di
  names(pair_id) = di_pairs
  di_subsel$drug_ind_index = pair_id[di_subsel$pair_names]
  
  evidence_top = c()
  index_top = c()
  size_top = c()
  
  evidence_other = c()
  index_other = c()
  size_other = c()
  
  
  for(ev in evidences){
    
    lmdf = di_subsel[,c("drug_ind_index",ev)] #,'drug_indication_category','phenotype_category','body_system')]
    lmdf = lmdf[!is.na(lmdf[,ev]),]
    lmdf = lmdf[!(lmdf[,ev]==0),]
    lmdf = lmdf[order(lmdf[,ev], decreasing=T),]
    
    other_df = lmdf[duplicated(lmdf$drug_ind_index),]
    lmdf = lmdf[!duplicated(lmdf$drug_ind_index),]
    lmdf = lmdf[order(lmdf$drug_ind_index),]
    other_df = other_df[order(other_df$drug_ind_index),]
    
    evidence_top = c(evidence_top, lmdf[,ev])
    index_top = c(index_top, lmdf$drug_ind_index)        
    
    evidence_other = c(evidence_other, other_df[,ev])
    index_other = c(index_other, other_df$drug_ind_index)        
    
    size_top = c(size_top, nrow(lmdf))
    size_other = c(size_other, nrow(other_df))
  }
  #hold_out = sample(di_pairs, as.integer(.2*length(di_pairs)),replace=F)
  #held_out = di_pairs %in% hold_out
  urows = unique(di_subsel[,c("drug_ind_index","drug_indicated_for_pheno",'pair_names')])
  urows = urows[order(urows$drug_ind_index),]
  
  
  
  
  stan_input = list(N_di = length(di_pairs),
                    N_top= length(evidence_top),
                    N_other = length(evidence_other),
                    N_evidence = length(evidences),
                    evidence_to_di_top = index_top,
                    evidence_to_di_other = index_other,
                    num_per_evidences_top = size_top,
                    num_per_evidences_other = size_other,  
                    cors_top = evidence_top,
                    cors_other = evidence_other,
                    #cor_sd_wid_other = 0.2,
                    phi = 0.003852488,  ## real "phi" in test set
                    lambda = 10,
                    run_estimation = 1,
                    sd_multiplier = 2.5,
                    cor_sd_manual = 0.025)
  
  sd_ev = apply(di_subsel[,evidences],FUN=sd, MARGIN=2, na.rm=T)
  nsplit = 5
  x = matrix(nrow=nsplit, ncol=length(sd_manuals)*length(omults))
  trainres = data.frame(x)
  testres = data.frame(x)
  
  
  
  splits = rep(1:nsplit, length(di_pairs)/nsplit + 1)[1:length(di_pairs)]
  test_est = rep(0, length(di_pairs))
  names(test_est) = di_pairs
  test_label = rep(0, length(di_pairs))
  names(test_label) = di_pairs
  #cat(sd_ev,"\n")
  for(trial in 1:nsplit){
    
    held_out = di_pairs %in% di_pairs[splits==trial] #di_pairs %in% xhold_out
    held = urows[held_out,]
    cat(trial,": ",sum(held$drug_indicated_for_pheno),"\n")
    
    # while(min(tab) < 23){
    #   hold_out = sample(di_pairs, as.integer(.2*length(di_pairs)),replace=F)
    #   held_out = di_pairs %in% hold_out
    #   tab =table(held_out, lmdf$drug_indicated_for_pheno)#[held_out]==1)
    # }
    stan_input$indicated = urows$drug_indicated_for_pheno[!held_out] #
    stan_input$N_label = sum(!held_out)
    stan_input$drug_ind_index_label = urows$drug_ind_index[!held_out] #
    dfcol = c()
    col = 1
    for(omult in omults){
      stan_input$cor_sd_manual = omult*sd_ev
      #cat("for:\n",sd_ev,"\n",omult, "\n",stan_input$cor_sd_manual)
      
      #dfcol = c()
      for(sm in sd_manuals){
        stan_input$sd_multiplier = sm
        fit_samp = sampling(mod_mu, data = stan_input,
                            chains = 4, 
                            cores = 4,
                            iter = 3000,
                            warmup = 2000,
                            pars = c('p_di'), #, 'cor_sd_other'),sd_multiplier
                            #control=list(adapt_delta=.95),
                            refresh=0, seed = 123)
        
        mats <- rstan::extract(fit_samp)
        #cat("training:\n") 
        train = roc(stan_input$indicated, colMeans(mats$p_di)[stan_input$drug_ind_index_label], quiet=T)
        #cat("testing:\n") 
        p_test = colMeans(mats$p_di)[held$drug_ind_index]
        cat('split=',trial,' p-sum=',sum(p_test),' p-first=',p_test[1:10],' names-first=',rownames(held)[1:10],'\n')
        test = roc(held$drug_indicated_for_pheno, p_test, quiet=T)
        test_est[held$pair_names] = p_test
        test_label[held$pair_names] = held$drug_indicated_for_pheno
        colnm = paste(omult, sm,collapse='-')
        dfcol = c(dfcol, colnm)
        trainres[trial, col] = as.numeric(auc(train))
        testres[trial, col] = as.numeric(auc(test))
        cat("...",colnm,"\ttrain=",as.numeric(auc(train)),"\ttest=",as.numeric(auc(test)), "\n")
        col = col +1
      }
    }
  }
  colnames(trainres) = dfcol
  colnames(testres) = dfcol
  test_p_lab = data.frame(p=test_est, label=test_label,split=splits)
  return(list(di_pairs, trainres, testres, test_p_lab))
}

#### stan side effects ####
library(pROC)
mult_split_eval_SE = function(di, subsel, sd_manuals, omults,di_pairs_in=c()){
  evidences = colnames(di)[7:ncol(di)]
  di$pair_names = paste(di$drug, di$phenotype, sep=":")
  mod_mu <- stan_model('src/06.12_multi_evidence.stan',save_dso=T)
  
  di_pairs = unique(di$pair_names) 
  nsplit = 2
  if(length(di_pairs_in) > 0){
    di_pairs = di_pairs_in
    cat("using saved di_pairs",di_pairs_in[1:5],"\n")
  }else if(subsel > 0){
    set.seed(11042024) # to always get the same pairs
    di_pairs = sample(di_pairs, size=subsel)
  }
  di_subsel = di[di$pair_names %in% di_pairs,]
  N_di = length(di_pairs) 
  pair_id = 1:N_di
  names(pair_id) = di_pairs
  di_subsel$drug_ind_index = pair_id[di_subsel$pair_names]
  
  evidence_top = c()
  index_top = c()
  size_top = c()
  
  evidence_other = c()
  index_other = c()
  size_other = c()
  
  
  for(ev in evidences){
    
    lmdf = di_subsel[,c("drug_ind_index",ev)] #,'drug_indication_category','phenotype_category','body_system')]
    lmdf = lmdf[!is.na(lmdf[,ev]),]
    lmdf = lmdf[!(lmdf[,ev]==0),]
    lmdf = lmdf[order(lmdf[,ev], decreasing=T),]
    
    other_df = lmdf[duplicated(lmdf$drug_ind_index),]
    lmdf = lmdf[!duplicated(lmdf$drug_ind_index),]
    lmdf = lmdf[order(lmdf$drug_ind_index),]
    other_df = other_df[order(other_df$drug_ind_index),]
    
    evidence_top = c(evidence_top, lmdf[,ev])
    index_top = c(index_top, lmdf$drug_ind_index)        
    
    evidence_other = c(evidence_other, other_df[,ev])
    index_other = c(index_other, other_df$drug_ind_index)        
    
    size_top = c(size_top, nrow(lmdf))
    size_other = c(size_other, nrow(other_df))
  }
  #hold_out = sample(di_pairs, as.integer(.2*length(di_pairs)),replace=F)
  #held_out = di_pairs %in% hold_out
  urows = unique(di_subsel[,c("drug_ind_index","drug_indicated_for_pheno",'pair_names')])
  urows = urows[order(urows$drug_ind_index),]
  
  
  
  
  stan_input = list(N_di = length(di_pairs),
                    N_top= length(evidence_top),
                    N_other = length(evidence_other),
                    N_evidence = length(evidences),
                    evidence_to_di_top = index_top,
                    evidence_to_di_other = index_other,
                    num_per_evidences_top = size_top,
                    num_per_evidences_other = size_other,  
                    cors_top = evidence_top,
                    cors_other = evidence_other,
                    #cor_sd_wid_other = 0.2,
                    phi = 0.05400814,  ## real "phi" in test set
                    lambda = 10,
                    run_estimation = 1,
                    sd_multiplier = 2.5,
                    cor_sd_manual = 0.025)
  
  sd_ev = apply(di_subsel[,evidences],FUN=sd, MARGIN=2, na.rm=T)
  nsplit = 5
  x = matrix(nrow=nsplit, ncol=length(sd_manuals)*length(omults))
  trainres = data.frame(x)
  testres = data.frame(x)
  
  
  
  splits = rep(1:nsplit, length(di_pairs)/nsplit + 1)[1:length(di_pairs)]
  test_est = rep(0, length(di_pairs))
  names(test_est) = di_pairs
  test_label = rep(0, length(di_pairs))
  names(test_label) = di_pairs
  #cat(sd_ev,"\n")
  for(trial in 1:nsplit){
    
    held_out = di_pairs %in% di_pairs[splits==trial] #di_pairs %in% xhold_out
    held = urows[held_out,]
    cat(trial,": ",sum(held$drug_indicated_for_pheno),"\n")
    
    # while(min(tab) < 23){
    #   hold_out = sample(di_pairs, as.integer(.2*length(di_pairs)),replace=F)
    #   held_out = di_pairs %in% hold_out
    #   tab =table(held_out, lmdf$drug_indicated_for_pheno)#[held_out]==1)
    # }
    stan_input$indicated = urows$drug_indicated_for_pheno[!held_out] #
    stan_input$N_label = sum(!held_out)
    stan_input$drug_ind_index_label = urows$drug_ind_index[!held_out] #
    dfcol = c()
    col = 1
    for(omult in omults){
      stan_input$cor_sd_manual = omult*sd_ev
      #cat("for:\n",sd_ev,"\n",omult, "\n",stan_input$cor_sd_manual)
      
      #dfcol = c()
      for(sm in sd_manuals){
        stan_input$sd_multiplier = sm
        fit_samp = sampling(mod_mu, data = stan_input,
                            chains = 4, 
                            cores = 4,
                            iter = 3000,
                            warmup = 2000,
                            pars = c('p_di'), #, 'cor_sd_other'),sd_multiplier
                            #control=list(adapt_delta=.95),
                            refresh=0, seed = 123)
        
        mats <- rstan::extract(fit_samp)
        #cat("training:\n") 
        train = roc(stan_input$indicated, colMeans(mats$p_di)[stan_input$drug_ind_index_label], quiet=T)
        #cat("testing:\n") 
        p_test = colMeans(mats$p_di)[held$drug_ind_index]
        cat('split=',trial,' p-sum=',sum(p_test),' p-first=',p_test[1:10],' names-first=',rownames(held)[1:10],'\n')
        test = roc(held$drug_indicated_for_pheno, p_test, quiet=T)
        test_est[held$pair_names] = p_test
        test_label[held$pair_names] = held$drug_indicated_for_pheno
        colnm = paste(omult, sm,collapse='-')
        dfcol = c(dfcol, colnm)
        trainres[trial, col] = as.numeric(auc(train))
        testres[trial, col] = as.numeric(auc(test))
        cat("...",colnm,"\ttrain=",as.numeric(auc(train)),"\ttest=",as.numeric(auc(test)), "\n")
        col = col +1
      }
    }
  }
  colnames(trainres) = dfcol
  colnames(testres) = dfcol
  test_p_lab = data.frame(p=test_est, label=test_label,split=splits)
  return(list(di_pairs, trainres, testres, test_p_lab))
}

#### subthreshold variant analysis ####
gensup_subthreshold = function(preds, drug_targets, gwas, pheno, genome_v) {
  
  # predicted pairs
  temp_preds = preds %>% 
    filter(phenotype == pheno) %>%
    dplyr::select(drug, PredProb) %>%
    distinct()
  temp_top_drugs = unique(temp_preds$drug)
  
  if (genome_v == 19) {
    temp_top_drugs_targets = drug_targets %>%
      filter(drug %in% temp_top_drugs) %>%
      dplyr::select(drug, target, chromosome_hg19, start_hg19, end_hg19) %>%
      filter(chromosome_hg19 %in% paste("chr", 1:22, sep = "")) %>%
      distinct() %>%
      na.omit() %>%
      left_join(temp_preds, by = "drug")
    temp_genes_gr = GRanges(
      seqnames = temp_top_drugs_targets$chromosome_hg19,
      ranges   = IRanges(start = temp_top_drugs_targets$start_hg19, end = temp_top_drugs_targets$end_hg19),
      gene_id  = temp_top_drugs_targets$target,
      pred_prob = temp_top_drugs_targets$PredProb,
      drug = temp_top_drugs_targets$drug
    )
  }
  if (genome_v == 38) {
    temp_top_drugs_targets = drug_targets %>%
      filter(drug %in% temp_top_drugs) %>%
      dplyr::select(drug, target, chromosome_hg38, start_hg38, end_hg38) %>%
      filter(chromosome_hg38 %in% paste("chr", 1:22, sep = "")) %>%
      distinct() %>%
      na.omit() %>%
      left_join(temp_preds, by = "drug")
    temp_genes_gr = GRanges(
      seqnames = temp_top_drugs_targets$chromosome_hg38,
      ranges   = IRanges(start = temp_top_drugs_targets$start_hg38, end = temp_top_drugs_targets$end_hg38),
      gene_id  = temp_top_drugs_targets$target,
      pred_prob = temp_top_drugs_targets$PredProb,
      drug = temp_top_drugs_targets$drug
    )
  }
  
  gwas = gwas %>% filter(chromosome %in% paste("chr", 1:22, sep = ""))
  temp_gwas_gr = GRanges(
    seqnames = gwas$chromosome,
    ranges   = IRanges(start = gwas$base_pair_location, end = gwas$base_pair_location),
    p_value  = gwas$p_value
  )
  hits_top = findOverlaps(temp_gwas_gr, temp_genes_gr)
  gwas_in_genes = cbind(
    gwas[queryHits(hits_top), ],
    gene_target = mcols(temp_genes_gr)$gene_id[subjectHits(hits_top)],
    pred_prob = mcols(temp_genes_gr)$pred_prob[subjectHits(hits_top)],
    drug = mcols(temp_genes_gr)$drug[subjectHits(hits_top)]
  )
  gwas_in_genes = gwas_in_genes %>%
    group_by(drug) %>%
    mutate(p_value = min(p_value)) %>%
    ungroup() %>%
    dplyr::select(drug, p_value, pred_prob) %>%
    distinct() %>%
    mutate(
      pval_less_5Eneg08 = if_else(p_value < 5e-08, 1, 0),
      pval_less_5Eneg06 = if_else(p_value < 5e-06, 1, 0),
      pval_less_5Eneg04 = if_else(p_value < 5e-04, 1, 0),
      pval_less_5Eneg02 = if_else(p_value < 0.05, 1, 0),
      pval_less_1Eneg01 = if_else(p_value < 0.1, 1, 0)
    ) %>%
    mutate(phenotype = pheno, .after = drug)
  return(gwas_in_genes)
}
