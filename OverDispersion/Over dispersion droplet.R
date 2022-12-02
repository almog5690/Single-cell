library(cowplot)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(BASiCS)
library(dplyr)

# reading the selection score data
select = read.delim("gnomad.v2.1.1.lof_metrics.by_gene.txt")
gene_name = toupper(select$gene) # making all the gene names in to upper case letters
gene_selection = select$pLI # the selection score vector
names(gene_selection) = gene_name # naming each score with the gene he belongs to

# reading the transcript length data
length_data_2 = read.delim("mart_export.txt")
glength_2 = length_data_2[,7] # gene length vector
names(glength_2) = toupper(length_data_2[,6]) # gene names
gene_name_2 = names(glength_2)


young_ages = "3m" # The ages of the young mice 
old_ages_1 = c("21m","24m") # The ages of the old mice
old_ages_2 = c("18m","24m") # secondary old mice age incase no mice in previous "old_ages"

DF_noise_drop = c()
DF_cor_selc_disp_drop = c()
DF_cor_len_disp_drop = c()
gene_OD_drop = c()
df_fc = c();df_fc_reg = c()
df_selc_reg = c();df_len_reg = c()
selc_OD_reg_data_sm_drop = c()
len_OD_reg_data_sm_drop = c()
for(i in 1:length(drop_organs)){
  SC = readRDS(file = paste(drop_organs[i],"drop","rds",sep = ".")) # current tissue seurat object
  young.ind = c(SC@meta.data$age %in% young_ages) # index for cells that came from young mouses
  old.ind = c(SC@meta.data$age %in% old_ages_1) # index for cells that came from old mouses
  old_ages = old_ages_1
  if(sum(old.ind) == 0){
    old.ind = c(SC@meta.data$age %in% old_ages_2)
    old_ages = old_ages_2
  }
  
  n_cell = SC@assays$RNA@counts@Dim[2] # Number of cells
  n_genes = SC@assays$RNA@counts@Dim[1] # Number of genes
  SC_gene_name = toupper(rownames(SC)) # Gene names (in upper case)

  cell_types = SC@meta.data$cell.ontology.class # Cell types vector
  cell_types_categories = meta.data.drop[[i]]$cell_ontology_class # Cell type names
  n_cell_types = max(cell_types) # Number of cell types
  
  earase = vector()
  
  for(ct_ind in 0:n_cell_types){
    # Filtering cell types with less than 100 cells
    # or cell types with less than 20 cells in either age groups.
    if(sum(cell_types==ct_ind)<100 | sum((cell_types==ct_ind)&(young.ind)) < 20 | sum((cell_types==ct_ind)&(old.ind)) < 20){
      earase = c(earase,ct_ind+1)
      next()
    }
  }
  
  gene_selc = gene_selection[gene_name %in% (SC_gene_name)] # Filtering the selection score to genes that are found in the current tissue
  cur_gene_name = gene_name[gene_name %in% (SC_gene_name)] # The names of the filtered genes
  
  gene_length_2 = glength_2[gene_name_2 %in% SC_gene_name] # Filtering gene length vector to genes that are found in the current tissue
  cur_gene_name_2 = gene_name_2[gene_name_2 %in% SC_gene_name] # The names of the filtered genes
  
  
  cells_ind = c(1:(n_cell_types+1))[-earase] # Vector of valid cell types.
  if(length(cells_ind) == 0) cells_ind = (1:(n_cell_types+1)) 

  for(k in cells_ind){
        
    # Differential over-dispersion test results file
    test_file = paste("/tmp/DVT/DVT",drop_organs[i],cell_types_categories[k],"drop 3-24 same-mean.RData")
    if(!file.exists(test_file)){ # Checking if the file exist
      load(test_file)
    } else {
      # Old Markov chain file name
      old_file = paste0("/tmp/chains","chain_",paste(drop_organs[i],cell_types_categories[k],"21-24m","drop"),".Rds")
      # Young Markov chain file name
      young_file = paste0("/tmp/chains","chain_",paste(drop_organs[i],cell_types_categories[k],"3m","drop"),".Rds")
      
      if(!(file.exists(old_file) & file.exists(young_file))) next() # Checking if both old and young Markov chain files exist.
      
      # Loading Marovchain files
      chain_old = BASiCS_LoadChain(RunName = paste(drop_organs[i],cell_types_categories[k],"21-24m","drop"),StoreDir = "/tmp/chains")
      chain_young = BASiCS_LoadChain(RunName = paste(drop_organs[i],cell_types_categories[k],"3m","drop"),StoreDir = "/tmp/chains")
      
      
      # Differential over-dispersion test - only on genes without significant difference in mean expression
      test = BASiCS_TestDE(Chain1 = chain_old,Chain2 = chain_young,GroupLabel1 = "Old",
                           GroupLabel2 = "Young",OffSet = T,PlotOffset = F,Plot = F,
                           EpsilonM = 0, EpsilonD = log2(1.5),
                           EpsilonR = log2(1.5)/log2(exp(1)),EFDR_M = 0.10, EFDR_D = 0.10)
      
    }
    
    # Differential over-dispersion test results table
    df = test@Results$Disp@Table
    df$GeneName = toupper(df$GeneName) # Uppercasing gene names
    df = df[!df$ResultDiffDisp %in% c("ExcludedFromTesting"),] # Filtering genes significant difference in mean expression
    
    Young_sign = sum(df$ResultDiffDisp == "Young+") # Number of genes with significant higher over-dispersion for youngs
    Old_sign = sum(df$ResultDiffDisp == "Old+") # Number of genes with significant higher over-dispersion for olds
    
    # Figure 3 data
    DF_noise_drop = rbind(DF_noise_drop,data.frame("Organs" = drop_organs[i],"Cell_type" = cell_types_categories[k],Old_sign,Young_sign))

    
    Disp = df$DispOverall # Overall over-dispersion (old and young together)
    names(Disp) = toupper(df$GeneName) 
    disp_log_fc = df$DispLog2FC # Over-dispersion log2 fold-change
    disp_fc = df$DispFC # Over-dispersion fold-change
    
    
    selc_disp_cor = cor(Disp,gene_selc[names(Disp)],use = "complete.obs") # selection and over-dispersion correlation (old and young together)
    selc_disp_log_fc_cor = cor(disp_log_fc,gene_selc[names(Disp)],use = "complete.obs",method = "spearman") # selection and over-dispersion log2 fold-change correlation
    selc_disp_fc_cor = cor(disp_fc,gene_selc[names(Disp)],use = "complete.obs",method = "spearman") # selection and over-dispersion fold-change correlation
    
    # corr per age group
    Disp_old = df$Disp1 # old over-dispersion vector
    Disp_young = df$Disp2 # young over-dispersion vector
    
    names(Disp_old) = names(Disp)
    names(Disp_young) = names(Disp)
    
    
    # over-dispersion-selection Spearman correlation for both age groups
    selc_disp_cor_old_spearman = cor(Disp_old,gene_selc[names(Disp)],use = "complete.obs",method = "spearman")
    selc_disp_cor_young_spearman = cor(Disp_young,gene_selc[names(Disp)],use = "complete.obs",method = "spearman")
    # Spearman's correlation p values
    p_val_old = cor.test(Disp_old,gene_selc[names(Disp)],use = "complete.obs",method = "spearman")$p.value 
    p_val_young = cor.test(Disp_young,gene_selc[names(Disp)],use = "complete.obs",method = "spearman")$p.value
    
    
    # data frame containig young and old over-dispersion-selection correlation for all cell-types
    DF_cor_selc_disp_drop = rbind(DF_cor_selc_disp_drop,data.frame("Organs" = drop_organs[i],"Cell_type" = cell_types_categories[k],
                                                                   "selc_disp_old_cor_spearman" = selc_disp_cor_old_spearman,"selc_disp_young_cor_spearman" = selc_disp_cor_young_spearman,
                                                                   "selc_disp_cor" = selc_disp_cor,"selc_disp_log_fc_cor" = selc_disp_log_fc_cor,"selc_disp_fc_cor" = selc_disp_fc_cor,
                                                                   p_val_old,p_val_young))
    
    
    len_disp_cor = cor(Disp,gene_length_2[names(Disp)],use = "complete.obs",method = "spearman") # Length and over-dispersion correlation (old and young together)
    
    # over-dispersion-length spearman correlation for both age groups
    len_disp_cor_old_spearman = cor(Disp_old,gene_length_2[names(Disp)],use = "complete.obs",method = "spearman")
    len_disp_cor_young_spearman = cor(Disp_young,gene_length_2[names(Disp)],use = "complete.obs",method = "spearman")
    # Spearman's correlation p values
    p_val_old = cor.test(Disp_old,gene_length_2[names(Disp)],use = "complete.obs",method = "spearman")$p.value 
    p_val_young = cor.test(Disp_young,gene_length_2[names(Disp)],use = "complete.obs",method = "spearman")$p.value
    
    len_disp_fc_cor = cor(disp_fc,gene_length_2[names(Disp)],use = "complete.obs",method = "spearman") # Length and over-dispersion fold-change correlation
    len_disp_log_fc_cor = cor(disp_log_fc,gene_length_2[names(Disp)],use = "complete.obs",method = "spearman") # Length and over-dispersion log fold-change correlation
    log_fc_pval = cor.test(disp_log_fc,gene_length_2[names(Disp)],use = "complete.obs",method = "spearman")$p.value # Length and over-dispersion log fold-change correlation p value
    
    
    DF_cor_len_disp_drop = rbind(DF_cor_len_disp_drop,data.frame("Organs" = drop_organs[i],"Cell_type" = cell_types_categories[k],len_disp_cor,
                                                                 "len_disp_old_cor_spearman" = len_disp_cor_old_spearman,"len_disp_young_cor_spearman" = len_disp_cor_young_spearman,
                                                                 p_val_old,p_val_young,len_disp_fc_cor,len_disp_log_fc_cor,log_fc_pval))
    
     
    Mean_df = test@Results$Mean@Table # diferential Mean expression test result table
    Mean_df$GeneName = toupper(Mean_df$GeneName) # Uppercasing gene names
    Mean_df = Mean_df[Mean_df$GeneName %in% df$GeneName,] # taking only the genes from the over-dispersion test
    Mean_all = Mean_df$MeanOverall # overall mean expression
    Mean_fc = Mean_df$MeanFC # mean expression fold change
    Mean_log_fc = Mean_df$MeanLog2FC # mean expression log2 fold change
    Mean_old = Mean_df$Mean1 # old mean
    Mean_young = Mean_df$Mean2 # young mean
    
    
    # linear regression: over-dispersion fold-change ~ selection + mean expression fold-change
    fc_lm = lm(disp_fc ~ gene_selc[Mean_df$GeneName] + Mean_fc) 
    
    # linear regression: over-dispersion log fold-change ~ selection + mean expression log fold-change (after scaling)
    log_fc_lm = lm(scale(disp_log_fc) ~ scale(gene_selc[Mean_df$GeneName]) + scale(Mean_log_fc))
    
    # above linear regressions data
    df_fc_reg = rbind(df_fc_reg,data.frame("Organs" = drop_organs[i],"Cell_type" = cell_types_categories[k],
                                           "beta_selc_fc" = fc_lm$coefficients[2],"pval_fc" = summary(fc_lm)$coef[2,4],
                                           "beta_selc_log_fc" = log_fc_lm$coefficients[2],"pval_log_fc" = summary(log_fc_lm)$coef[2,4]))
    # linear regression: log over-dispersion ~ selection + log mean expression (after scaling)
    selc_disp_lm = lm(scale(log2(Disp)) ~ scale(gene_selc[Mean_df$GeneName]) + scale(log2(Mean_all)))
    
    # linear regression: over-dispersion ~ selection + mean expression (For old)
    selc_disp_lm_old = lm(Disp_old ~ gene_selc[Mean_df$GeneName] + Mean_old)
    # linear regression: over-dispersion ~ selection + mean expression (For young)
    selc_disp_lm_young = lm(Disp_young ~ gene_selc[Mean_df$GeneName] + Mean_young)
    
    # above linear regression data
    df_selc_reg = rbind(df_selc_reg,data.frame("Organs" = drop_organs[i],"Cell_type" = cell_types_categories[k],
                                               "beta_selc_all" = selc_disp_lm$coefficients[2],"pval_all" = summary(selc_disp_lm)$coef[2,4],
                                               "beta_selc_old" = selc_disp_lm_old$coefficients[2],"pval_old" = summary(selc_disp_lm_old)$coef[2,4],
                                               "beta_selc_young" = selc_disp_lm_young$coefficients[2],"pval_young" = summary(selc_disp_lm_young)$coef[2,4]))
    
    # linear regression: log over-dispersion ~ length + log mean expression (after scaling)
    len_disp_lm = lm(scale(log2(Disp)) ~ scale(gene_length_2[Mean_df$GeneName]) + scale(log2(Mean_all)))
    # linear regression: over-dispersion ~ length + mean expression (For old)
    len_disp_lm_old = lm(Disp_old ~ gene_length_2[Mean_df$GeneName] + Mean_old)
    # linear regression: over-dispersion ~ length + mean expression (For young)
    len_disp_lm_young = lm(Disp_young ~ gene_length_2[Mean_df$GeneName] + Mean_young)
    # linear regression: over-dispersion log fold-change ~ length + mean expression log fold-change (after scaling)
    log_fc_lm_len = lm(scale(disp_log_fc) ~ scale(gene_length_2[Mean_df$GeneName]) + scale(Mean_log_fc))
    
    # above linear regression data
    df_len_reg = rbind(df_len_reg,data.frame("Organs" = drop_organs[i],"Cell_type" = cell_types_categories[k],
                                             "beta_len_all" = len_disp_lm$coefficients[2],"pval_all" = summary(len_disp_lm)$coef[2,4],
                                             "beta_len_old" = len_disp_lm_old$coefficients[2],"pval_old" = summary(len_disp_lm_old)$coef[2,4],
                                             "beta_len_young" = len_disp_lm_young$coefficients[2],"pval_young" = summary(len_disp_lm_young)$coef[2,4],
                                             "beta_len_fc" = log_fc_lm_len$coefficients[2],"pval_fc" = summary(log_fc_lm_len)$coef[2,4]))
    
    
    
    
    # Linear regression: selection ~ old over-dispersion + old mean  + length
    selc_OD_old_reg = lm(gene_selc[names(Disp)] ~ Disp_old + Mean_old + gene_length_2[names(Disp)])
    # Linear regression: selection ~ young over-dispersion + young mean  + length
    selc_OD_young_reg = lm(gene_selc[names(Disp)] ~ Disp_young + Mean_young + gene_length_2[names(Disp)])
    
    p_val_old = summary(selc_OD_old_reg)$coef[2,4]
    p_val_young = summary(selc_OD_young_reg)$coef[2,4]
    
    # selection and over-dispersion regression data
    selc_OD_reg_data_sm_drop = rbind(selc_OD_reg_data_sm_drop,data.frame("Organs" = drop_organs[i],"Cell_type" = cell_types_categories[k],
                                                                         "selc_OD_old_reg" = selc_OD_old_reg$coef[2],"selc_OD_young_reg" = selc_OD_young_reg$coef[2],
                                                                         p_val_old,p_val_young))
    
    
    # Linear regression: Length ~ old over-dispersion + old mean  + selection
    len_OD_old_reg = lm(gene_length_2[names(Disp)] ~ Disp_old + Mean_old + gene_selc[names(Disp)])
    # Linear regression: Length ~ young over-dispersion + young mean  + selection
    len_OD_young_reg = lm(gene_length_2[names(Disp)] ~ Disp_young + Mean_young + gene_selc[names(Disp)])
    
    p_val_old = summary(len_OD_old_reg)$coef[2,4]
    p_val_young = summary(len_OD_young_reg)$coef[2,4]
    
    
    # Length and over-dispersion regression data
    len_OD_reg_data_sm_drop = rbind(len_OD_reg_data_sm_drop,data.frame("Organs" = drop_organs[i],"Cell_type" = cell_types_categories[k],
                                                                       "len_OD_old_reg" = len_OD_old_reg$coef[2],"len_OD_young_reg" = len_OD_young_reg$coef[2],
                                                                       "len_log_OD_old_reg" = len_log_OD_old_reg$coef[2],"len_log_OD_young_reg" = len_log_OD_young_reg$coef[2],p_val_old,p_val_young,p_val_log_old,p_val_log_young))
    
    
  }
}
