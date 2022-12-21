library(Seurat)
library(BASiCS)


overdispersion_droplet <- function(data.type, force.rerun = FALSE,attribute = "selection") {
  disp.analysis.outfile <- paste0(analysis.results.dir, 'OverDisperion.analysis.', data.type, '.RData')
  disp.selc.analysis.outfile <- paste0(analysis.results.dir, 'OverDisperion.selection.analysis.', data.type, '.RData')
  disp.len.analysis.outfile <- paste0(analysis.results.dir, 'OverDisperion.length.analysis.', data.type, '.RData')
  
  samples <- get_tissue_file_names(data.type)
  
  meta.data = get_meta_data(data.type)
  
  if(attribute == "OD")
  {
    if(file.exists(disp.analysis.outfile) & (force.rerun==FALSE))
    {
      load(disp.analysis.outfile)
      return(DF_OD_res)
    }
  }
  
  
  
  
  if(attribute == "selection")
  {
    if(file.exists(disp.selc.analysis.outfile) & (force.rerun==FALSE))
    {
      load(disp.selc.analysis.outfile)
      return(DF_OD_res)
    }
  # reading the selection score data
  gene_att = read_gene_features("selection")
  gene_name = names(gene_att)
  } 
  if(attribute == "gene.len")
  {
    if(file.exists(disp.len.analysis.outfile) & (force.rerun==FALSE))
    {
      load(disp.len.analysis.outfile)
      return(DF_OD_res)
    }
    # reading the transcript length data
    gene_att = read_gene_features("gene.len")
    gene_name = names(gene_att)
  } 
  
  
  
  
  if(data.type == "TM.droplet")
  {
    young_ages = "3m" # The ages of the young mice 
    old_ages_1 = c("21m","24m") # The ages of the old mice
    old_ages_2 = c("18m","24m") # secondary old mice age in case no mice in previous "old_ages"
  }  
  if(data.type == "TM.facs")  # Set young/old ages here ! 
  {
    young_ages = "3m"
    old_ages_1 = c("18m","21m","24m")
  }
  
  DF_OD_res = c()
  for(i in 1:length(samples$organs)){
    read.file <- paste0(processed.data.dir, '/', samples$organs[i], ".", processed.files.str[data.type], ".rds")
    print(read.file)
    print(paste0("Read file ", i, " out of: ", length(samples$organs)))
    SC = readRDS(file = paste0(processed.data.dir, '/', samples$organs[i], ".", processed.files.str[data.type], ".rds")) # Current tissue seurat object
    young.ind = c(SC@meta.data$age %in% young_ages) # index for cells that came from young mouses
    old.ind = c(SC@meta.data$age %in% old_ages_1) # index for cells that came from old mouses
    old_ages = old_ages_1
    if(sum(old.ind) == 0){
      old.ind = c(SC@meta.data$age %in% old_ages_2)
      old_ages = old_ages_2
    }
    
    n_cell = SC@assays$RNA@counts@Dim[2] # Number of cells
    n_genes = SC@assays$RNA@counts@Dim[1] # Number of genes
    SC_gene_name = toupper(rownames(SC)) # Gene names (in uppercase)
    
    cell_types = SC@meta.data$cell.ontology.class # Cell types vector
    cell_types_categories = meta.data[[i]]$cell_ontology_class # Cell type names
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
    
    
    if(attribute != "OD"){
      genes_attribute = gene_att[gene_name %in% (SC_gene_name)] # Filtering the selection score to genes that are found in the current tissue
      cur_gene_name = gene_name[gene_name %in% (SC_gene_name)] # The names of the filtered genes
      
    }
       
    cells_ind = c(1:(n_cell_types+1))[-earase] # Vector of valid cell types.
    if(length(cells_ind) == 0) cells_ind = (1:(n_cell_types+1)) 
    
    for(k in cells_ind){
      
      # Differential over-dispersion test results file
      test_file = paste0(basics.dir,paste("/DVT/DVT",samples$organs[i],cell_types_categories[k],"drop 3-24 same-mean.RData"))
      if(!file.exists(test_file)){ # Checking if the file exist
        load(test_file)
      } else {
        # Old Markov chain file name
        old_file = paste0(basics.dir,paste0("/chains","chain_",paste(samples$organs[i],cell_types_categories[k],"21-24m","drop"),".Rds"))
        # Young Markov chain file name
        young_file = paste0(basics.dir,paste0("/chains","chain_",paste(samples$organs[i],cell_types_categories[k],"3m","drop"),".Rds"))
        
        if(!(file.exists(old_file) & file.exists(young_file))) next() # Checking if both old and young Markov chain files exist.
        
        # Loading Marovchain files
        chain_old = BASiCS_LoadChain(RunName = paste(samples$organs[i],cell_types_categories[k],"21-24m","drop"),StoreDir = paste0(basics.dir,"/chains"))
        chain_young = BASiCS_LoadChain(RunName = paste(samples$organs[i],cell_types_categories[k],"3m","drop"),StoreDir = paste0(basics.dir,"/chains"))
        
        
        # Differential over-dispersion test - only on genes without significant difference in mean expression
        test = BASiCS_TestDE(Chain1 = chain_old,Chain2 = chain_young,GroupLabel1 = "Old",
                             GroupLabel2 = "Young",OffSet = T,PlotOffset = F,Plot = F,
                             EpsilonM = 0, EpsilonD = log2(1.5),
                             EpsilonR = log2(1.5)/log2(exp(1)),EFDR_M = 0.10, EFDR_D = 0.10)
        
        save(test,file = paste0(basics.dir,paste("/DVT/DVT",samples$organs[i],cell_types_categories[k],"drop 3-24 same-mean.RData")))
      }
      
      # Differential over-dispersion test results table
      df = test@Results$Disp@Table
      df$GeneName = toupper(df$GeneName) # Uppercasing gene names
      df = df[!df$ResultDiffDisp %in% c("ExcludedFromTesting"),] # Filtering genes significant difference in mean expression
      
      
      if(attribute == "OD"){
        Young_sign = sum(df$ResultDiffDisp == "Young+") # Number of genes with significant higher over-dispersion for youngs
        Old_sign = sum(df$ResultDiffDisp == "Old+") # Number of genes with significant higher over-dispersion for olds
        
        # Figure 3 data
        DF_OD_res = rbind(DF_OD_res,data.frame("Organs" = samples$organs[i],"Cell_type" = cell_types_categories[k],Old_sign,Young_sign))
        
      } else {
      
        Disp = df$DispOverall # Overall over-dispersion (old and young together)
        names(Disp) = toupper(df$GeneName) 
        disp_log_fc = df$DispLog2FC # Over-dispersion log2 fold-change
        disp_fc = df$DispFC # Over-dispersion fold-change
        
        Mean_df = test@Results$Mean@Table # diferential Mean expression test result table
        Mean_df$GeneName = toupper(Mean_df$GeneName) # Uppercasing gene names
        Mean_df = Mean_df[Mean_df$GeneName %in% df$GeneName,] # taking only the genes from the over-dispersion test
        Mean_all = Mean_df$MeanOverall # overall mean expression
        Mean_fc = Mean_df$MeanFC # mean expression fold change
        Mean_log_fc = Mean_df$MeanLog2FC # mean expression log2 fold change
        Mean_old = Mean_df$Mean1 # old mean
        Mean_young = Mean_df$Mean2 # young mean
        
        
        # corr per age group
        Disp_old = df$Disp1 # old over-dispersion vector
        Disp_young = df$Disp2 # young over-dispersion vector
        
        names(Disp_old) = names(Disp)
        names(Disp_young) = names(Disp)
        
        
        # over-dispersion-selection Spearman correlation for both age groups
        att_disp_cor_old_spearman = cor(Disp_old,gene_attribute[names(Disp)],use = "complete.obs",method = "spearman")
        att_disp_cor_young_spearman = cor(Disp_young,gene_attribute[names(Disp)],use = "complete.obs",method = "spearman")
        # Spearman's correlation p values
        p_val_old = cor.test(Disp_old,gene_attribute[names(Disp)],use = "complete.obs",method = "spearman")$p.value 
        p_val_young = cor.test(Disp_young,gene_attribute[names(Disp)],use = "complete.obs",method = "spearman")$p.value
        
            
        # linear regression: over-dispersion log fold-change ~ selection + mean expression log fold-change (after scaling)
        log_fc_lm = lm(scale(disp_log_fc) ~ scale(gene_attribute[Mean_df$GeneName]) + scale(Mean_log_fc))
        
        beta_log_fc = log_fc_lm$coefficients[2] 
        p_val_fc = summary(log_fc_lm)$coef[2,4]
        
        # linear regression: log over-dispersion ~ selection + log mean expression (after scaling)
        att_disp_lm = lm(scale(log2(Disp)) ~ scale(gene_attribute[Mean_df$GeneName]) + scale(log2(Mean_all)))
        
        beta_att_lm = att_disp_lm$coefficients[2]
        p_val_att_lm =  summary(att_disp_lm)$coef[2,4]
        
        # data frame containig young and old over-dispersion-selection correlation for all cell-types
        DF_OD_res = rbind(DF_OD_res,data.frame("Organs" = samples$organs[i],"Cell_type" = cell_types_categories[k],
                                               paste0(attribute,"_disp_old_cor_spearman") = att_disp_cor_old_spearman,paste0(attribute,"_disp_young_cor_spearman") = att_disp_cor_young_spearman,
                                               p_val_old,p_val_young,
                                               paste0("beta",attribute,"all",sep = "_") = beta_att_lm,"pval_all" = p_val_att_lm,
                                               paste0("beta",attribute,"log_fc",sep = "_") = beta_log_fc,"pval_log_fc" = p_val_fc))
      
      }
      
    }
  }  # End loop on tissues 
  if(attribute == "OD"){
    print("Saving over-dispersion analysis results!")
    # Save dataframe to file: 
    save(DF_cor_res, file=disp.analysis.outfile)
  }
  if(attribute == "selection"){
    print("Saving over-dispersion - selection analysis results!")
    # Save dataframe to file: 
    save(DF_cor_res, file=disp.selc.analysis.outfile)
  }
  if(attribute == "gene.len"){
    print("Saving over-dispersion - length analysis results!")
    # Save dataframe to file: 
    save(DF_cor_res, file=disp.len.analysis.outfile)
  }
  return(DF_OD_res)
} # End overdispersion_droplet function 


overdispersion_facs <- function(data.type, force.rerun = FALSE,attribute = "selection") {
  
  disp.analysis.outfile <- paste0(analysis.results.dir, 'OverDisperion.analysis.', data.type, '.RData')
  disp.selc.analysis.outfile <- paste0(analysis.results.dir, 'OverDisperion.selection.analysis.', data.type, '.RData')
  disp.len.analysis.outfile <- paste0(analysis.results.dir, 'OverDisperion.length.analysis.', data.type, '.RData')
  
  samples <- get_tissue_file_names(data.type)
  
  meta.data = get_meta_data(data.type)
  
  if(attribute == "OD")
  {
    if(file.exists(disp.analysis.outfile) & (force.rerun==FALSE))
    {
      load(disp.analysis.outfile)
      return(DF_OD_res)
    }
  }
  
  
  
  
  if(attribute == "selection")
  {
    if(file.exists(disp.selc.analysis.outfile) & (force.rerun==FALSE))
    {
      load(disp.selc.analysis.outfile)
      return(DF_OD_res)
    }
    # reading the selection score data
    gene_att = read_gene_features("selection")
    gene_name = names(gene_att)
  } 
  if(attribute == "gene.len")
  {
    if(file.exists(disp.len.analysis.outfile) & (force.rerun==FALSE))
    {
      load(disp.len.analysis.outfile)
      return(DF_OD_res)
    }
    # reading the transcript length data
    gene_att = read_gene_features("gene.len")
    gene_name = names(gene_att)
  } 
  
  
  
  
  if(data.type == "TM.droplet")
  {
    young_ages = "3m" # The ages of the young mice 
    old_ages_1 = c("21m","24m") # The ages of the old mice
    old_ages_2 = c("18m","24m") # secondary old mice age in case no mice in previous "old_ages"
  }  
  if(data.type == "TM.facs")  # Set young/old ages here ! 
  {
    young_ages = "3m"
    old_ages_1 = c("18m","21m","24m")
  }
  
  
  
  # Over dispersion differential variability test
  DF_OD_res = c()
  for(i in 1:length(organs)){
    read.file <- paste0(processed.data.dir, '/', samples$organs[i], ".", processed.files.str[data.type], ".rds")
    print(read.file)
    print(paste0("Read file ", i, " out of: ", length(samples$organs)))
    SC = readRDS(file = paste0(processed.data.dir, '/', samples$organs[i], ".", processed.files.str[data.type], ".rds")) # Current tissue seurat object
    young.ind = c(SC@meta.data$age %in% young_ages) # index for cells that came from young mouses
    old.ind = c(SC@meta.data$age %in% old_ages_1) # index for cells that came from old mouses
    
    n_cell = SC@assays$RNA@counts@Dim[2] # Number of cells
    n_genes = SC@assays$RNA@counts@Dim[1] # Number of genes
    SC_gene_name = toupper(rownames(SC)) # Gene names (in uppercase)
    
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
    
    if(attribute != "OD"){
      genes_attribute = gene_att[gene_name %in% (SC_gene_name)] # Filtering the selection score to genes that are found in the current tissue
      cur_gene_name = gene_name[gene_name %in% (SC_gene_name)] # The names of the filtered genes
      
    }    
    
    for(k in cells_ind){
      
      # Differential over-dispersion test results file
      test_file = paste0(basics.dir,paste("/DVT/DVT",samples$organs[i],cell_types_categories[k],"same-mean.RData"))
      if(!file.exists(test_file)){ # Checking if the file exist
        load(test_file)
      } else {
        # Old Markov chain file name
        old_file = paste0(basics.dir,paste0("/chains","chain_",paste(samples$organs[i],cell_types_categories[k],"old"),".Rds"))
        # Young Markov chain file name
        young_file = paste0(basics.dir,paste0("/chains","chain_",paste(samples$organs[i],cell_types_categories[k],"young"),".Rds"))
        
        if(!(file.exists(old_file) & file.exists(young_file))) next() # Checking if both old and young Markov chain files exist.
        
        # Loading Marovchain files
        chain_old = BASiCS_LoadChain(RunName = paste(samples$organs[i],cell_types_categories[k],"old"),StoreDir = paste0(basics.dir,"/chains"))
        chain_young = BASiCS_LoadChain(RunName = paste(samples$organs[i],cell_types_categories[k],"young"),StoreDir = paste0(basics.dir,"/chains"))
        
        
        # Differential over-dispersion test - only on genes without significant difference in mean expression
        test = BASiCS_TestDE(Chain1 = chain_old,Chain2 = chain_young,GroupLabel1 = "Old",
                             GroupLabel2 = "Young",OffSet = T,PlotOffset = F,Plot = F,
                             EpsilonM = 0, EpsilonD = log2(1.5),
                             EpsilonR = log2(1.5)/log2(exp(1)),EFDR_M = 0.10, EFDR_D = 0.10)
        
        save(test,file = paste0(basics.dir,paste("/DVT/DVT",samples$organs[i],cell_types_categories[k],"same-mean.RData")))
        
      }
      
      
      # Differential over-dispersion test results table
      df = test@Results$Disp@Table
      df$GeneName = toupper(df$GeneName) # Uppercasing gene names
      df = df[!df$ResultDiffDisp %in% c("ExcludedFromTesting"),] # Filtering genes significant difference in mean expression
      
     
      if(attribute == "OD"){
        Young_sign = sum(df$ResultDiffDisp == "Young+") # Number of genes with significant higher over-dispersion for youngs
        Old_sign = sum(df$ResultDiffDisp == "Old+") # Number of genes with significant higher over-dispersion for olds
        
        # Figure 3 data
        DF_OD_res = rbind(DF_OD_res,data.frame("Organs" = samples$organs[i],"Cell_type" = cell_types_categories[k],Old_sign,Young_sign))
        
      } else {
      
        Disp = df$DispOverall # Overall over-dispersion (old and young together)
        names(Disp) = toupper(df$GeneName) 
        disp_log_fc = df$DispLog2FC # Over-dispersion log2 fold-change
        disp_fc = df$DispFC # Over-dispersion fold-change
        
        Mean_df = test@Results$Mean@Table # diferential Mean expression test result table
        Mean_df$GeneName = toupper(Mean_df$GeneName) # Uppercasing gene names
        Mean_df = Mean_df[Mean_df$GeneName %in% df$GeneName,] # taking only the genes from the over-dispersion test
        Mean_all = Mean_df$MeanOverall # overall mean expression
        Mean_fc = Mean_df$MeanFC # mean expression fold change
        Mean_log_fc = Mean_df$MeanLog2FC # mean expression log2 fold change
        Mean_old = Mean_df$Mean1 # old mean
        Mean_young = Mean_df$Mean2 # young mean
        
        
        # corr per age group
        Disp_old = df$Disp1 # old over-dispersion vector
        Disp_young = df$Disp2 # young over-dispersion vector
        
        names(Disp_old) = names(Disp)
        names(Disp_young) = names(Disp)
        
        
        # over-dispersion-selection Spearman correlation for both age groups
        att_disp_cor_old_spearman = cor(Disp_old,gene_attribute[names(Disp)],use = "complete.obs",method = "spearman")
        att_disp_cor_young_spearman = cor(Disp_young,gene_attribute[names(Disp)],use = "complete.obs",method = "spearman")
        # Spearman's correlation p values
        p_val_old = cor.test(Disp_old,gene_attribute[names(Disp)],use = "complete.obs",method = "spearman")$p.value 
        p_val_young = cor.test(Disp_young,gene_attribute[names(Disp)],use = "complete.obs",method = "spearman")$p.value
        
        
        # linear regression: over-dispersion log fold-change ~ selection + mean expression log fold-change (after scaling)
        log_fc_lm = lm(scale(disp_log_fc) ~ scale(gene_attribute[Mean_df$GeneName]) + scale(Mean_log_fc))
        
        beta_log_fc = log_fc_lm$coefficients[2] 
        p_val_fc = summary(log_fc_lm)$coef[2,4]
        
        # linear regression: log over-dispersion ~ selection + log mean expression (after scaling)
        att_disp_lm = lm(scale(log2(Disp)) ~ scale(gene_attribute[Mean_df$GeneName]) + scale(log2(Mean_all)))
        
        beta_att_lm = att_disp_lm$coefficients[2]
        p_val_att_lm =  summary(att_disp_lm)$coef[2,4]
        
        # data frame containig young and old over-dispersion-selection correlation for all cell-types
        DF_OD_res = rbind(DF_OD_res,data.frame("Organs" = samples$organs[i],"Cell_type" = cell_types_categories[k],
                                               paste0(attribute,"_disp_old_cor_spearman") = att_disp_cor_old_spearman,paste0(attribute,"_disp_young_cor_spearman") = att_disp_cor_young_spearman,
                                               p_val_old,p_val_young,
                                               paste0("beta",attribute,"all",sep = "_") = beta_att_lm,"pval_all" = p_val_att_lm,
                                               paste0("beta",attribute,"log_fc",sep = "_") = beta_log_fc,"pval_log_fc" = p_val_fc))
        
      }      
    }
  } # end loop on tissues 
  
  if(attribute == "OD"){
    print("Saving over-dispersion analysis results!")
    # Save dataframe to file: 
    save(DF_cor_res, file=disp.analysis.outfile)
  }
  if(attribute == "selection"){
    print("Saving over-dispersion - selection analysis results!")
    # Save dataframe to file: 
    save(DF_cor_res, file=disp.selc.analysis.outfile)
  }
  if(attribute == "gene.len"){
    print("Saving over-dispersion - length analysis results!")
    # Save dataframe to file: 
    save(DF_cor_res, file=disp.len.analysis.outfile)
  }
  
  return(DF_OD_res)
} # end overdispersion_facs function


