library(Seurat)



# Perform mean expression analysis for TM droplet dataset 
mean_expression_analysis <- function(data.type, feature.type = "selection", covariates = NULL, force.rerun = FALSE) {
  mean.analysis.outfile <- paste0(analysis.results.dir, 'mean.analysis.', data.type, '.RData')
  samples <- get_tissue_file_names(data.type)
  
  meta.data = get_meta_data(data.type)
  
  # Set parameters for filtering (may depend on the data type?): 
  min.count <- 10 
  min.cells.total <- 100
  min.cells.per.age <- 20
  
  if(file.exists(mean.analysis.outfile) & (force.rerun==FALSE))
  {
    load(mean.analysis.outfile)
    return(DF_cor)
  }
  
  # reading the selection score data
  gene_features = read_gene_features(feature.type)
  gene_name = names(gene_features)
  
  
  if(data.type == "TM.droplet")
  {
    young_ages = "3m" # The ages of the young mice 
    old_ages_1 = c("21m","24m") # The ages of the old mice
    old_ages_2 = c("18m","24m") # secondary old mice age in case no mice in previous "old_ages"
  }  
  if(data.type == "TM.facs")  # Set young/old ages here ! 
  {
    young_ages = "3m" # The ages of the young mice 
    old_ages_1 = c("18m", "21m", "24m")
    old_ages_2 = c()
  }
  if(data.type == "CR.Rat")  # Set young/old ages here ! 
  {
    young_ages = "Y" # The ages of the young mice 
    old_ages_1 = "O"
    old_ages_2 = c()
  }
  DF_cor = c()
  for(i in 1:length(samples$organs)){
    read.file <- paste0(processed.data.dir, '/', samples$organs[i], ".", processed.files.str[data.type], ".rds")
    print(read.file)
    print(paste0("Read file ", i, " out of: ", length(samples$organs)))
    SC = readRDS(file = paste0(processed.data.dir, '/', samples$organs[i], ".", processed.files.str[data.type], ".rds")) # Current tissue seurat object
    counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the current tissue
    young.ind = c(SC@meta.data$age %in% young_ages) # index for cells that came from 3 month old mouses
    old.ind = c(SC@meta.data$age %in% old_ages_1) # index for cells that came from old mouses
    if(sum(old.ind) == 0){ # Empty
      old.ind = c(SC@meta.data$age %in% old_ages_2)
    }
    if((data.type == "TM.droplet") & (i==6)){ # special tissue for droplet (which?). Should move this line
      old.ind = SC@meta.data$age %in% c("18m","21m")
    }
    
    n_cell = SC@assays$RNA@counts@Dim[2] # Number of cells
    n_genes = SC@assays$RNA@counts@Dim[1] # Number of genes
    SC_gene_name = toupper(rownames(SC)) # Gene names (in uppercase)
    rownames(counts.mat) = SC_gene_name 
    all.ind = rep(TRUE, n_cell)
    
    cell_types = SC@meta.data$cell.ontology.class # Cell types vector
    cell_types_categories = meta.data[[i]]$cell_ontology_class # Cell type names. Missing variable meta.data.drop
    if(data.type == "CR.Rat"){
      cell_types_categories = levels(SC$cell_types)
    }
    n_cell_types = max(cell_types) # Number of cell types
    earase = vector()
    
    for(ct_ind in 0:n_cell_types){
      #filtering cell types with less then 100 cells or less the 20 cells in each the age groups
      if(sum(cell_types==ct_ind)<min.cells.total|sum(cell_types==ct_ind & young.ind)<min.cells.per.age|sum(cell_types==ct_ind & old.ind)<min.cells.per.age){
        earase = c(earase,ct_ind+1)
        next()
      }
    }
    
    cur_gene_features = gene_features[gene_name %in% (SC_gene_name)] # filtering the selection score to genes that are found in the current tissue
    cur_gene_name = gene_name[gene_name %in% (SC_gene_name)] # the names of the filtered genes
    cur_counts.mat <- counts.mat[cur_gene_name,] # Take only relevant genes
    
    cells_ind = c(1:(n_cell_types+1))[-earase]
    if(length(cells_ind) == 0) cells_ind = (1:(n_cell_types+1))
    
    print("Start loop cell types")
    for(k in cells_ind){ # loop on cell types 
      #      cur_gene_features = gene_features[gene_name %in% (SC_gene_name)] # filtering the selection score to genes that are found in the current tissue
      
      
      # 1. Pearson correlations for all, young, old
      gene_mean = rowMeans(cur_counts.mat[,cell_types==k-1])
      for(age.group in c("all", "young", "old"))
      {
        cur.ind = switch(age.group, 
                         "all" = all.ind, 
                         "young" = young.ind, 
                         "old" = old.ind)
        
        # Filtering genes with low expression for old/young (less than 10 counts)
        cur.gene.ind = rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&(!cur.ind)]) > min.count # Filter 
        names(cur.gene.ind) = toupper(names(cur.gene.ind))
        cur.gene.ind = cur.gene.ind[cur_gene_name]
        cur.gene.mean <- gene_mean[cur.ind]
        cur.gene.mean <- cur.gene.mean[cur.gene.ind]
        # genes mean vectors: all, young, old
        list2env(setNames(cor.test(gene_mean[cur.ind], cur_gene_features[cur.ind], use = "complete.obs", method = "spearman")[3:4], 
                          c(paste0("mean_", feature.type, "_", age.group, "_cor"), paste0("mean_", feature.type, "_", age.group, "_pval"))), envir = .GlobalEnv)
      }
      
      # genes mean vectors: all, young, old
      #      gene_mean = rowMeans(cur_counts.mat[,cell_types==k-1])
      #      gene_mean_young = rowMeans(cur_counts.mat[,(cell_types==k-1)&(young.ind)])
      #      gene_mean_old = rowMeans(cur_counts.mat[,(cell_types==k-1)&(old.ind)])
      
      # Mean expression and selection correlation (old and young together)
      #      list2env(setNames(cor.test(gene_mean, cur_gene_features, use = "complete.obs", method = "spearman")[3:4], 
      #               c("pval_all", "mean_selc_cor")), envir = .GlobalEnv)
      
      # 2. Fold-change       
      # Filtering genes with low expression for old (less than 10 counts)
      genes_ind = rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&(!young.ind)]) > min.count
      names(genes_ind) = toupper(names(genes_ind))
      genes_ind_old = genes_ind[cur_gene_name]
      
      # Filtering genes with low expression for young (less than 10 counts)
      genes_ind = rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&(young.ind)]) > min.count
      names(genes_ind) = toupper(names(genes_ind))
      genes_ind_young = genes_ind[cur_gene_name]
      
      genes_fc_ind = rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&(!young.ind)]) > min.count & rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&(!young.ind)]) > min.count
      
      mean_young = gene_mean_young[genes_fc_ind] # filtered young mean expression vector
      mean_old = gene_mean_old[genes_fc_ind] # filtered young mean expression vector
      
      # Mean expression fold-change and selection correlation
      #      list2env(setNames(cor.test((mean_old/mean_young), cur_gene_features[genes_ind_young], 
      #                                 use = "complete.obs", method = "spearman")[3:4], 
      list2env(setNames(cor.test((mean_old/mean_young), cur_gene_features[genes_fc_ind],
                                 use = "complete.obs",method = "spearman")[3:4], 
                        c("pval_fc", "fc_cor")), envir = .GlobalEnv)
      
      
      gene_selc_old = cur_gene_features[genes_ind_old]
      gene_mean_old = gene_mean_old[genes_ind_old]
      gene_selc_young = cur_gene_features[genes_ind_young]
      gene_mean_young = gene_mean_young[genes_ind_young]
      
      # mean-selection Spearman correlation for both age groups
      old_cor_spearman = cor(gene_selc_old,gene_mean_old, use = "complete.obs",method = "spearman")
      young_cor_spearman = cor(gene_selc_young,gene_mean_young, use = "complete.obs",method = "spearman")
      # mean-selection Spearman correlation P values for both age groups
      p_val_old = cor.test(gene_selc_old, gene_mean_old, use = "complete.obs", method = "spearman")$p.value
      p_val_young = cor.test(gene_selc_young, gene_mean_young, use = "complete.obs", method = "spearman")$p.value
      
      
      # data frame containing young and old mean-selection correlation for all cell-types
      DF_cor = rbind(DF_cor, data.frame("Organs" = samples$organs[i],"Cell_type" = cell_types_categories[k],mean_selc_cor,pval_all,
                                        "selc_mean_old_cor_spearman" = old_cor_spearman,"selc_mean_young_spearman" = young_cor_spearman,
                                        p_val_old,p_val_young,fc_cor,pval_fc))
      
    }  # end loop on cell-types in tissue
    print("End loop cell types")
    
  }  # end loop on tissues
  
  print("Saving DF core!")
  # Save dataframe to file: 
  save(DF_cor, file=mean.analysis.outfile)
  # Save also to excel? 
  return(DF_cor)
  
} # End function mean expression droplet 



