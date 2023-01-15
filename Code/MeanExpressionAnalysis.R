library(Seurat)


# Perform mean expression analysis for single cell RNA seq dataset
mean_expression_analysis <- function(data.type, expression.stat.y = c("overdispersion"), # dependent variable
                                     expression.stat.x = c("mean"), feature.types = c("selection"), # covariates
                                     force.rerun = FALSE) {
  # Set data directories: 
  set_data_dirs(data.type)
  mean.analysis.outfile <- paste0(analysis.results.dir, 'mean.analysis.', data.type, '.', 
                                  paste0( feature.types, collapse="_"), '.RData')
  samples <- get_tissue_file_names(data.type)
  meta.data = get_meta_data(data.type)
  if(file.exists(mean.analysis.outfile) & (force.rerun==FALSE))
  {
    load(mean.analysis.outfile)
    return(DF_cor)
  }
  
  # reading the selection score data
  print("Read gene features")
  n.features <- length(feature.types)
  gene_features = read_gene_features(feature.types)  # Problem: we need the same number of genes !!! 
  gene_name <- vector("list", n.features)
  names(gene_name) <- feature.types
  
  #  gene_name = names(gene_features)
  groups = dataset_to_age_groups(data.type)
  
  # Set column names: 
  col.names <- c("Organs", "Cell_type")
  for(feature.type in feature.types)
  {
    gene_name[[feature.type]] <- names(gene_features[[feature.type]])
    # Marginal coefficients (Spearman)    
    for(age.group in c("all", "young", "old"))
    {
      col.names <- c(col.names, c(paste0(feature.type, "_", age.group, "_cor"), 
                                  paste0(feature.type, "_", age.group, "_pval")))
    }
    col.names <- c(col.names, c(paste0(feature.type, "_fc_cor"), paste0(feature.type, "_fc_pval"), 
                                paste0(feature.type, "_fc_abs_cor"), paste0(feature.type, "_fc_abs_pval")))
  }    
  beta.inds <- seq(length(col.names)+1, length(col.names)+n.features*8, 2)
  beta.pvals.inds <- seq(length(col.names)+2, length(col.names)+n.features*8, 2)
  for(feature.type in feature.types)     # Multiple regression coefficients     
  {
    for(age.group in c("all", "young", "old"))
    {
      col.names <- c(col.names, c(paste0(feature.type, "_", age.group, "_beta"), 
                                  paste0(feature.type, "_", age.group, "_beta_pval")))
    }
    col.names <- c(col.names, c(paste0(feature.type, "_fc_beta"), paste0(feature.type, "_fc_beta_pval"), 
                                paste0(feature.type, "_fc_abs_beta"), paste0(feature.type, "_fc_abs_beta_pval")))
  }
  DF_cor <- data.frame(matrix(ncol = length(col.names), nrow = 0, dimnames=list(NULL, col.names)))
  cell.type.ctr <- 1
  for(i in 1:length(samples$organs)){  #
    # New: use utility to extract statistics: 
    print(paste0("Read expression stats ", i, " out of ", length(samples$organs), ": ", basename(read.file)))
    expression.features.by.cell.type.and.age <- extract_expression_statistics(data.type, samples$organs[i], 
                                                                              expression.stats = c("mean"), SeuratOutput=c()) # extract means

    read.file <- paste0(processed.data.dir, samples$organs[i], ".", processed.files.str[data.type], ".rds")
    print(paste0("Read file ", i, " out of ", length(samples$organs), ": ", basename(read.file)))
    SC = readRDS(file = read.file) # Current tissue seurat object
    counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the current tissue
    list2env(tissue_to_age_inds(data.type, samples$organs[i], groups, SC@meta.data), env=environment()) # set specific ages for all age groups in all datasets
    
    #    n_cell = SC@assays$RNA@counts@Dim[2] # Number of cells
    #    n_genes = SC@assays$RNA@counts@Dim[1] # Number of genes
    SC_gene_name = toupper(rownames(SC)) # Gene names (in uppercase)
    if(length(SC_gene_name) != dim(counts.mat)[1]) { # For Rats, names are already in the matrix
      SC_gene_name = toupper(rownames(counts.mat)) 
    } 
    rownames(counts.mat) = SC_gene_name # make sure upper 
    
    if(data.type == "CR.Rat"){ # Convert to dummy variables 
      cell_types = as.numeric(SC@meta.data$cell_types) # Cell types vector
      cell_types_categories = levels(SC$cell_types)
    } else
    {
      cell_types = SC@meta.data$cell.ontology.class # Cell types vector
      cell_types_categories = meta.data[[i]]$cell_ontology_class # Cell type names. Missing variable meta.data.drop
    }
    
    cells_ind = filter_cells(cell_types, young.ind, old.ind, filter.params)
    
    cur_gene_features <- vector("list", n.features)
    cur_gene_name <- vector("list", n.features)  # genes that both have the feature and are present in the sc-RNA-seq data of the tissue
    names(cur_gene_features) <- feature.types
    names(cur_gene_name) <- feature.types
    
    for(feature.type in feature.types)
    {
      cur_gene_features[[feature.type]] = gene_features[[feature.type]][gene_name[[feature.type]] %in% (SC_gene_name)] # filtering the selection score to genes that are found in the current tissue
      cur_gene_name[[feature.type]] = gene_name[[feature.type]][gene_name[[feature.type]] %in% (SC_gene_name)] # the names of the filtered genes in this tissue
      #      cur_counts.mat <- counts.mat[cur_gene_name,] # Take only relevant genes
    }
    print("Cells inds:")
    print(cells_ind)
    for(k in cells_ind){ # loop on cell types . 
      print(paste0("Analyze cell type: ", cell_types_categories[k]))
      for(feature.type in feature.types)
      {
        # 1. Pearson correlations for all, young, old
#        gene.mean.by.age.group <- data.frame(matrix(ncol = 3, nrow = length(cur_gene_name[[feature.type]])))
#        colnames(gene.mean.by.age.group) <- c("all", "young", "old")

        for(age.group in c("all", "young", "old"))  # colnames(gene.mean.by.age.group))
        {
#          cur.ind = switch(age.group, # Indices of cells in each group
#                           "all" = all.ind, 
#                           "young" = young.ind, 
#                           "old" = old.ind)  & (cell_types==k-1) # take only cell type 
#          cur.gene.ind = rowSums(SC@assays$RNA@counts[,cur.ind]) > filter.params$min.count # Filter genes with low expression for old/young (less than 10 counts). Keep indices of genes
#          names(cur.gene.ind) = toupper(names(cur.gene.ind))
#          cur.gene.ind = cur.gene.ind[cur_gene_name[[feature.type]]]
#          gene.mean.by.age.group[age.group] <- rowMeans(counts.mat[cur_gene_name[[feature.type]], cur.ind])  # Take only filtered cells
#          DF_cor[cell.type.ctr, c(paste0(feature.type, "_", age.group, "_pval"), 
#                                  paste0(feature.type, "_", age.group, "_cor"))] <- 
#            cor.test(gene.mean.by.age.group[[age.group]][cur.gene.ind], 
#                     cur_gene_features[[feature.type]][cur.gene.ind], 
#                     use = "complete.obs", method = "spearman")[3:4]
          
          # New: use the extracted means 
#          print(cur.gene.ind)
#          print("Length, sum:")
#          print(length(cur.gene.ind))
#          print(sum(cur.gene.ind))
#          print(length(names(cur.gene.ind)))
          
#          print("Dimensions:")
#          print(   length(gene.mean.by.age.group[[age.group]][cur.gene.ind])  )
#          print(  length( expression.features.by.cell.type.and.age[[k]][cur.gene.ind, paste0("mean_", age.group)]  ) )
#          print(  length(cur_gene_features[[feature.type]][cur.gene.ind]) )
          print("Feature type: ")
          print(feature.type)
          print(age.group)
          print(k)
          print(length(expression.features.by.cell.type.and.age))
          cur.gene.ind <- expression.features.by.cell.type.and.age[[k]][, paste0("mean_", age.group)] > -1  # Set current gene inds 
          
          print("cur.gene.ind: ")
          print(sum(cur.gene.ind))
          DF_cor[cell.type.ctr, c(paste0(feature.type, "_", age.group, "_pval"), 
                                  paste0(feature.type, "_", age.group, "_cor"))] <- 
            cor.test(expression.features.by.cell.type.and.age[[k]][which(cur.gene.ind), paste0("mean_", age.group)], 
                     cur_gene_features[[feature.type]][cur.gene.ind], 
                     use = "complete.obs", method = "spearman")[3:4]
          
          print("Computed cor ")
          
        }
        
        print("Finished Age groups")
        # 2. Fold-change       
        # Filtering genes with low expression for old AND for young (less than 10 counts)
#        fc.gene.ind = rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&young.ind]) > filter.params$min.count |
#          rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&old.ind]) > filter.params$min.count  # choose genes for fold-change. Filter on both young and old
#        names(fc.gene.ind) = toupper(names(fc.gene.ind))
#        fc.gene.ind = fc.gene.ind[cur_gene_name[[feature.type]]]
#        mean_young = gene.mean.by.age.group[["young"]][fc.gene.ind]  # filtered young mean expression vector
#        mean_old = gene.mean.by.age.group[["old"]][fc.gene.ind]  # filtered young mean expression vector
        
        
        
        # Mean expression vs. fold-change and selection correlation
        fc.gene.ind <- (expression.features.by.cell.type.and.age[[k]][, "mean_young"] > -1) |
          (expression.features.by.cell.type.and.age[[k]][, "mean_old"] > -1)# Set current gene inds 
        
        print("Compute cor fold-change")
        DF_cor[cell.type.ctr, c(paste0(feature.type, "_fc_pval"), paste0(feature.type, "_fc_cor"))] <- 
          cor.test(log(expression.features.by.cell.type.and.age[[k]][fc.gene.ind, "mean_old"] / 
                         expression.features.by.cell.type.and.age[[k]][fc.gene.ind, "mean_young"]), 
                   cur_gene_features[[feature.type]][fc.gene.ind],  use = "complete.obs", method = "spearman")[3:4]  # Take log of fold-change. Maybe take difference? (they're after log)

        # Mean expression vs. Absolute fold-change and selection correlation
#        DF_cor[cell.type.ctr, c(paste0(feature.type, "_fc_abs_pval"), paste0(feature.type, "_fc_abs_cor"))] <- 
#          cor.test(abs(log(mean_old/mean_young)), cur_gene_features[[feature.type]][fc.gene.ind],  use = "complete.obs", method = "spearman")[3:4]  # Take log of fold-change. Maybe take difference? (they're after log)

        print("Compute cor absolute fold-change")
        DF_cor[cell.type.ctr, c(paste0(feature.type, "_fc_abs_pval"), paste0(feature.type, "_fc_abs_cor"))] <- 
          cor.test(abs(log(expression.features.by.cell.type.and.age[[k]][fc.gene.ind, "mean_old"] / 
                         expression.features.by.cell.type.and.age[[k]][fc.gene.ind, "mean_young"])), 
                   cur_gene_features[[feature.type]][fc.gene.ind],  use = "complete.obs", method = "spearman")[3:4]  # Take log of fold-change. Maybe take difference? (they're after log)
        
        # Finally, set other fields      
        print("Set other fields")
        DF_cor[cell.type.ctr, "Organs"] <- samples$organs[i]
        DF_cor[cell.type.ctr, "Cell_type"] <- cell_types_categories[k]

      } # end loop on feature type
      
      # Add multiple linear regression with all features together! (for each age group separately!)
      reg.ctr = 1
      all.features.gene.names <- SC_gene_name  # genes that both appear in the tissue, and have ALL features 
      for(feature.type in feature.types)  # get few in intersection
        all.features.gene.names <- intersect(all.features.gene.names, names(gene_features[[feature.type]]))
      n.gene.with.all.features <- length(all.features.gene.names)
      cur_gene_features_mat <- matrix(data = 0, nrow = n.gene.with.all.features, ncol = n.features)         # Intersect with data gene names 
      colnames(cur_gene_features_mat) <- feature.types
      for(feature.type in feature.types)
        cur_gene_features_mat[, feature.type] <- gene_features[[feature.type]][all.features.gene.names]
      cur_gene_features_mat <- as.data.frame(cur_gene_features_mat)
      row.names(cur_gene_features_mat) <- all.features.gene.names
      
      print("Do linear regression:") # Fix also regression: 
      for(age.group in c("all", "young", "old")) # colnames(gene.mean.by.age.group))
      {
#        cur.ind = switch(age.group, # Indices of cells in each group
#                         "all" = all.ind, 
#                         "young" = young.ind, 
#                         "old" = old.ind)  & (cell_types==k-1) # take only cell type 
#        cur.gene.ind = rowSums(SC@assays$RNA@counts[,cur.ind]) > filter.params$min.count # Filter genes with low expression for old/young (less than 10 counts). Keep indices of genes
#        names(cur.gene.ind) = toupper(names(cur.gene.ind))
#        no.features.gene.names = setdiff(names(cur.gene.ind), all.features.gene.names)        
#        cur.gene.ind[no.features.gene.names] = FALSE # take only genes with features 
#        
#        gene.mean.by.age.group.reg <- rowMeans(counts.mat[cur.gene.ind, cur.ind])  # Take only filtered cells
        
        # New: get extracted expression features
        cur.gene.ind <- expression.features.by.cell.type.and.age[[k]][, paste0("mean_", age.group)] > -1  # Set current gene inds 
        reg.model <- lm(expression.features.by.cell.type.and.age[[k]][which(cur.gene.ind), paste0("mean_", age.group)] ~ ., 
                        data = as.data.frame(cur_gene_features_mat[cur.gene.ind,]))  # Take log of fold-change. Maybe take difference? (they're after log)
        DF_cor[cell.type.ctr, beta.inds[((reg.ctr-1)*n.features+1):(reg.ctr*n.features)]] <- reg.model$coefficients[-1] # get p-values (excluding intercept)
        DF_cor[cell.type.ctr, beta.pvals.inds[((reg.ctr-1)*n.features+1):(reg.ctr*n.features)]] <- summary(reg.model)$coefficients[-1,4]  # get p-values (excluding intercept?)
        reg.ctr = reg.ctr + 1
      }
      
      # Do regression for fold-change 
      fc.gene.ind = rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&young.ind]) > filter.params$min.count &  # require both !! 
        rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&old.ind]) > filter.params$min.count  # choose genes for fold-change. Filter on both young and old
      names(fc.gene.ind) = toupper(names(fc.gene.ind))
      no.features.gene.names = setdiff(names(fc.gene.ind), all.features.gene.names)        
      fc.gene.ind[no.features.gene.names] = FALSE # take only genes with features 
      mean_young = rowMeans(counts.mat[fc.gene.ind, young.ind])   # filtered young mean expression vector
      mean_old = rowMeans(counts.mat[fc.gene.ind, old.ind])   # filtered young mean expression vector
      
      fc.reg.model <- lm(log(mean_old / mean_young) ~ ., data = as.data.frame(cur_gene_features_mat[fc.gene.ind,]))  # Take log of fold-change. Maybe take difference? (they're after log)
      DF_cor[cell.type.ctr, beta.inds[((reg.ctr-1)*n.features+1):(reg.ctr*n.features)]] <- fc.reg.model$coefficients[-1] # get p-values (excluding intercept)
      DF_cor[cell.type.ctr, beta.pvals.inds[((reg.ctr-1)*n.features+1):(reg.ctr*n.features)]] <- summary(fc.reg.model)$coefficients[-1,4]  # get p-values (excluding intercept?)
      reg.ctr = reg.ctr + 1
      
      cell.type.ctr = cell.type.ctr+1  # update counter
    }  # end loop on cell-types in tissue
  }  # end loop on tissues
  
  save(DF_cor, file=mean.analysis.outfile)   # Save dataframe to file ( Save also to excel? ) 
  return(DF_cor)
} # End function mean expression droplet 
