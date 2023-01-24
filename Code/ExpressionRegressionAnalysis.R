library(Seurat)


# Perform mean expression analysis for single cell RNA seq dataset
# Input: 
# data.type - dataset name
# expression.stat.y - dependent expression variable in the regression
# expression.stat.x - independent expression variable in the regression
# feature.types - gene features variables in the regression
# force.rerun - compute again, do not read from file 
# Output: 
# DF_cor - data-frame with results, also saved to file 
expression_regression_analysis <- function(data.type, expression.stat.y = c("mean"), # dependent variable
                                     expression.stat.x = c(), feature.types = c("selection"), # covariates
                                     force.rerun = FALSE) {
  # Set data directories: 
  set_data_dirs(data.type)
  reg.analysis.outfile <- paste0(analysis.results.dir, 'expr.reg.analysis.', data.type, '.', expression.stat.y, '.vs.', 
                                  paste0( expression.stat.x, collapse="_"), rep('.', min(length(expression.stat.x), 1)),
                                  paste0( feature.types, collapse="_"), '.RData')
  samples <- get_tissue_file_names(data.type)
  meta.data = get_meta_data(data.type)
  if(file.exists(reg.analysis.outfile) & (force.rerun==FALSE))
  {
    load(reg.analysis.outfile)
    return(DF_cor)
  }
  
  # reading the selection score data
  print("Read gene features")
  n.expr.features <- length(expression.stat.x)  # additional expression features (may allow later to take squares, combinations ...)
  n.features <- length(feature.types)
  gene_features = read_gene_features(feature.types)  # Problem: we need the same number of genes !!! 
  gene_name <- vector("list", n.features)
  names(gene_name) <- feature.types
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
  DF_cor <- data.frame(matrix(ncol = length(col.names), nrow = 0, dimnames=list(NULL, col.names))) # start with empty data-frame. Add cell types that work
  cell.type.ctr <- 1
  expression.stats <- c(expression.stat.y,  expression.stat.x) # load expression statistics
  for(i in 1:length(samples$organs)){  # SKIN, SCAT, .. Many files in facs couldn't be read. Magin number error (?!) 
    # New: use utility to extract statistics: can be very heavy 
    expr.stats <- extract_expression_statistics(data.type, samples$organs[i], 
                                                expression.stats = expression.stats, SeuratOutput=c(), force.rerun = FALSE) # extract means
    
    
    read.file <- paste0(processed.data.dir, samples$organs[i], ".", processed.files.str[data.type], ".rds")
    print(paste0("Read file ", i, " out of ", length(samples$organs), ": ", basename(read.file)))
    cells_ind = which(!unlist(lapply(expr.stats$DF.expr.stats, is.null)))
    SC_gene_name = toupper( rownames(expr.stats$DF.expr.stats[[cells_ind[1]]]) )

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

    for(k in cells_ind){ # loop on cell types . 
      print(paste0("Analyze cell type: ", expr.stats$cell_types_categories[k]))
      for(feature.type in feature.types)
      {
        # 1. Pearson correlations for all, young, old
        for(age.group in c("all", "young", "old"))  # colnames(gene.mean.by.age.group))
        {
          cur.gene.ind <- expr.stats$DF.expr.stats[[k]][, paste0(expression.stat.y, "_", age.group)] > -1  # Set current gene inds 
#          print(paste0("cur.gene.ind: ", length(cur.gene.ind), ", ", sum(cur.gene.ind)))
#          print("unique:")
#          print(unique(cur.gene.ind))
          if(length(which(cur.gene.ind)) == 0)  # no genes, probably bad data
          {
#            print(paste0("Skipping age group ", age.group, ", no data!"))
            next
          }
#          print(cur.gene.ind[1:10])
#          print("NA cells:")
#          print(which(is.na(cur.gene.ind)))
#          
#          print(paste0("Cor x:", length( expr.stats$DF.expr.stats[[k]][which(cur.gene.ind), paste0(expression.stat.y, "_", age.group)]  )))
#          print(paste0("Cor y:", length( cur_gene_features[[feature.type]][cur.gene.ind]  )))
          DF_cor[cell.type.ctr, c(paste0(feature.type, "_", age.group, "_pval"), 
                                  paste0(feature.type, "_", age.group, "_cor"))] <- 
            cor.test(expr.stats$DF.expr.stats[[k]][which(cur.gene.ind), paste0(expression.stat.y, "_", age.group)], 
                     cur_gene_features[[feature.type]][which(cur.gene.ind)], 
                     use = "complete.obs", method = "spearman")[3:4]
        }
        
        # 2. Fold-change       
        # Filtering genes with low expression for old AND for young (less than 10 counts)
        # Mean expression vs. fold-change and selection correlation
        fc.gene.ind <- (expr.stats$DF.expr.stats[[k]][, paste0(expression.stat.y, "_young")] > -1) |
          (expr.stats$DF.expr.stats[[k]][,  paste0(expression.stat.y, "_old")] > -1)# Set current gene inds 
        
        if(length(which(fc.gene.ind)) == 0)  # no genes, probably bad data
          next
        DF_cor[cell.type.ctr, c(paste0(feature.type, "_fc_pval"), paste0(feature.type, "_fc_cor"))] <- 
          cor.test(log(expr.stats$DF.expr.stats[[k]][fc.gene.ind,  paste0(expression.stat.y, "_old")] / 
                         expr.stats$DF.expr.stats[[k]][fc.gene.ind,  paste0(expression.stat.y, "_young")]), 
                   cur_gene_features[[feature.type]][fc.gene.ind],  use = "complete.obs", method = "spearman")[3:4]  # Take log of fold-change. Maybe take difference? (they're after log)

        # Mean expression vs. Absolute fold-change and selection correlation
        DF_cor[cell.type.ctr, c(paste0(feature.type, "_fc_abs_pval"), paste0(feature.type, "_fc_abs_cor"))] <- 
          cor.test(abs(log(expr.stats$DF.expr.stats[[k]][fc.gene.ind, paste0(expression.stat.y, "_old")] / 
                         expr.stats$DF.expr.stats[[k]][fc.gene.ind, paste0(expression.stat.y, "_young")])), 
                   cur_gene_features[[feature.type]][fc.gene.ind],  use = "complete.obs", method = "spearman")[3:4]  # Take log of fold-change. Maybe take difference? (they're after log)
        
        # Finally, set other fields      
        DF_cor[cell.type.ctr, "Organs"] <- samples$organs[i]
        DF_cor[cell.type.ctr, "Cell_type"] <- expr.stats$cell_types_categories[k]
      } # end loop on feature type
      
      if(length(which(fc.gene.ind)) == 0)  # no genes, probably bad data
      {
        print(paste0("Skipping cell, no data!"))
        next
      }
      
      
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
      
      # New: Add possible expression covariates to gene features: 
      
      for(age.group in c("all", "young", "old")) # colnames(gene.mean.by.age.group))
      {

        # New: get extracted expression features
        cur.gene.ind <- expr.stats$DF.expr.stats[[k]][, paste0(expression.stat.y, "_", age.group)] > -1  # Set current gene inds 
        names(cur.gene.ind) <- rownames(expr.stats$DF.expr.stats[[k]])

        if(sum(cur.gene.ind, na.rm=TRUE) < filter.params$min.genes.for.reg) # require minimal number of genes !!!
          next # skip this cell type and age group
        
        cur_expr_covariates = unlist( lapply(expression.stat.x, paste0, "_", age.group) )
        cur_reg_covariates_mat <- cbind(cur_gene_features_mat[which(cur.gene.ind),], 
                                        expr.stats$DF.expr.stats[[k]][which(cur.gene.ind),cur_expr_covariates])
        reg.model <- lm(expr.stats$DF.expr.stats[[k]][which(cur.gene.ind), paste0(expression.stat.y, "_", age.group)] ~ ., 
                        data = as.data.frame(cur_reg_covariates_mat))
#                          as.data.frame(cur_gene_features_mat[which(cur.gene.ind),]))  # Take log of fold-change. Maybe take difference? (they're after log)

        DF_cor[cell.type.ctr, beta.inds[((reg.ctr-1)*n.features+1):(reg.ctr*n.features)]] <- reg.model$coefficients[-1] # get coefficients (excluding intercept)
        DF_cor[cell.type.ctr, beta.pvals.inds[((reg.ctr-1)*n.features+1):(reg.ctr*n.features)]] <- summary(reg.model)$coefficients[-1,4]  # get p-values (excluding intercept?)
        reg.ctr = reg.ctr + 1
      }
      
      # Do regression for fold-change 
      fc.gene.ind <- expr.stats$DF.expr.stats[[k]][, paste0(expression.stat.y, "_fc")] > -1  # Set current gene inds 
      cur_reg_covariates_mat <- cbind(cur_gene_features_mat[which(fc.gene.ind),], 
                                      expr.stats$DF.expr.stats[[k]][which(fc.gene.ind),cur_expr_covariates])
      

      fc.reg.model <- lm(expr.stats$DF.expr.stats[[k]][which(fc.gene.ind), paste0(expression.stat.y, "_fc")]  ~ ., 
                         data = as.data.frame(cur_reg_covariates_mat))
#                         data = as.data.frame(cur_gene_features_mat[fc.gene.ind,]))  # Take log of fold-change. Maybe take difference? (they're after log)
      DF_cor[cell.type.ctr, beta.inds[((reg.ctr-1)*n.features+1):(reg.ctr*n.features)]] <- fc.reg.model$coefficients[-1] # get p-values (excluding intercept)
      DF_cor[cell.type.ctr, beta.pvals.inds[((reg.ctr-1)*n.features+1):(reg.ctr*n.features)]] <- summary(fc.reg.model)$coefficients[-1,4]  # get p-values (excluding intercept?)
      reg.ctr = reg.ctr + 1
      
      cell.type.ctr = cell.type.ctr+1  # update counter
    }  # end loop on cell-types in tissue
  }  # end loop on tissues
  
  save(DF_cor, file=reg.analysis.outfile)   # Save dataframe to file ( Save also to excel? ) 
  return(DF_cor)
} # End function mean expression droplet 
