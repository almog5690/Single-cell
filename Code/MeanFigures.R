library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)

source("scRNA_seq_utilities.R")

# Add saving figures to file in code!!! (use ggsave !!! )

draw_mean_figures <- function(data.types, fig.num, feature.types = c("selection", "gene.len"), 
                              tissue = "Lung", cell_type = "type II pneumocyte") {  # default tissue + cell type 
  n.datas <- length(data.types)
  p_feature_vs_mean_bar <- DF_cors <- vector("list", n.datas)
  for (i in 1:n.datas) {  # Either perform the analysis or read the output file 
    set_data_dirs(data.types[i])
    mean.analysis.outfile <- paste0(analysis.results.dir, 'mean.analysis.', data.types[i], # '.RData') 
                                    ".", paste0( feature.types, collapse="_"), '.RData') # should include also features 
    load(mean.analysis.outfile)
    DF_cors[[i]] <- DF_cor     #    DF_cors[[i]] <- mean_expression_analysis(data.types[i]) # don't run again 
  }
  
  if(fig.num %in% c(1, 11))  #### Figure 1 (all) or 11 (fold-change)
  {
    group.str = if(fig.num == 1) "all" else "fc"
    for(feature.type in feature.types) # Here plot each feature vs. mean expression, not just selection
    {
      cor.col <- paste0(feature.type, "_", group.str, "_cor")
      pval.col <- paste0(feature.type, "_", group.str, "_pval")
      # adding cell type variable
      for (i in 1:n.datas) {
        DF_cors[[i]]$CT = interaction(DF_cors[[i]]$Organs, DF_cors[[i]]$Cell_type, sep = ":")
        DF_cors[[i]]$fill.plot <- unlist(-log10(DF_cors[[i]][pval.col]))
        DF_cors[[i]]$y.plot <- reorder(DF_cors[[i]]$CT, unlist(DF_cors[[i]][cor.col]))
        DF_cors[[i]]$x.plot <- unlist(DF_cors[[i]][cor.col])

        
        
        # Joint multiple plots !!! 
        p_feature_vs_mean_bar[[i]] = ggplot(DF_cors[[i]], aes(y = y.plot, x = x.plot, fill = fill.plot)) + 
          geom_bar(stat = "identity") + theme_classic() +
          scale_fill_gradient2(high = "red", mid = "white",low = "blue",midpoint = -log10(0.05)) +
          labs(title = processed.files.str[data.types[i]], x = "correlation",y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
          theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),
                plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
        
        
      }
      # merging facs and droplet plots into figure 1
      multi.plot <- ggarrange(plotlist = p_feature_vs_mean_bar, nrow = 1, ncol = n.datas) # all in the same row 
      annotate_figure(multi.plot, top = text_grob(paste0("Mean ", group.str, " vs. ", feature.type, " correlation across all cell types"), 
                                                  color = "black", face = "bold", size = 14))
      ggsave(paste0(analysis.figures.dir, "Mean.Figure", fig.num, '.', paste(data.types, collapse = "_"), 
                    '.', feature.type, '.png'), height = 6,width = 9)  # Modify name to get figure  
    }    
  }  
  
  #### Figure 2
  if(fig.num == 2){
    if(is.null(tissue) | is.null(cell_type)){
      stop("Tissue or cell type not provided")
    }
    gene_features = read_gene_features(feature.types) 
    for(feature.type in feature.types) # Here plot each feature vs. mean expression, not just selection
    {
      gene_name <- names(gene_features[[feature.type]])
      p3 <- p_denst <- vector("list", n.datas)
      for(i in 1:n.datas) {   
        meta.data = get_meta_data(data.types[[i]])
        samples <- get_tissue_file_names(data.types[[i]])
        groups <- dataset_to_age_groups(data.types[[i]])
        
        # Choosing cell to highlight (default: Lung's type 2 pneumocyte)
        highlight_cell = which(DF_cors[[i]]$Organs == tissue & DF_cors[[i]]$Cell_type == cell_type)
        DF_cors[[i]]$highlight <- ifelse((1:nrow(DF_cors[[i]])) == highlight_cell, "highlight", "normal")

        cor.young.col <- paste0(feature.type, "_young_cor")
        cor.old.col <- paste0(feature.type, "_old_cor")
        DF_cors[[i]]$x.plot <- unlist(DF_cors[[i]][cor.young.col])
        DF_cors[[i]]$y.plot <- unlist(DF_cors[[i]][cor.old.col])
        textdf <- DF_cors[[i]][highlight_cell, ]
        
        print("p3")
        # mean expression vs. selection correlation coefficients for both old and young
        p3[[i]] = ggplot(DF_cors[[i]], aes(x.plot, y.plot)) + 
          geom_point(color = "blue") + 
          geom_text(data = textdf, aes(x = x.plot * 1.12, y = y.plot * 0.91, label = cell_type)) +
          geom_abline(slope = 1, intercept = 0, col = "red") +
          geom_point(data = textdf, aes(x = x.plot , y = y.plot), color = "orangered", fill = "orange", size = 3) +
          labs(title = processed.files.str[data.types[i]], x = "Young cor",y = "Old cor") +
          theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8), axis.title = element_text(size = 10),legend.position = "none")

        ## highlighted cell (Lung Pneumocyte cell) type Mean vs selection for both age groups
        tissue.ind = which(samples$organs == tissue) # Lung index
        if(data.type == "CR.Rat"){ # Convert to dummy variables 
          cell_types_categories = levels(SC$cell_types)
        } else
            cell_types_categories = meta.data[[i]]$cell_ontology_class  # the names of the different cell-types
        k = which(cell_types_categories == cell_type) # type II pneumocyte cell index
        print("Read SC")
        
        SC = readRDS(file = paste0(processed.data.dir, samples$organs[tissue.ind], ".", processed.files.str[data.types[[i]]], ".rds")) # Current (Lung) tissue seurat object
        list2env(tissue_to_age_inds(data.types[[i]], samples$organs[tissue.ind], groups, SC@meta.data), env=environment()) # set specific ages for all age groups in all datasets
        counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the Lng tissue
        SC_gene_name = toupper(rownames(SC)) # genes names in upper case letters
        if(length(SC_gene_name) != dim(counts.mat)[1]) { # For Rats, names are already in the matrix
          SC_gene_name = toupper(rownames(counts.mat)) 
        } 
        rownames(counts.mat) = SC_gene_name # make sure upper 

        if(data.type == "CR.Rat"){ # Convert to dummy variables 
          cell_types = as.numeric(SC@meta.data$cell_types) # Cell types vector
        } else
          cell_types = SC@meta.data$cell.ontology.class # Cell types vector

        gene_feat = gene_features[[feature.type]][gene_name %in% (SC_gene_name)] # filtering the selection score to genes that are found in the current tissue
        cur_gene_name = gene_name[gene_name %in% (SC_gene_name )] # the names of the filtered genes
        
        # genes mean vector
        gene_mean = rowMeans(counts.mat[,cell_types==k-1])
        gene_mean = gene_mean[cur_gene_name]
        
        # genes mean vector-old
        gene_mean_old = rowMeans(counts.mat[,(cell_types==k-1)&(!young.ind)])
        gene_mean_old = gene_mean_old[cur_gene_name]
        
        # genes mean vector-young
        gene_mean_young = rowMeans(counts.mat[,(cell_types==k-1)&(young.ind)])
        gene_mean_young = gene_mean_young[cur_gene_name]
        
        # Filtering genes with less then 10 reads in current cell type
        genes_ind = rowSums(SC@assays$RNA@counts[,cell_types==k-1]) > filter.params$min.count # 10
        names(genes_ind) = toupper(names(genes_ind))
        genes_ind = genes_ind[cur_gene_name]
        
        # filtering selection and old and young mean expressions
        gene_feat = gene_feat[genes_ind]
        gene_mean_old = gene_mean_old[genes_ind]
        gene_mean_young = gene_mean_young[genes_ind]
        
        # data-frame contains old and young genes mean expression and selection for the plots
        df = data.frame("Mean" = c(gene_mean_old, gene_mean_young),
                        "Selection" = c(gene_feat, gene_feat),
                        "Age" = rep(c("Old","Young"), each = length(gene_feat)))
        
        # ranking selection and old and young mean expressions
        g = na.omit(gene_feat)
        selc_rank = rank(g,ties.method = "average")
        mean_young_rank = rank(gene_mean_young[names(g)],ties.method = "average")
        mean_old_rank = rank(gene_mean_old[names(g)],ties.method = "average")
        
        # data frame contains old and young genes mean expression and selection ranks for the plots
        print("df4")
        df_4 = data.frame("Mean" = c(mean_old_rank,mean_young_rank),
                          "gene.feature" = c(selc_rank, selc_rank),
                          "Age" = rep(c("Old","Young"),each = length(selc_rank)))
        print("age names")
        # Should replace this with actual computed correlation. No need for p-value
        age_name = c("Old" = paste0("Old: ","\u03c1","=0.18"), # ,p<2.2e-16"),
                     "Young" = paste0("Young: ","\u03c1","=0.21")) # ,p<2.2e-16"))
        
        # getting the 2D density of selection and mean for the plots
        print("df4 density")
        df_4$density =  get_density(df_4$Mean, df_4$gene.feature, n = 100) # problem here: df_4 is empty!!!! 
        print("p_denst")
        
        # Selection rank vs mean expression rank for both young and old for the Lung Pneumocyte cell type
        p_denst[[i]] = ggplot(df_4) +
          geom_point(aes(x = Mean, y = gene.feature, fill = density),color = "white", 
                     alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
          scale_fill_gradientn(colors = matlab.like(100)) + 
          facet_wrap(~Age, labeller = labeller(Age = age_name)) + 
          geom_smooth(method = "lm",se = F,data = df_4, aes(x = Mean, y = gene.feature, color = Age)) +
          scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
          labs(title = paste(samples$organs[i],cell_types_categories[k], sep = ": "), x = "Mean expression rank", y = paste0(feature.type, " rank")) + 
          theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold")) +
          guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
      } # end loop on data types 
      
      # merging all datasets plots into figure 2
      multi.plot <- ggarrange(plotlist = c(rbind(p_denst, p3)), nrow = 2, ncol = n.datas) # all in the same row 
      annotate_figure(multi.plot, top = text_grob(paste0("Genes ", feature.type, " and mean correlation"), 
                                                  color = "black", face = "bold", size = 14))
      ggsave(paste0(analysis.figures.dir, "Mean.Figure", fig.num, '.', paste(data.types, collapse = "_"),
                    '.',  feature.type, '.png'), height = 6,width = 9)  # Modify name to get figure  
    }  # loop on explanatory features 
  }  # if figure 2
  #    ggsave(paste(analysis.figures.dir,"mean FC vs selection correlation across all cell types.png",sep = "/"),height = 6,width = 9)
  
}
