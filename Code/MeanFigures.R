library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)

source("scRNA_seq_utilities.R")

# Add saving figures to file in code!!! (use ggsave !!! )

draw_mean_figures <- function(data.types, fig.num, tissue = "Lung", cell_type = "type II pneumocyte") {  # default tissue + cell type 
  n.datas <- length(data.types)
  DF_cors <- vector("list", n.datas)
  p_selc_mean_bar <- vector("list", n.datas)
  for (i in 1:n.datas) {  # Either perform the analysis or read the output file 
    mean.analysis.outfile <- paste0(analysis.results.dir, 'mean.analysis.', data.types[i], '.RData') 
    #                                    paste0( feature.types, collapse="_"), '.RData') # should include also features 
    load(mean.analysis.outfile)
    DF_cors[[i]] <- DF_cor     #    DF_cors[[i]] <- mean_expression_analysis(data.types[i]) # don't run again 
  }
  
  if(fig.num == 1)  #### Figure 1
  {
    ## facs
    # adding cell type variable
    for (i in 1:n.datas) {
      DF_cors[[i]]$CT = interaction(DF_cors[[i]]$Organs, DF_cors[[i]]$Cell_type,sep = ":")
      # Joint multiple plots !!! 
      p_selc_mean_bar[[i]] = ggplot(DF_cors[[i]],aes(y = reorder(CT,mean_selc_cor),x = mean_selc_cor, fill = -log10(pval_all))) + 
        geom_bar(stat = "identity") + 
        theme_classic() +
        scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
        labs(title = processed.files.str[data.types[i]], x = "correlation",y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
        theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),
              plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    }
    
    # merging facs and droplet plots into figure 1
    multi.plot <- ggarrange(plotlist = p_selc_mean_bar,nrow = 1, ncol = n.datas) # all in the same row 
    annotate_figure(multi.plot, top = text_grob("Mean vs. selection correlation across all cell types", 
                                                color = "black", face = "bold", size = 14))
    
    
  }  
  
  #### Figure 2
  if(fig.num == 2){
    if(is.null(tissue) | is.null(cell_type)){
      stop("Tissue or cell type not provided")
    }
    
    gene_features = read_gene_features(c("selection"))  # Problem: we need the same number of genes !!! 
    gene_selection = gene_features[[1]]
    gene_name <- names(gene_selection)
    
    p3 <- p_denst <- vector("list", n.datas)
    for(i in 1:n.datas) {   
      meta.data = get_meta_data(data.types[[i]])
      samples <- get_tissue_file_names(data.types[[i]])
      

      # Choosing cell to highlight (default: Lung's type 2 pneumocyte)
      highlight_cell = which(DF_cors[[i]]$Organs == tissue & DF_cors[[i]]$Cell_type == cell_type)
      DF_cors[[i]]$highlight <- ifelse((1:nrow(DF_cors[[i]])) == highlight_cell, "highlight", "normal")
      textdf <- DF_cors[[i]][highlight_cell, ]
      mycolours <- c("highlight" = "green", "normal" = "blue")
      
      # mean expression vs selection correlation coefficients for both old and young
      p3[[i]] = ggplot(DF_cors[[i]], aes(selc_mean_young_spearman, selc_mean_old_cor_spearman)) + 
        geom_point(color = "blue") + 
        geom_text(data = textdf, aes(x = selc_mean_young_spearman * 1.12, y = selc_mean_old_cor_spearman * 0.91, label = cell_type)) +
        geom_abline(slope = 1,intercept = 0,col = "red") +
        geom_point(data = textdf,aes(x = selc_mean_young_spearman , y = selc_mean_old_cor_spearman),color = "orangered",fill = "orange",size = 3) +
        labs(title = processed.files.str[data.types[i]], x = "Young cor",y = "Old cor") +
        theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),legend.position = "none")
      
      ## highlighted cell (Lung Pneumocyte cell) type Mean vs selection for both age groups
      tissue.ind = which(samples$organs == tissue) # Lung index
      cell_types_categories = meta.data[[tissue.ind]]$cell_ontology_class # the names of the different cell-types
      k = which(cell_types_categories == cell_type) # type II pneumocyte cell index
      
      SC = readRDS(file = paste0(processed.data.dir, samples$organs[tissue.ind], ".", processed.files.str[data.type], ".rds")) # Current (Lung) tissue seurat object
      print(paste0("Read SC ", i))
      
      counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the Lng tissue
      young.ind = c(SC@meta.data$age %in% c("3m")) # index for cells that came from young mice
      SC_gene_name = toupper(rownames(SC)) # genes names in upper case letters
      rownames(counts.mat) = SC_gene_name 
      cell_types = SC@meta.data$cell.ontology.class # Cell types vector
      
      
      gene_selc = gene_selection[gene_name %in% (SC_gene_name)] # filtering the selection score to genes that are found in the current tissue
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
      genes_ind = rowSums(SC@assays$RNA@counts[,cell_types==k-1]) > 10
      names(genes_ind) = toupper(names(genes_ind))
      genes_ind = genes_ind[cur_gene_name]
      
      # filtering selection and old and young mean expressions
      gene_selc = gene_selc[genes_ind]
      gene_mean_old = gene_mean_old[genes_ind]
      gene_mean_young = gene_mean_young[genes_ind]
      print(paste0("Data-frame ", i))
      
      # data frame contains old and young genes mean expression and selection for the plots
      df = data.frame("Mean" = c(gene_mean_old,gene_mean_young),"Selection" = c(gene_selc,gene_selc),"Age" = rep(c("Old","Young"),each = length(gene_selc)))
      
      # ranking selection and old and young mean expressions
      g = na.omit(gene_selc)
      selc_rank = rank(g,ties.method = "average")
      mean_young_rank = rank(gene_mean_young[names(g)],ties.method = "average")
      mean_old_rank = rank(gene_mean_old[names(g)],ties.method = "average")
      
      # data frame contains old and young genes mean expression and selection ranks for the plots
      df_4 = data.frame("Mean" = c(mean_old_rank,mean_young_rank),"Selection" = c(selc_rank,selc_rank),"Age" = rep(c("Old","Young"),each = length(selc_rank)))
      
      # Should replace this with actual computed correlation. No need for p-value
      age_name = c("Old" = paste0("Old: ","\u03c1","=0.18"), # ,p<2.2e-16"),
                   "Young" = paste0("Young: ","\u03c1","=0.21")) # ,p<2.2e-16"))
      
      # getting the 2D density of selection and mean for the plots
      df_4$density =  get_density(df_4$Mean, df_4$Selection, n = 100)
      print(paste0("ggplot-again ", i))
      
      # Selection rank vs mean expression rank for both young and old for the Lung Pneumocyte cell type
      p_denst[[i]] = ggplot(df_4) +
        geom_point(aes(x = Mean, y = Selection, fill = density),color = "white", 
                   alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
        scale_fill_gradientn(colors = matlab.like(100)) + 
        facet_wrap(~Age,labeller = labeller(Age = age_name)) + 
        geom_smooth(method = "lm",se = F,data = df_4,aes(x = Mean, y = Selection,color = Age)) +
        scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
        labs(title = paste(samples$organs[i],cell_types_categories[k],sep = ": "),x = "Mean expression rank", y = "Selection rank") + 
        theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold")) +
        guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
    } # end loop on data types 
    
    # merging all datasets plots into figure 2
    multi.plot <- ggarrange(plotlist = c(rbind(p_denst, p3)), nrow = 2, ncol = n.datas) # all in the same row 
    annotate_figure(multi.plot, top = text_grob("Genes selection and mean correlation", 
                                                color = "black", face = "bold", size = 14))
    
    
    
#    title = ggdraw() + draw_label("Genes selection and mean correlation") # plot title
#    p = plot_grid(p_denst,p3,p_denst_drop,p3_drop,labels = LETTERS[1:4]) 
#    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
#    ggsave(paste(analysis.figures.dir,"selection and mean corr gene-filter density and cor plot.png",sep = "/"),height = 6,width = 9)
    
  }
  
  if(fig.num == 11) #### Figure 11
  {
    ## facs
    # adding cell type variable
    DF_cors[[1]]$CT = interaction(DF_cors[[1]]$Organs,DF_cors[[1]]$Cell_type,sep = ":")
    
    p_selc_mean_fc_bar = ggplot(DF_cors[[1]],aes(y = reorder(CT, fc_cor),x = fc_cor,fill = -log10(pval_fc))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = processed.files.str[data.types[1]],x = "correlation",y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    ## droplet
    # adding cell type variable
    DF_cors[[2]]$CT = interaction(DF_cors[[2]]$Organs,DF_cors[[2]]$Cell_type,sep = ":")
    
    p_selc_mean_fc_bar_drop = ggplot(DF_cors[[2]],aes(y = reorder(CT,fc_cor),x = fc_cor,fill = -log10(pval_fc))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = processed.files.str[data.types[2]],x = "correlation",y = "Cell types",fill = "-log10(P_val)") +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    # merging facs and droplet plots into figure 11
    title = ggdraw() + draw_label("Mean FC vs selection correlation across all cell types") # plot title
    p = plot_grid(p_selc_mean_fc_bar,p_selc_mean_fc_bar_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  } 
  ggsave(paste0(analysis.figures.dir, "Mean.Figure", fig.num, '.png'), height = 6,width = 9)  # Modify name to get figure  
  #    ggsave(paste(analysis.figures.dir,"mean FC vs selection correlation across all cell types.png",sep = "/"),height = 6,width = 9)
  
}
