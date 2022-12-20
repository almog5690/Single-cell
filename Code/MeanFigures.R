library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)

source("scRNA_seq_utilities.R")

# Add saving figures to file in code!!! (use ggsave !!! )

draw_mean_figures <- function(DF_cor, DF_cor_drop, fig.num, tissue = NULL, cell_type = NULL) {
  if(fig.num == 1)  #### Figure 1
  {
    ## facs
    # adding cell type variable
    DF_cor$CT = interaction(DF_cor$Organs, DF_cor$Cell_type,sep = ":")
    
    p_selc_mean_bar = ggplot(DF_cor,aes(y = reorder(CT,mean_selc_cor),x = mean_selc_cor,fill = -log10(pval_all))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "facs",x = "correlation",y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    ## droplet
    # adding cell type variable
    DF_cor_drop$CT = interaction(DF_cor_drop$Organs,DF_cor_drop$Cell_type,sep = ":")
    
    p_selc_mean_bar_drop = ggplot(DF_cor_drop,aes(y = reorder(CT,mean_selc_cor),x = mean_selc_cor,fill = -log10(pval_all))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "droplet",x = "correlation",y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    # merging facs and droplet plots into figure 1
    title = ggdraw() + draw_label("Mean vs selection correlation across all cell types") # plot title
    p = plot_grid(p_selc_mean_bar,p_selc_mean_bar_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
    ggsave(paste(analysis.figure.dir,"mean vs selection correlation across all cell types.png",sep = "/"),height = 6,width = 9)
    
  }  
  
  #### Figure 2
  if(fig.num == 2){
    if(is.null(tissue) | is.null(cell_type)){
      stop("Tissue or cell type not provided")
    }
    ## facs
    # Choosing Lung's type 2 pneumocyte cell type to highlight
    highlight_cell = which(DF_cor$Organs == tissue & DF_cor$Cell_type == cell_type)
    DF_cor$highlight <- ifelse((1:nrow(DF_cor)) == highlight_cell, "highlight", "normal")
    textdf <- DF_cor[highlight_cell, ]
    mycolours <- c("highlight" = "green", "normal" = "blue")
    
    # facs mean expression vs selection correlation coefficients for both old and young
    p3 = ggplot(DF_cor,aes(selc_mean_young_cor_spearman,selc_mean_old_cor_spearman)) + 
      geom_point(color = "blue") + 
      geom_text(data = textdf, aes(x = selc_mean_young_spearman * 1.12, y = selc_mean_old_cor_spearman * 0.91, label = cell_type)) +
      geom_abline(slope = 1,intercept = 0,col = "red") +
      geom_point(data = textdf,aes(x = selc_mean_young_spearman , y = selc_mean_old_cor_spearman),color = "orangered",fill = "orange",size = 3) +
      labs(title = "facs",x = "Young cor",y = "Old cor") +
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),legend.position = "none")
    
    ## facs - Lung Pneumocyte cell type Mean vs selection for both age groups
    i = which(organs == tissue) # Lung index
    cell_types_categories = meta.data[[i]]$cell_ontology_class # the names of the different cell-types
    k = which(cell_types_categories == cell_type) # type II pneumocyte cell index
    
    SC = readRDS(file = paste(processed.data.dir,paste(organs[i],"rds",sep = "."),sep = "/")) # Lung tissue seurat object
    counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the Lng tissue
    young.ind = c(SC@meta.data$age %in% c("3m")) # index for cells that came from young mouses
    SC_gene_name = toupper(rownames(SC)) # genes names in upper case letters
    rownames(counts.mat) = SC_gene_name 
    
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
    
    # data frame contains old and young genes mean expression and selection for the plots
    df = data.frame("Mean" = c(gene_mean_old,gene_mean_young),"Selection" = c(gene_selc,gene_selc),"Age" = rep(c("Old","Young"),each = length(gene_selc)))
    
    # ranking selection and old and young mean expressions
    g = na.omit(gene_selc)
    selc_rank = rank(g,ties.method = "average")
    mean_young_rank = rank(gene_mean_young[names(g)],ties.method = "average")
    mean_old_rank = rank(gene_mean_old[names(g)],ties.method = "average")
    
    # data frame contains old and young genes mean expression and selection ranks for the plots
    df_4 = data.frame("Mean" = c(mean_old_rank,mean_young_rank),"Selection" = c(selc_rank,selc_rank),"Age" = rep(c("Old","Young"),each = length(selc_rank)))
    
    age_name = c("Old" = paste0("Old: ","\u03c1","=0.18,p<2.2e-16"),"Young" = paste0("Young: ","\u03c1","=0.21 ,p<2.2e-16"))
    
    # getting the 2D density of selection and mean for the plots
    df_4$density =  get_density(df_4$Mean, df_4$Selection, n = 100)
    # Selection rank vs mean expression rank for both young and old for the Lung Pneumocyte cell type
    p_denst = ggplot(df_4) +
      geom_point(aes(x = Mean, y = Selection, fill = density),color = "white", alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
      scale_fill_gradientn(colors = matlab.like(100)) + 
      facet_wrap(~Age,labeller = labeller(Age = age_name)) + 
      geom_smooth(method = "lm",se = F,data = df_4,aes(x = Mean, y = Selection,color = Age)) +
      scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
      labs(title = paste(organs[i],cell_types_categories[k],sep = ": "),x = "Mean expression rank",y = "Selection rank") + 
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold")) +
      guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
    
    
    ## droplet
    # Choosing Lung's type 2 pneumocyte cell type to highlight
    highlight_cell = which(DF_cor_drop$Organs == tissue & DF_cor_drop$Cell_type == cell_type)
    DF_cor_drop$highlight <- ifelse((1:nrow(DF_cor_drop)) == highlight_cell, "highlight", "normal")
    textdf <- DF_cor_drop[highlight_cell, ]
    mycolours <- c("highlight" = "green", "normal" = "blue")
    
    # droplet mean expression vs selection correlation coefficients for both old and young
    p3_drop = ggplot(DF_cor_drop,aes(selc_mean_young_spearman,selc_mean_old_cor_spearman)) + 
      geom_point(color = "blue") +
      geom_text(data = textdf, aes(x = selc_mean_young_spearman * 1.03, y = selc_mean_old_cor_spearman * 0.89, label = cell_type)) +
      geom_abline(slope = 1,intercept = 0,col = "red") +
      geom_point(data = textdf,aes(x = selc_mean_young_spearman , y = selc_mean_old_cor_spearman),color = "orangered",fill = "orange",size = 3) +
      labs(title = "droplet",x = "Young cor",y = "Old cor") +
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),legend.position = "none")
    
    
    ## droplet - Lung Pneumocyte cell type Mean vs selection for both age groups
    i = which(drop_organs == tissue) # Lung index
    cell_types_categories = meta.data.drop[[i]]$cell_ontology_class # the names of the different cell-types
    k = which(cell_types_categories == cell_type) # type II pneumocyte cell index
    old_ages_1 = c("21m","24m") # The ages of the old mice
    
    SC = readRDS(file = paste(drop_organs[i],"rds",sep = ".")) # Lung tissue seurat object
    counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the Lng tissue
    young.ind = c(SC@meta.data$age %in% c("3m")) # index for cells that came from young mouses
    old.ind = c(SC@meta.data$age %in% old_ages_1) # index for cells that came from old mouses
    SC_gene_name = toupper(rownames(SC)) # genes names in upper case letters
    rownames(counts.mat) = SC_gene_name 
    
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
    
    # data frame contains old and young genes mean expression and selection for the plots
    df = data.frame("Mean" = c(gene_mean_old,gene_mean_young),"Selection" = c(gene_selc,gene_selc),"Age" = rep(c("Old","Young"),each = length(gene_selc)))
    
    # ranking selection and old and young mean expressions
    g = na.omit(gene_selc)
    selc_rank = rank(g,ties.method = "average")
    mean_young_rank = rank(gene_mean_young[names(g)],ties.method = "average")
    mean_old_rank = rank(gene_mean_old[names(g)],ties.method = "average")
    
    # data frame contains old and young genes mean expression and selection ranks for the plots
    df_4 = data.frame("Mean" = c(mean_old_rank,mean_young_rank),"Selection" = c(selc_rank,selc_rank),"Age" = rep(c("Old","Young"),each = length(selc_rank)))
    
    age_name = c("Old" = paste0("Old: ","\u03c1","=0.093,p=1.2e-13"),"Young" = paste0("Young: ","\u03c1","=0.16 ,p<2.2e-16"))
    
    # getting the 2D density of selection and mean for the plots
    df_4$density =  get_density(df_4$Mean, df_4$Selection, n = 100)
    # Selection rank vs mean expression rank for both young and old for the Lung Pneumocyte cell type
    p_denst_drop = ggplot(df_4) +
      geom_point(aes(x = Mean, y = Selection, fill = density),color = "white", alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
      scale_fill_gradientn(colors = matlab.like(100)) + 
      facet_wrap(~Age,labeller = labeller(Age = age_name)) + 
      geom_smooth(method = "lm",se = F,data = df_4,aes(x = Mean, y = Selection,color = Age)) +
      scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
      labs(title = paste(drop_organs[i],cell_types_categories[k],sep = ": "),x = "Mean expression rank",y = "Selection rank") + 
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold")) +
      guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
    
    # merging facs and droplet plots into figure 2
    title = ggdraw() + draw_label("Genes selection and mean correlation") # plot title
    p = plot_grid(p_denst,p3,p_denst_drop,p3_drop,labels = LETTERS[1:4]) 
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
    ggsave(paste(analysis.figure.dir,"selection and mean corr gene-filter density and cor plot.png",sep = "/"),height = 6,width = 9)
    
  }
  
  if(fig.num == 11) #### Figure 11
  {
    ## facs
    # adding cell type variable
    DF_cor$CT = interaction(DF_cor$Organs,DF_cor$Cell_type,sep = ":")
    
    p_selc_mean_fc_bar = ggplot(DF_cor,aes(y = reorder(CT,fc_cor),x = fc_cor,fill = -log10(pval_fc))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "facs",x = "correlation",y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    
    ## droplet
    # adding cell type variable
    DF_cor_drop$CT = interaction(DF_cor_drop$Organs,DF_cor_drop$Cell_type,sep = ":")
    
    p_selc_mean_fc_bar_drop = ggplot(DF_cor_drop,aes(y = reorder(CT,fc_cor),x = fc_cor,fill = -log10(pval_fc))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "droplet",x = "correlation",y = "Cell types",fill = "-log10(P_val)") +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    # merging facs and droplet plots into figure 11
    title = ggdraw() + draw_label("Mean FC vs selection correlation across all cell types") # plot title
    p = plot_grid(p_selc_mean_fc_bar,p_selc_mean_fc_bar_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
    ggsave(paste(analysis.figure.dir,"mean FC vs selection correlation across all cell types.png",sep = "/"),height = 6,width = 9)
    
  }  
}
