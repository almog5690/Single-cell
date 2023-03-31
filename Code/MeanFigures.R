library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)
library(reshape2)
library(ggcorrplot)
library(stringr)

source("scRNA_seq_utilities.R")

# Add saving figures to file in code!!! (use ggsave !!! )
draw_expr_reg_figures <- function(data.types, expr.stat.y = "mean", expr.stat.x = c(), fig.num, feature.types = c("selection", "gene.len"), 
                              tissue = "Lung", cell_type = "type II pneumocyte", n.features = 2) {  # default tissue + cell type 
  n.datas <- length(data.types)
  p_feature_vs_mean_bar <- DF_cors <- vector("list", n.datas)
  num.cell.types <- rep(0, n.datas)
  for (i in 1:n.datas) {  # Either perform the analysis or read the output file 
    set_data_dirs(data.types[i])
    expr.reg.analysis.outfile <- paste0(analysis.results.dir, 'expr.reg.analysis.', data.types[i], # '.RData') 
                                    ".", expr.stat.y, ".vs.", 
                                    paste0( expr.stat.x, collapse="_"), rep('.', min(length(expr.stat.x), 1)), 
                                    paste0( feature.types, collapse="_"), '.RData') # should include also features 
    print(expr.reg.analysis.outfile)
    load(expr.reg.analysis.outfile)
    DF_cors[[i]] <- DF_cor     #    DF_cors[[i]] <- mean_expression_analysis(data.types[i]) # don't run again 
    num.cell.types[i] = dim(DF_cors[[i]])[1]
  }

  if(fig.num %in% c(1, 11, 111))  #### Figure 1 (all) or 11 (fold-change) or 111 (abs-fold-change)
    draw_cor_bars_figure(fig.num, data.types, feature.types, DF_cors, analysis.figures.dir)

  #### Figure 2: Detailed correlations plot
  if(fig.num %in% c(2, 22)){
    if(is.null(tissue) | is.null(cell_type)){
      stop("Tissue or cell type not provided")
    }
    draw_cor_scatters_figure(fig.num, data.types, feature.types, DF_cors, analysis.figures.dir, num.cell.types, 
                                         tissue = "Lung", cell_type = "type II pneumocyte")
  }  # end if figure 2 or 22
  if(fig.num %in% c(66,666))  # show overview boxplots 
    draw_boxplot_cor_overview_figure(fig.num, data.types, feature.types, DF_cors, analysis.figures.dir, n.features)
  

  if(fig.num == 99)  # New: draw heatmap of correlations of all gene features. Need to compute all pairwise correlations across tissues and cell types 
    draw_heatmap_cor_figure(fig.num, data.types, feature.types, analysis.figures.dir,  tissue = "Lung")
  
  #    ggsave(paste(analysis.figures.dir,"mean FC vs selection correlation across all cell types.png",sep = "/"),height = 6,width = 9)
}



draw_heatmap_cor_figure <- function(fig.num, data.types, feature.types, analysis.figures.dir, tissue = "Lung")
{
  print("Read gene features:")
  gene.features <- read_gene_features(feature.names=feature.types)
  print("Extract expr. stat:")
  expr.stats <- extract_expression_statistics(data.types[1], tissue, expression.stats = "mean") # extract means
  n.cell.types <- length(expr.stats$cells_ind)
  mean.expr.list = vector("list", n.cell.types)
  print("Loop on cell types:")
  for(i in 1:n.cell.types)
  {
    mean.expr.list[[i]] = expr.stats$DF.expr.stats[[expr.stats$cells_ind[i]]][,"mean_all"]
    names(mean.expr.list[[i]]) = rownames(expr.stats$DF.expr.stats[[expr.stats$cells_ind[i]]])
  }
#  print("Names:")
  names(mean.expr.list) <- expr.stats$cell_types_categories[expr.stats$cells_ind]
#  print("list.to::")
  df.features.and.mean.expr <- list_to_common_dataframe(c(mean.expr.list, gene.features))
  
#  print("Dim:")
#  print(dim(df.features.and.mean.expr))
#  print("cor::")
#  save("df.features.and.mean.expr", file="tmp_cor.RData")
  features.and.mean.expr.cor.mat = cor(df.features.and.mean.expr, use = "complete.obs") 
#  print("corrplot::")
  
  ggcorrplot(features.and.mean.expr.cor.mat)
  ggsave(paste0(analysis.figures.dir, "Mean.Figure", fig.num, '.', paste(data.types[1], collapse = "_"), 
                '.', paste(feature.types, collapse = "_"), '.', tissue, '.png'), height = 6,width = 9)  # Modify name to get figure  
}

# Draw figure of bars of correlation or beta regression coefficients with color showing log p-values
draw_cor_bars_figure <- function(fig.num, data.types, feature.types, DF_cors, analysis.figures.dir)
{
  n.datas <- length(data.types)
  p_feature_vs_mean_bar <- vector("list", n.datas)
  
  group.str = switch(as.character(fig.num), "1" = "all", "11" = "fc", "111" = "fc_abs")  #  == 1) "all" else "fc"
  for(feature.type in feature.types) # Here plot each feature vs. mean expression, not just selection
  {
    cor.col <- paste0(feature.type, "_", group.str, "_cor")
    pval.col <- paste0(feature.type, "_", group.str, "_pval")
    # adding cell type variable
    for (i in 1:n.datas) {
      if(!(cor.col %in% colnames(DF_cors[[i]])) | !(pval.col %in% colnames(DF_cors[[i]]))) # missing data
      {
        print(paste0("Error! Missing ", cor.col, " ", pval.col, " data! aborting!"))
        return(-1)
      }
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


# Detailed regression of one cell type with two-dim density, and other young vs. old. scatter plots 
# Plot colorful scatter 
draw_cor_scatters_figure <- function(fig.num, data.types, feature.types, DF_cors, analysis.figures.dir, num.cell.types,
                                     tissue = "Lung", cell_type = "type II pneumocyte")
{
  n.datas <- length(data.types)
  gene_features = read_gene_features(feature.types) 
  for(feature.type in feature.types) # Here plot each feature vs. mean expression, not just selection
  {
    gene_name <- names(gene_features[[feature.type]])
    p3 <- p_denst <- vector("list", n.datas)
    for(i in 1:n.datas) {   
      meta.data = get_meta_data(data.types[i])
      samples <- get_tissue_file_names(data.types[i])
      groups <- dataset_to_age_groups(data.types[i])
      
      # Choosing cell to highlight (default: Lung's type 2 pneumocyte)
      highlight_cell = which(DF_cors[[i]]$Organs == tissue & DF_cors[[i]]$Cell_type == cell_type)
      DF_cors[[i]]$highlight <- ifelse((1:nrow(DF_cors[[i]])) == highlight_cell, "highlight", "normal")
      
      if(fig.num == 2)  # Spearman correlations
      {
        fig.str = "cor"
      } else # beta coefficients of multiple regression 
      {
        fig.str = "beta"
      } 
      cor.young.col <- paste0(feature.type, "_young_", fig.str)
      cor.old.col <- paste0(feature.type, "_old_", fig.str)
      DF_cors[[i]]$x.plot <- unlist(DF_cors[[i]][cor.young.col])
      DF_cors[[i]]$y.plot <- unlist(DF_cors[[i]][cor.old.col])
      textdf <- DF_cors[[i]][highlight_cell, ]
      
      # Signed binomial test
      num.old.bigger.young = sum(DF_cors[[i]]$y.plot > DF_cors[[i]]$x.plot)
      num.old.bigger.young.pval = min(1, 2 * pbinom(min(num.old.bigger.young, num.cell.types[i]-num.old.bigger.young), 
                                                    num.cell.types[i], 0.5))  # two sided test 
      
      # mean expression vs. selection correlation coefficients for both old and young
      p3[[i]] = ggplot(DF_cors[[i]], aes(x.plot, y.plot)) + 
        geom_point(color = "blue") + 
        geom_text(data = textdf, aes(x = x.plot * 1.12, y = y.plot * 0.91, label = cell_type)) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        geom_point(data = textdf, aes(x = x.plot , y = y.plot), color = "orangered", fill = "orange", size = 3) +
        labs(title = paste0(processed.files.str[data.types[i]], 
                            ", Sign test: ", num.old.bigger.young, "/", num.cell.types[i], 
                            " Pval=", signif(num.old.bigger.young.pval, 3)), 
             x = paste0("Young ", fig.str), y = paste0("Old ", fig.str)) +
        theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8), axis.title = element_text(size = 10),legend.position = "none")
      
      ## highlighted cell (Lung Pneumocyte cell) type Mean vs selection for both age groups
      tissue.ind = which(samples$organs == tissue) # Lung index
      if(data.types[i] == "CR.Rat"){ # Convert to dummy variables 
        cell_types_categories = levels(SC$cell_types)
      } else
        cell_types_categories = meta.data[[tissue.ind]]$cell_ontology_class  # the names of the different cell-types
      k = which(cell_types_categories == cell_type) # type II pneumocyte cell index
      SC = readRDS(file = paste0(processed.data.dir, samples$organs[tissue.ind], ".", processed.files.str[data.types[i]], ".rds")) # Current (Lung) tissue seurat object
      list2env(tissue_to_age_inds(data.types[i], samples$organs[tissue.ind], groups, SC@meta.data), env=environment()) # set specific ages for all age groups in all datasets
      counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the Lung tissue
      SC_gene_name = toupper(rownames(SC)) # genes names in upper case letters
      if(length(SC_gene_name) != dim(counts.mat)[1]) { # For Rats, names are already in the matrix
        SC_gene_name = toupper(rownames(counts.mat)) 
      } 
      rownames(counts.mat) = SC_gene_name # make sure upper 
      
      if(data.types[i] == "CR.Rat"){ # Convert to dummy variables 
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
      
      # Filtering genes with less then 10 reads in ALL cells of current cell type
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
      selc_rank = rank(g, ties.method = "average")
      mean_young_rank = rank(gene_mean_young[names(g)], ties.method = "average")
      mean_old_rank = rank(gene_mean_old[names(g)], ties.method = "average")
      
      # data frame contains old and young genes mean expression and selection ranks for the plots
      df_4 = data.frame("Mean" = c(mean_old_rank, mean_young_rank),
                        "gene.feature" = c(selc_rank, selc_rank),
                        "Age" = rep(c("Old","Young"),each = length(selc_rank)))
      # Display correlation (no need for p-value)
      age_name = c("Old" = paste0("Old: ","\u03c1","=", round(DF_cors[[i]][highlight_cell,cor.old.col], 3)), # ,p<2.2e-16"),
                   "Young" = paste0("Young: ","\u03c1","=", round(DF_cors[[i]][highlight_cell,cor.young.col], 3))) # ,p<2.2e-16"))
      
      # getting the 2D density of selection and mean for the plots
      df_4$density =  get_density(df_4$Mean, df_4$gene.feature, n = 100) # problem here: df_4 is empty!!!! 
      
      
      # Selection rank vs mean expression rank for both young and old for the Lung Pneumocyte cell type
      p_denst[[i]] = ggplot(df_4) +
        geom_point(aes(x = Mean, y = gene.feature, fill = density), color = "white", 
                   alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
        scale_fill_gradientn(colors = matlab.like(100)) + 
        facet_wrap(~Age, labeller = labeller(Age = age_name)) + 
        geom_smooth(method = "lm",se = F,data = df_4, aes(x = Mean, y = gene.feature, color = Age)) +
        scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
        labs(title = paste(samples$organs[i],cell_types_categories[k], sep = ": "), 
             x = "Mean expression rank", y = paste0(feature.type, " rank")) + 
        theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),
              axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold")) +
        guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
    } # end loop on data types 
    
    # merging all datasets plots into figure 2
    multi.plot <- ggarrange(plotlist = c(rbind(p_denst, p3)), nrow = 2, ncol = n.datas) # all in the same row 
    annotate_figure(multi.plot, top = text_grob(paste0("Genes ", feature.type, " and mean correlation"), 
                                                color = "black", face = "bold", size = 14))
    
    if(fig.num == 2)  # marginal correlation
      ggsave(paste0(analysis.figures.dir, "Mean.Figure", fig.num, '.', paste(data.types, collapse = "_"),
                  '.',  feature.type, '.png'), height = 6,width = 9)  # Modify name to get figure  
    else # regression coefficient, context matters
      ggsave(paste0(analysis.figures.dir, "Mean.Figure", fig.num, '.', paste(data.types, collapse = "_"),
                    '.', paste(feature.types, collapse = "_"), 
                    '._beta_', feature.type, '.png'), height = 6,width = 9)  # Modify name to get figure  
    
  }  # loop on explanatory features 
} # End function for plotting fig. 2,22


# new: Boxplots for all correlations: 
# First top will be spearman correlation, second bottom will be beta regression coefficient
# Separate figures for mean and overdispersion
# 1. selection
# 2. length
#
draw_boxplot_cor_overview_figure <- function(fig.num, data.types, feature.types, DF_cors, analysis.figures.dir, n.features = 2)
{
  expr.stat.y <- if (fig.num == 66) "mean" else "overdispersion"
  n.datas <- length(data.types)
  meta.data <- vector("list", n.datas)
  names(meta.data) <- data.types
  num.cell.types <- rep(0, n.datas)
  for (i in 1:n.datas) {  # Either perform the analysis or read the output file 
    set_data_dirs(data.types[i])
    num.cell.types[i] = dim(DF_cors[[i]])[1]
    meta.data[[data.type]] = get_meta_data(data.types[i])
  }
  
  
#  n.features <- 5 #  n.features <- length(feature.types[1:5])
  all.show.columns <- c()
  for(s in feature.types[1:n.features])
    for(age.group in c(age.groups, "fc", "deltaYO"))
      for(stat.type in c("cor", "beta"))
        all.show.columns <- c(all.show.columns, paste0(s, "_", age.group, "_", stat.type))  

  
  df.s <- vector("list", n.datas)
  for (i in 1:n.datas)   # Either perform the analysis or read the output file 
  {
    for(s in feature.types[1:n.features])
      for(stat.type in c("cor", "beta"))
        DF_cors[[i]][, paste0(s, "_deltaYO_", stat.type)] <- DF_cors[[i]][, paste0(s, "_old_", stat.type)] - 
          DF_cors[[i]][, paste0(s, "_young_", stat.type)]
    # Create a new column for difference
    df.s[[i]] <- stack(DF_cors[[i]][1:40,all.show.columns])  %>% mutate(dataset = data.types[[i]])
  }
  df.box <- df.s[[1]]
  for (i in 2:n.datas)
    df.box <- rbind(df.box, df.s[[i]])
#  df.box <- df.box %.% mutate(covariate = strsplit(as.character(ind), "_")[[1]][1])
  df.box$covariate <- str_extract(as.character(df.box$ind),  "[^_]+") 
#  df.box$ind.only <- str_extract(as.character(df.box$ind),  "[_^]+")  # (?<=_)
  df.box$ind.only <- as.factor(sub("^[^_]*_", "", as.character(df.box$ind)))
  df.box$stat.type <- as.factor(sub("^[^_]*_", "", as.character(df.box$ind.only)))
  df.box$covariate.and.stat <- paste(df.box$covariate, df.box$stat.type)
  df.box$ind.only <- str_extract(as.character(df.box$ind.only),  "[^_]+")
  
  df.box$ind <- as.factor(df.box$ind) 
  facet.levels <- unique(df.box$covariate.and.stat)[c(seq(1, 2*n.features, 2), seq(2, 2*n.features, 2))]
    ggplot(df.box, aes(x = ind.only, y = values)) + 
      geom_boxplot(aes(fill = dataset)) + 
      facet_wrap( ~ factor(covariate.and.stat, levels = facet.levels), 
                  scales = "free", ncol = n.features) + # , labeller = label_both) +
      xlab(paste0(expr.stat.y, " analysis")) + 
    theme(axis.text.x = element_text(angle = 90,hjust=0))

  ggsave(paste0(analysis.figures.dir, "Mean.Figure", fig.num, '.', paste(data.types, collapse = "_"),
                '.', paste(feature.types, collapse = "_"), 
                '.png'), height = 6,width = 9)  # Modify name to get figure  
  
}

signif_plot <- function(beta_young,beta_old,p_val_young,p_val_old,p_val_diff
                        ,corr = TRUE,data.type = "facs")
{
  DF_signif = data.frame(beta_young,beta_old,"corr_diff_pval" = p_val_diff)
  
  x = p_val_old
  y = p_val_young
  p_thresh = 0.1
  DF_signif$which_signif = case_when(x <= p_thresh & y <= p_thresh ~ "Both",x > p_thresh & y <= p_thresh ~ "Young",
                                            x <= p_thresh & y > p_thresh ~ "Old",x > p_thresh & y > p_thresh ~ "Neither")
  
  
  DF_signif$cor_diff_pval_bh = p.adjust(DF_signif$corr_diff_pval,method = "BH")
  DF_signif$cor_diff_signif = ifelse(DF_signif$cor_diff_pval_bh < p_thresh,"Significant","Not_significant")
  
  if(corr)
  {
    p_signif = ggplot(DF_signif,aes(beta_young,beta_young)) + 
      geom_point(aes(shape = which_signif,size = which_signif,color = cor_diff_signif)) + 
      geom_abline(slope = 1,intercept = 0,col = "black",linetype = "dashed") +
      labs(title = data.type,x = "Young cor",y = "Old cor") +
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) + 
      scale_shape_manual(values = c("Both" = "+","Neither" = "•","Old" = "|","Young" = "—"),
                         name = TeX(r'($\rho$ Signif.)')) +
      scale_size_manual(values = c(4.5,5,2.7,2),name = TeX(r'($\rho$ Signif.)')) +
      scale_color_manual(values = c("gray50","red"),name = TeX(r'($\Delta\rho$ Signif.)'),breaks = c("Not_significant","Significant") ,labels = c("False","True")) + 
      geom_hline(yintercept = 0,color = "black") +
      geom_vline(xintercept = 0,color = "black")
  } else 
  {
    p_signif = ggplot(DF_signif,aes(beta_young,beta_young)) + 
      geom_point(aes(shape = which_signif,size = which_signif,color = cor_diff_signif)) + 
      geom_abline(slope = 1,intercept = 0,col = "black",linetype = "dashed") +
      labs(title = data.type,x = TeX(r'(Young $\beta$)'),y = TeX(r'(Old $\beta$)')) +
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) + 
      scale_shape_manual(values = c("Both" = "+","Neither" = "•","Old" = "|","Young" = "—"),
                         name = TeX(r'($\beta$ Signif.)')) +
      scale_size_manual(values = c(4.5,5,2.7,2),name = TeX(r'($\beta$ Signif.)')) +
      scale_color_manual(values = c("gray50","red"),name = TeX(r'($\Delta\beta$ Signif.)'),breaks = c("Not_significant","Significant") ,labels = c("False","True")) + 
      geom_hline(yintercept = 0,color = "black") +
      geom_vline(xintercept = 0,color = "black")
  }
  plot(p_signif)
  return(p_signif)
}
