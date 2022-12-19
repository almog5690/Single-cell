library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)

source("scRNA_seq_utilities.R")

# Add saving figures to file in code!!! (use ggsave !!! )

draw_mean_figures <- function(DF_cor, fig.num) {
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
  }  
  
  #### Figure 2
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
  }  
}
