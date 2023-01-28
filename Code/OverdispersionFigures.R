library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)

source("scRNA_seq_utilities.R")


draw_var_figures <- function(DF_core, fig.num) {
  if(fig.num ==3)
  {
    #### Figure 3
    ## facs
    DF_noise = DF_noise[order(DF_noise$Old_sign - DF_noise$Young_sign,decreasing = T),]  # sorting the data frame by "Old_sign" - "Young_sign"
    DF_noise$ct = interaction(DF_noise$Organs,DF_noise$Cell_type) # creating cell type variable
    DF_noise$ct = factor(DF_noise$ct,levels = DF_noise$ct)
    # creating data frame for the figure
    df_3 = data.frame("Organs" = rep(DF_noise$Organs,2),"Cell_type" = rep(DF_noise$Cell_type,2),"num_gene" = c(DF_noise$Old_sign,DF_noise$Young_sign),
                      "Age" = rep(c("Old","Young"),each = nrow(DF_noise)),"ct" = c(DF_noise$ct,DF_noise$ct))
    # facs figure 3 plot
    p_noise = ggplot(df_3,aes(x = ct,y = num_gene,fill = Age)) + 
      geom_bar(stat = "identity",position = position_dodge()) + 
      theme_classic() + 
      theme(axis.text.x = element_blank(),axis.ticks = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) +
      labs(y = "# significant genes",title = "facs",x = "Cell type") + 
      scale_fill_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) 
    
    N = sum(DF_noise$Old_sign > DF_noise$Young_sign) # Number of cell types with "Old_sign" > "Young_sign"
    p_val_noise = 1-pbinom(N,nrow(DF_noise),0.5) # One sided binomial test P value.
    
    ## droplet
    DF_noise_drop = DF_noise_drop[order(DF_noise_drop$Old_sign - DF_noise_drop$Young_sign,decreasing = T),] # sorting the data frame by "Old_sign" - "Young_sign"
    DF_noise_drop$ct = interaction(DF_noise_drop$Organs,DF_noise_drop$Cell_type) # creating cell type variable
    DF_noise_drop$ct = factor(DF_noise_drop$ct,levels = DF_noise_drop$ct)
    
    # creating data frame for the figure
    df_3_drop = data.frame("Organs" = rep(DF_noise_drop$Organs,2),"Cell_type" = rep(DF_noise_drop$Cell_type,2),"num_gene" = c(DF_noise_drop$Old_sign,DF_noise_drop$Young_sign),
                           "Age" = rep(c("Old","Young"),each = nrow(DF_noise_drop)),"ct" = c(DF_noise_drop$ct,DF_noise_drop$ct))
    # droplet figure 3 plot
    p_noise_drop = ggplot(df_3_drop,aes(x = ct,y = num_gene,fill = Age)) + 
      geom_bar(stat = "identity",position = position_dodge()) + 
      theme_classic() + 
      theme(axis.text.x = element_blank(),axis.ticks = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) +
      labs(y = "# significant genes",title = "droplet",x = "Cell type") + 
      scale_fill_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) 
    
    N = sum(DF_noise_drop$Old_sign > DF_noise_drop$Young_sign) # Number of cell types with "Old_sign" > "Young_sign"
    p_val_noise = 1-pbinom(N,nrow(DF_noise_drop),0.5) # One sided binomial test P value.
    
    # merging facs and droplet plots into figure 3
    title = ggdraw() + draw_label("Differential Over-dispersion tests") # plot title
    p = plot_grid(p_noise,p_noise_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  }  
  
  if(fig.num ==4)
  {
    #### Figure 4
    ## facs
    # Choosing Marrow's late pro-B cell cell type to highlight
    highlight_cell = which(DF_cor_len_disp$Organs == "Marrow" & DF_cor_len_disp$Cell_type == "late pro-B cell")
    DF_cor_len_disp$highlight <- ifelse((1:nrow(DF_cor_len_disp)) == highlight_cell, "highlight", "normal")
    textdf <- DF_cor_len_disp[highlight_cell, ]
    mycolours <- c("highlight" = "green", "normal" = "blue")
    
    # facs over-dispersion vs length correlation coefficients for both old and young
    p_len_disp = ggplot(DF_cor_len_disp,aes(len_disp_young_cor_spearman,len_disp_old_cor_spearman)) + 
      geom_point(color = "blue")+
      geom_abline(slope = 1,intercept = 0,col = "red") +
      labs(title = "facs",x = "Young cor",y = "Old cor") +
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),legend.position = "none") +
      geom_text(data = textdf, aes(x = len_disp_young_spearman * 1.01, y = len_disp_old_cor_spearman * 1.13, label = "B cell")) +
      geom_point(data = textdf,aes(x = len_disp_young_spearman , y = len_disp_old_cor_spearman),color = "orangered",fill = "orange",size = 3) + 
      geom_hline(yintercept = 0,lwd = 0.1) + 
      geom_vline(xintercept = 0,lwd = 0.1)
    
    
    ## facs Marrow's late pro-B cell cell type over-dispersion rank vs selection rank plot
    i = which(organs == "Marrow") # Marrow index
    cell_types_categories = meta.data[[i]]$cell_ontology_class # the names of the different cell-types
    k = which(cell_types_categories == "late pro-B cell") # precursor B cell index
    
    # Loading Marrow's late pro-B cell BASiCS differential over-dispersion test file
    load(file = paste("/tmp/DVT/DVT",organs[i],cell_types_categories[k],"same-mean.RData"))
    
    df_od = test@Results$Disp@Table # Differential over-dispersion test results table
    df_od = df_od[df_od$ResultDiffDisp != "ExcludedFromTesting",] # Filtering genes significant difference in mean expression
    
    Disp = df_od$DispOverall # Overall over-dispersion (old and young together)
    names(Disp) = toupper(df_od$GeneName)
    
    # corr per age group
    Disp_old = df_od$Disp1 # old over-dispersion vector
    Disp_young = df_od$Disp2 # young over-dispersion vector
    names(Disp_old) = names(Disp_young) = names(Disp)
    
    g = na.omit(glength_2[names(Disp)]) # removing na
    len_rank = rank(g,ties.method = "average") # length rank
    Disp_young_rank = rank(Disp_young[names(g)],ties.method = "average",na.last = NA) # young over-dispersion rank
    Disp_old_rank = rank(Disp_old[names(g)],ties.method = "average",na.last = NA) # old over-dispersion rank
    
    # Over-dispersion and length ran data frame for the plot
    df_4 = data.frame("disp" = c(Disp_old_rank,Disp_young_rank),"Length" = c(len_rank,len_rank),"Age" = rep(c("Old","Young"),each = length(len_rank)))
    # old and young plot titles
    age_name = c("Old" = paste0("Old: ","\u03c1","=0.26,p<2.2e-16"),"Young" = paste0("Young: ","\u03c1","=0.12 ,p=1.4e-7"))
    # Marrow's late pro-B cell cell type over-dispersion rank vs length rank plot
    df_4$density =  get_density(df_4$disp, df_4$Length, n = 100)
    p_denst_len = ggplot(df_4) +
      geom_point(aes(x = disp, y = Length, fill = density),color = "white", alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
      scale_fill_gradientn(colors = matlab.like(100)) + 
      facet_wrap(~Age,labeller = labeller(Age = age_name)) +
      geom_smooth(method = "lm",se = F,data = df_4,aes(x = disp, y = Length,color = Age)) +
      scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
      labs(title = paste(organs[i],cell_types_categories[k],sep = ": "),x = "Over-dispersion rank",y = "Length rank") + 
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold")) 
    
    ## droplet
    # Choosing Marrow's late pro-B cell cell type to highlight
    highlight_cell = which(DF_cor_len_disp_drop$Organs == "Marrow" & DF_cor_len_disp_drop$Cell_type == "late pro-B cell")
    DF_cor_len_disp_drop$highlight <- ifelse((1:nrow(DF_cor_len_disp_drop)) == highlight_cell, "highlight", "normal")
    textdf <- DF_cor_len_disp_drop[highlight_cell, ]
    mycolours <- c("highlight" = "green", "normal" = "blue")
    
    p_len_disp_drop = ggplot(DF_cor_len_disp_drop,aes(len_disp_young_cor_spearman,len_disp_old_cor_spearman)) + 
      geom_point(color = "blue")+
      geom_abline(slope = 1,intercept = 0,col = "red") +
      labs(title = "droplet",x = "Young cor",y = "Old cor") +
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),legend.position = "none") + 
      geom_text(data = textdf, aes(x = len_disp_young_spearman * 1.01, y = len_disp_old_cor_spearman * 1.08, label = "B cell")) +
      geom_point(data = textdf,aes(x = len_disp_young_spearman , y = len_disp_old_cor_spearman),color = "orangered",fill = "orange",size = 3) + 
      geom_hline(yintercept = 0,lwd = 0.1) + 
      geom_vline(xintercept = 0,lwd = 0.1)
    
    ## droplet Marrow's late pro-B cell cell type over-dispersion rank vs selection rank plot
    i = which(drop_organs == "Marrow") # Marrow index
    cell_types_categories = meta.data.drop[[i]]$cell_ontology_class # the names of the different cell-types
    k = which(cell_types_categories == "late pro-B cell") # precursor B cell index
    
    # Loading Marrow's precursor B cell BASiCS differential over-dispersion test file
    load(file = paste("/tmp/DVT/DVT",drop_organs[i],cell_types_categories[k],"drop 3-24 same-mean.RData"))
    
    df_od = test@Results$Disp@Table # Differential over-dispersion test results table
    df_od = df_od[df_od$ResultDiffDisp != "ExcludedFromTesting",] # Filtering genes significant difference in mean expression
    
    Disp = df_od$DispOverall # Overall over-dispersion (old and young together)
    names(Disp) = toupper(df_od$GeneName)
    
    # corr per age group
    Disp_old = df_od$Disp1 # old over-dispersion vector
    Disp_young = df_od$Disp2 # young over-dispersion vector
    names(Disp_old) = names(Disp_young) = names(Disp)
    
    g = na.omit(glength_2[names(Disp)]) # removing na
    len_rank = rank(g,ties.method = "average") # length rank
    Disp_young_rank = rank(Disp_young[names(g)],ties.method = "average",na.last = NA) # young over-dispersion rank
    Disp_old_rank = rank(Disp_old[names(g)],ties.method = "average",na.last = NA) # old over-dispersion rank
    
    # Over-dispersion and length ran data frame for the plot
    df_4 = data.frame("disp" = c(Disp_old_rank,Disp_young_rank),"Length" = c(len_rank,len_rank),"Age" = rep(c("Old","Young"),each = length(len_rank)))
    # old and young plot titles
    age_name = c("Old" = paste0("Old: ","\u03c1","=0.26,p<2.2e-16"),"Young" = paste0("Young: ","\u03c1","=0.12 ,p=1.4e-7"))
    # Marrow's late pro-B cell cell type over-dispersion rank vs length rank plot
    df_4$density =  get_density(df_4$disp, df_4$Length, n = 100)
    p_denst_len_drop = ggplot(df_4) +
      geom_point(aes(x = disp, y = Length, fill = density),color = "white", alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
      scale_fill_gradientn(colors = matlab.like(100)) + 
      facet_wrap(~Age,labeller = labeller(Age = age_name)) +
      geom_smooth(method = "lm",se = F,data = df_4,aes(x = disp, y = Length,color = Age)) +
      scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
      labs(title = paste(drop_organs[i],cell_types_categories[k],sep = ": "),x = "Over-dispersion rank",y = "Length rank") + 
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold")) 
    
    # merging facs and droplet plots into figure 4
    title = ggdraw() + draw_label("Genes Length and over-dispersion correlation") # plot title
    p = plot_grid(p_denst_len,p_len_disp,p_denst_len_drop,p_len_disp_drop,labels = LETTERS[1:4])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  } # end Figure 4
  
  if(fig.num ==5)
  {
    #### Figure 5
    ## facs
    # Choosing Marrow's precursor B cell cell type to highlight
    highlight_cell = which(DF_cor_selc_disp$Organs == "Marrow" & DF_cor_selc_disp$Cell_type == "precursor B cell")
    DF_cor_selc_disp$highlight <- ifelse((1:nrow(DF_cor_selc_disp)) == highlight_cell, "highlight", "normal")
    textdf <- DF_cor_selc_disp[highlight_cell, ]
    mycolours <- c("highlight" = "green", "normal" = "blue")
    
    # facs over-dispersion vs selection correlation coefficients for both old and young
    p_selc_disp = ggplot(DF_cor_selc_disp,aes(selc_disp_young_cor_spearman,selc_disp_old_cor_spearman)) + 
      geom_point(color = "blue")+
      geom_abline(slope = 1,intercept = 0,col = "red") +
      labs(title = "facs",x = "Young cor",y = "Old cor") +
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),legend.position = "none") + 
      geom_text(data = textdf, aes(x = selc_disp_young_spearman * 1.11, y = selc_disp_old_cor_spearman * 1.0, label = "B cell")) +
      geom_point(data = textdf,aes(x = selc_disp_young_spearman , y = selc_disp_old_cor_spearman),color = "orangered",fill = "orange",size = 3) + 
      geom_hline(yintercept = 0,lwd = 0.1) + 
      geom_vline(xintercept = 0,lwd = 0.1)
    
    ## facs Marrow's precursor B cell cell type over-dispersion rank vs selection rank plot
    i = which(organs == "Marrow") # Marrow index
    cell_types_categories = meta.data[[i]]$cell_ontology_class # the names of the different cell-types
    k = which(cell_types_categories == "precursor B cell") # precursor B cell index
    
    # Loading Marrow's precursor B cell BASiCS differential over-dispersion test file
    load(file = paste("/tmp/DVT/DVT",organs[i],cell_types_categories[k],"same-mean.RData"))
    
    df_od = test@Results$Disp@Table # Differential over-dispersion test results table
    df_od = df_od[df_od$ResultDiffDisp != "ExcludedFromTesting",] # Filtering genes significant difference in mean expression
    
    Disp = df_od$DispOverall # Overall over-dispersion (old and young together)
    names(Disp) = toupper(df_od$GeneName)
    
    # corr per age group
    Disp_old = df_od$Disp1 # old over-dispersion vector
    Disp_young = df_od$Disp2 # young over-dispersion vector
    names(Disp_old) = names(Disp_young) = names(Disp)
    
    g = na.omit(gene_selection[names(Disp)]) # removing na
    
    selc_rank = rank(g,ties.method = "average",na.last = NA) # selection rank
    Disp_young_rank = rank(Disp_young[names(g)],ties.method = "average",na.last = NA) # young over-dispersion rank
    Disp_old_rank = rank(Disp_old[names(g)],ties.method = "average",na.last = NA) # old over-dispersion rank
    
    # Over-dispersion and selection ran data frame for the plot
    df_4 = data.frame("disp" = c(Disp_old_rank,Disp_young_rank),"Selection" = c(selc_rank,selc_rank),"Age" = rep(c("Old","Young"),each = length(selc_rank)))
    # old and young plot titles
    age_name = c("Old" = paste0("Old: ","\u03c1","=-0.18,p<2.2e-16"),"Young" = paste0("Young: ","\u03c1","=-0.25 ,p<2.2e-16"))
    
    # Marrow's precursor B cell cell type over-dispersion rank vs selection rank plot
    df_4$density =  get_density(df_4$disp, df_4$Selection, n = 100)
    p_denst = ggplot(df_4) +
      geom_point(aes(x = disp, y = Selection, fill = density),color = "white", alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
      scale_fill_gradientn(colors = matlab.like(100)) + 
      facet_wrap(~Age,labeller = labeller(Age = age_name)) + 
      geom_smooth(method = "lm",se = F,data = df_4,aes(x = disp, y = Selection,color = Age)) +
      scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
      labs(title = paste(organs[i],cell_types_categories[k],sep = ": "),x = "Over-dispersion rank",y = "Selection rank") + 
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold"),axis.text.x = element_text(size = 8)) 
    
    ## droplet
    # Choosing Marrow's precursor B cell cell type to highlight
    highlight_cell = which(DF_cor_selc_disp_drop$Organs == "Marrow" & DF_cor_selc_disp_drop$Cell_type == "precursor B cell")
    DF_cor_selc_disp_drop$highlight <- ifelse((1:nrow(DF_cor_selc_disp_drop)) == highlight_cell, "highlight", "normal")
    textdf <- DF_cor_selc_disp_drop[highlight_cell, ]
    mycolours <- c("highlight" = "green", "normal" = "blue")
    
    # droplet over-dispersion vs selection correlation coefficients for both old and young
    p_selc_disp_drop = ggplot(DF_cor_selc_disp_drop,aes(selc_disp_young_cor_spearman,selc_disp_old_cor_spearman)) + 
      geom_point(color = "blue")+
      geom_abline(slope = 1,intercept = 0,col = "red") +
      labs(title = "droplet",x = "Young cor",y = "Old cor") +
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),legend.position = "none") + 
      geom_text(data = textdf, aes(x = selc_disp_young_spearman * 1.17, y = selc_disp_old_cor_spearman * 1.0, label = "B cell")) +
      geom_point(data = textdf,aes(x = selc_disp_young_spearman , y = selc_disp_old_cor_spearman),color = "orangered",fill = "orange",size = 3) + 
      geom_hline(yintercept = 0,lwd = 0.1) + 
      geom_vline(xintercept = 0,lwd = 0.1)
    
    ## droplet Marrow's precursor B cell cell type over-dispersion rank vs selection rank plot
    i = which(drop_organs == "Marrow") # Marrow index
    cell_types_categories = meta.data.drop[[i]]$cell_ontology_class # the names of the different cell-types
    k = which(cell_types_categories == "precursor B cell") # precursor B cell index
    
    # Loading Marrow's precursor B cell BASiCS differential over-dispersion test file
    load(file = paste("/tmp/DVT/DVT",drop_organs[i],cell_types_categories[k],"drop 3-24 same-mean.RData"))
    
    df_od = test@Results$Disp@Table # Differential over-dispersion test results table
    df_od = df_od[df_od$ResultDiffDisp != "ExcludedFromTesting",] # Filtering genes significant difference in mean expression
    
    Disp = df_od$DispOverall # Overall over-dispersion (old and young together)
    names(Disp) = toupper(df_od$GeneName)
    
    # corr per age group
    Disp_old = df_od$Disp1 # old over-dispersion vector
    Disp_young = df_od$Disp2 # young over-dispersion vector
    names(Disp_old) = names(Disp_young) = names(Disp)
    
    g = na.omit(gene_selection[names(Disp)]) # removing na
    
    selc_rank = rank(g,ties.method = "average",na.last = NA) # selection rank
    Disp_young_rank = rank(Disp_young[names(g)],ties.method = "average",na.last = NA) # young over-dispersion rank
    Disp_old_rank = rank(Disp_old[names(g)],ties.method = "average",na.last = NA) # old over-dispersion rank
    
    # Over-dispersion and selection ran data frame for the plot
    df_4 = data.frame("disp" = c(Disp_old_rank,Disp_young_rank),"Selection" = c(selc_rank,selc_rank),"Age" = rep(c("Old","Young"),each = length(selc_rank)))
    # old and young plot titles
    age_name = c("Old" = paste0("Old: ","\u03c1","=-0.18,p<2.2e-16"),"Young" = paste0("Young: ","\u03c1","=-0.25 ,p<2.2e-16"))
    
    # Marrow's precursor B cell cell type over-dispersion rank vs selection rank plot
    df_4$density =  get_density(df_4$disp, df_4$Selection, n = 100)
    p_denst_drop = ggplot(df_4) +
      geom_point(aes(x = disp, y = Selection, fill = density),color = "white", alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
      scale_fill_gradientn(colors = matlab.like(100)) + 
      facet_wrap(~Age,labeller = labeller(Age = age_name)) + 
      geom_smooth(method = "lm",se = F,data = df_4,aes(x = disp, y = Selection,color = Age)) +
      scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
      labs(title = paste(drop_organs[i],cell_types_categories[k],sep = ": "),x = "Over-dispersion rank",y = "Selection rank") + 
      theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold"),axis.text.x = element_text(size = 8)) 
    
    # merging facs and droplet plots into figure 5
    title = ggdraw() + draw_label("Genes selection and over-dispersion correlation") # plot title
    p = plot_grid(p_denst,p_selc_disp,p_denst_drop,p_selc_disp_drop,labels = LETTERS[1:4])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  } # End figure 5
  
  if(fig.num ==8)
  {
    #### Figure 8
    ## facs
    p_log_OD_len_sm = ggplot(len_OD_reg_data_sm,aes(len_log_OD_young_reg,len_log_OD_old_reg)) + 
      geom_point(color = "blue") + 
      geom_abline(slope = 1,intercept = 0,col = "red") + 
      labs(title = "facs",x = TeX(r'($\beta_{o}^{(t,young)}$)'),y = TeX(r'($\beta_{o}^{(t,old)}$)')) +
      theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8)) + 
      geom_hline(yintercept = 0,lwd = 0.1) + 
      geom_vline(xintercept = 0,lwd = 0.1)
    
    ## droplet
    p_log_OD_len_sm_drop = ggplot(len_OD_reg_data_sm_drop,aes(len_log_OD_young_reg,len_log_OD_old_reg)) + 
      geom_point(color = "blue") + 
      geom_abline(slope = 1,intercept = 0,col = "red") + 
      labs(title = "droplet",x = TeX(r'($\beta_{o}^{(t,young)}$)'),y = TeX(r'($\beta_{o}^{(t,old)}$)')) +
      theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8)) +
      geom_hline(yintercept = 0,lwd = 0.1) + 
      geom_vline(xintercept = 0,lwd = 0.1)
    
    # merging facs and droplet plots into figure 8
    title = ggdraw() + draw_label("Length vs Over-dispersion - regressing mean and selection") # plot title
    p = plot_grid(p_log_OD_len_sm,p_log_OD_len_sm_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  } # End figure 8
  
  if(fig.num ==9)
  {
    #### Figure 9
    ## facs
    p_log_OD_selc_reg_sm = ggplot(selc_OD_reg_data_sm,aes(selc_log_OD_young_reg,selc_log_OD_old_reg)) + 
      geom_point(color = "blue") + 
      geom_abline(slope = 1,intercept = 0,col = "red") + 
      labs(title = "facs",x = TeX(r'($\beta_{o}^{(t,young)}$)'),y = TeX(r'($\beta_{o}^{(t,old)}$)')) +
      theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8)) + 
      geom_hline(yintercept = 0,lwd = 0.1) + 
      geom_vline(xintercept = 0,lwd = 0.1)
    
    ## droplet
    p_log_OD_selc_reg_sm_drop = ggplot(selc_OD_reg_data_sm_drop,aes(selc_log_OD_young_reg,selc_log_OD_old_reg)) + 
      geom_point(color = "blue") + 
      geom_abline(slope = 1,intercept = 0,col = "red") + 
      labs(title = "droplet",x = TeX(r'($\beta_{o}^{(t,young)}$)'),y = TeX(r'($\beta_{o}^{(t,old)}$)')) +
      theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8)) + 
      geom_hline(yintercept = 0,lwd = 0.1) + 
      geom_vline(xintercept = 0,lwd = 0.1)
    
    # merging facs and droplet plots into figure 9
    title = ggdraw() + draw_label("Selection vs Over-dispersion - regressing mean and length") # plot title
    p = plot_grid(p_log_OD_selc_reg_sm,p_log_OD_selc_reg_sm_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  } # End figure 9
  
  
  if(fig.num ==12)
  {
    #### Figure 12
    ## facs
    df_selc_reg_facs$CT = interaction(df_selc_reg_facs$Organs,df_selc_reg_facs$Cell_type,sep = ":") # adding cell type variable
    
    p_selc_disp_scale_reg = ggplot(df_selc_reg_facs,aes(y = reorder(CT,beta_selc_all),x = beta_selc_all,fill = -log10(pval_all))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "facs",x = TeX(r'($\beta^{(t)}_{s}$)'),y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    ## droplet
    df_selc_reg$CT = interaction(df_selc_reg$Organs,df_selc_reg$Cell_type,sep = ":") # adding cell type variable
    
    p_selc_disp_scale_reg_drop = ggplot(df_selc_reg,aes(y = reorder(CT,beta_selc_all),x = beta_selc_all,fill = -log10(pval_all))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "droplet",x = TeX(r'($\beta^{(t)}_{s}$)'),y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    # merging facs and droplet plots into figure 12
    title = ggdraw() + draw_label("Over-dispersion vs selection regressing mean") # plot title
    p = plot_grid(p_selc_disp_scale_reg,p_selc_disp_scale_reg_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  } # End figure 12  
  
  if(fig.num ==13)
  {
    #### Figure 13
    ## facs
    df_fc_reg_facs$CT = interaction(df_fc_reg_facs$Organs,df_fc_reg_facs$Cell_type,sep = ":") # adding cell type variable
    
    p_selc_disp_scale_log_fc_reg = ggplot(df_fc_reg_facs,aes(y = reorder(CT,beta_selc_log_fc),x = beta_selc_log_fc,fill = -log10(pval_log_fc))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "facs",x = TeX(r'($\beta^{(t)}_{s}$)'),y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    ## droplet
    df_fc_reg$CT = interaction(df_fc_reg$Organs,df_fc_reg$Cell_type,sep = ":") # adding cell type variable
    
    p_selc_disp_scale_log_fc_reg_drop = ggplot(df_fc_reg,aes(y = reorder(CT,beta_selc_log_fc),x = beta_selc_log_fc,fill = -log10(pval_log_fc))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "droplet",x = TeX(r'($\beta^{(t)}_{s}$)'),y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    # merging facs and droplet plots into figure 13
    title = ggdraw() + draw_label("Over-dispersion FC vs selection regressing mean") # plot title
    p = plot_grid(p_selc_disp_scale_log_fc_reg,p_selc_disp_scale_log_fc_reg_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  } # end figure 13
  
  if(fig.num ==14)
  {
    #### Figure 14
    ## facs
    df_len_reg_facs$CT = interaction(df_len_reg_facs$Organs,df_len_reg_facs$Cell_type,sep = ":") # adding cell type variable
    
    p_len_disp_scale_reg = ggplot(df_len_reg_facs,aes(y = reorder(CT,beta_len_all),x = beta_len_all,fill = -log10(pval_all))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "facs",x = TeX(r'($\beta^{(t)}_{l}$)'),y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    ## droplet
    df_len_reg$CT = interaction(df_len_reg$Organs,df_len_reg$Cell_type,sep = ":") # adding cell type variable
    
    p_len_disp_scale_reg_drop = ggplot(df_len_reg,aes(y = reorder(CT,beta_len_all),x = beta_len_all,fill = -log10(pval_all))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "droplet",x = TeX(r'($\beta^{(t)}_{l}$)'),y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    # merging facs and droplet plots into figure 14
    title = ggdraw() + draw_label("Over-dispersion vs length regressing means") # plot title
    p = plot_grid(p_len_disp_scale_reg,p_len_disp_scale_reg_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  } # End figure 14
  
  if(fig.num ==15)
  {
    #### Figure 15
    ## facs
    p_len_disp_scale_reg_fc = ggplot(df_len_reg_facs,aes(y = reorder(CT,beta_len_fc),x = beta_len_fc,fill = -log10(pval_fc))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "facs",x = TeX(r'($\beta^{(t)}_{l}$)'),y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    ## droplet
    p_len_disp_scale_reg_fc_drop = ggplot(df_len_reg,aes(y = reorder(CT,beta_len_fc),x = beta_len_fc,fill = -log10(pval_fc))) + 
      geom_bar(stat = "identity") + 
      theme_classic() +
      scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
      labs(title = "droplet",x = TeX(r'($\beta^{(t)}_{l}$)'),y = "Cell types",fill = TeX(r'($-log_{10}(\it{Pval})$)')) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 
    
    # merging facs and droplet plots into figure 15
    title = ggdraw() + draw_label("Over-dispersion FC vs length regressing means") # plot title
    p = plot_grid(p_len_disp_scale_reg_fc,p_len_disp_scale_reg_fc_drop,labels = LETTERS[1:2])
    plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))
  } # End figure 15
  
} # End function 






# Draw figure of bars of correlation or beta regression coefficients with color showing log p-values
draw_var_cor_bars_figure <- function(fig.num, data.types, feature.types, DF_cors, analysis.figures.dir)
{
  n.datas <- length(data.types)
  p_feature_vs_mean_bar <- vector("list", n.datas)
  
  group.str = switch(as.character(fig.num), "1" = "all", "11" = "fc", "111" = "fc_abs")  #  == 1) "all" else "fc"
  for(feature.type in feature.types) # Here plot each feature vs. mean expression, not just selection
  {
    print("Feature:")
    print(feature.type)
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
      print("ggplot:")
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
    ggsave(paste0(analysis.figures.dir, "OD.Figure", fig.num, '.', paste(data.types, collapse = "_"), 
                  '.', feature.type, '.png'), height = 6,width = 9)  # Modify name to get figure  
  }    
}


# Detailed regression of one cell type with two-dim density, and other young vs. old. scatter plots 
# Plot colorful scatter 
draw_var_cor_scatters_figure <- function(fig.num, data.types, feature.types, DF_cors, analysis.figures.dir, num.cell.types,
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
                            ", Sign test: ", num.old.bigger.young, "/", num.cell.types[i], " Pval=", signif(num.old.bigger.young.pval, 3)), 
             x = paste0("Young", fig.str), y = paste0("Old ", fig.str)) +
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
      # counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the Lung tissue
      # SC_gene_name = toupper(rownames(SC)) # genes names in upper case letters
      # if(length(SC_gene_name) != dim(counts.mat)[1]) { # For Rats, names are already in the matrix
      #   SC_gene_name = toupper(rownames(counts.mat)) 
      # } 
      # rownames(counts.mat) = SC_gene_name # make sure upper 

      if(data.types[i] == "CR.Rat"){ # Convert to dummy variables 
        cell_types = as.numeric(SC@meta.data$cell_types) # Cell types vector
      } else
        cell_types = SC@meta.data$cell.ontology.class # Cell types vector
      
      # gene_feat = gene_features[[feature.type]][gene_name %in% (SC_gene_name)] # filtering the selection score to genes that are found in the current tissue
      # cur_gene_name = gene_name[gene_name %in% (SC_gene_name )] # the names of the filtered genes
      
      DVT.file.name = get_DVT_file_name(data.type, samples$organs[tissue.ind],cell_types_categories[k])
      
      # Loading Marrow's late pro-B cell BASiCS differential over-dispersion test file
      load(DVT.file.name)
      
      df_od = test@Results$Disp@Table # Differential over-dispersion test results table
      df_od = df_od[df_od$ResultDiffDisp != "ExcludedFromTesting",] # Filtering genes significant difference in mean expression
      
      Disp = df_od$DispOverall # Overall over-dispersion (old and young together)
      names(Disp) = toupper(df_od$GeneName)
      
      # corr per age group
      Disp_old = df_od$Disp1 # old over-dispersion vector
      Disp_young = df_od$Disp2 # young over-dispersion vector
      names(Disp_old) = names(Disp_young) = names(Disp)
      
      g = na.omit(gene_feat)
      selc_rank = rank(g[names(Disp_old)], ties.method = "average")
      Disp_young_rank = rank(Disp_young,ties.method = "average",na.last = NA) # young over-dispersion rank
      Disp_old_rank = rank(Disp_old,ties.method = "average",na.last = NA) # old over-dispersion rank
      
      
      # data frame contains old and young genes mean expression and selection ranks for the plots
      df_4 = data.frame("Disp" = c(Disp_old_rank, Disp_young_rank),
                        "gene.feature" = c(selc_rank, selc_rank),
                        "Age" = rep(c("Old","Young"),each = length(selc_rank)))
      # Display correlation (no need for p-value)
      age_name = c("Old" = paste0("Old: ","\u03c1","=", round(DF_cors[[i]][highlight_cell,cor.old.col], 3)), # ,p<2.2e-16"),
                   "Young" = paste0("Young: ","\u03c1","=", round(DF_cors[[i]][highlight_cell,cor.young.col], 3))) # ,p<2.2e-16"))
      
      # getting the 2D density of selection and mean for the plots
      df_4$density =  get_density(df_4$Disp, df_4$gene.feature, n = 100) # problem here: df_4 is empty!!!! 
      
      
      # Selection rank vs mean expression rank for both young and old for the Lung Pneumocyte cell type
      p_denst[[i]] = ggplot(df_4) +
        geom_point(aes(x = Disp, y = gene.feature, fill = density),color = "white", 
                   alpha = 1, size = 1.8,  shape = 21, show.legend = T) +
        scale_fill_gradientn(colors = matlab.like(100)) + 
        facet_wrap(~Age, labeller = labeller(Age = age_name)) + 
        geom_smooth(method = "lm",se = F,data = df_4, aes(x = Disp, y = gene.feature, color = Age)) +
        scale_color_manual(values = c("Old"="#00ba38", "Young"="#f8766d")) +
        labs(title = paste(samples$organs[i],cell_types_categories[k], sep = ": "), x = "Over-dispersion expression rank", y = paste0(feature.type, " rank")) + 
        theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold")) +
        guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
    } # end loop on data types 
    
    # merging all datasets plots into figure 2
    multi.plot <- ggarrange(plotlist = c(rbind(p_denst, p3)), nrow = 2, ncol = n.datas) # all in the same row 
    annotate_figure(multi.plot, top = text_grob(paste0("Genes ", feature.type, " and Over-dispersion correlation"), 
                                                color = "black", face = "bold", size = 14))
    ggsave(paste0(analysis.figures.dir, "OD.Figure", fig.num, '.', paste(data.types, collapse = "_"),
                  '.',  feature.type, '.png'), height = 6,width = 9)  # Modify name to get figure  
  }  # loop on explanatory features 
} # End function for plotting fig. 2,22

