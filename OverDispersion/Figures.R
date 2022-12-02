library(ggplot2)
library(cowplot)
library(latex2exp)

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


#### Figure 4

#### Figure 5

#### Figure 8
## facs
p_OD_len_sm = ggplot(len_OD_reg_data_sm,aes(len_OD_young_reg,len_OD_old_reg)) + 
  geom_point(color = "blue") + 
  geom_abline(slope = 1,intercept = 0,col = "red") + 
  labs(title = "facs",x = "\u03B2 young",y = "\u03B2 old") +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8))

## droplet
p_OD_len_sm_drop = ggplot(len_OD_reg_data_sm_drop,aes(len_OD_young_reg,len_OD_old_reg)) + 
  geom_point(color = "blue") + 
  geom_abline(slope = 1,intercept = 0,col = "red") + 
  labs(title = "facs",x = "\u03B2 young",y = "\u03B2 old") +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8))

# merging facs and droplet plots into figure 8
title = ggdraw() + draw_label("Length vs Over-dispersion - regressing mean and selection") # plot title
p = plot_grid(p_OD_len_sm,p_OD_len_sm_drop,labels = LETTERS[1:2])
plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))

#### Figure 9
## facs
p_OD_selc_reg_sm = ggplot(selc_OD_reg_data_sm,aes(selc_OD_young_reg,selc_OD_old_reg)) + 
  geom_point(color = "blue") + 
  geom_abline(slope = 1,intercept = 0,col = "red") + 
  labs(title = "facs",x = "\u03B2 young",y = "\u03B2 old") +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8))

## droplet
p_OD_selc_reg_sm_drop = ggplot(selc_OD_reg_data_sm_drop,aes(selc_OD_young_reg,selc_OD_old_reg)) + 
  geom_point(color = "blue") + 
  geom_abline(slope = 1,intercept = 0,col = "red") + 
  labs(title = "facs",x = "\u03B2 young",y = "\u03B2 old") +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8))

# merging facs and droplet plots into figure 9
title = ggdraw() + draw_label("Selection vs Over-dispersion - regressing mean and length") # plot title
p = plot_grid(p_OD_selc_reg_sm,p_OD_selc_reg_sm_drop,labels = LETTERS[1:2])
plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))


#### Figure 12
## facs
df_selc_reg_facs$CT = interaction(df_selc_reg_facs$Organs,df_selc_reg_facs$Cell_type,sep = ":") # adding cell type variable

p_selc_disp_scale_reg = ggplot(df_selc_reg_facs,aes(y = reorder(CT,beta_selc_all),x = beta_selc_all,fill = -log10(pval_all))) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
  labs(title = "facs",x = "\u03B2 selection",y = "Cell types",fill = "-log10(P_val)") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 

## droplet
df_selc_reg$CT = interaction(df_selc_reg$Organs,df_selc_reg$Cell_type,sep = ":") # adding cell type variable

p_selc_disp_scale_reg_drop = ggplot(df_selc_reg,aes(y = reorder(CT,beta_selc_all),x = beta_selc_all,fill = -log10(pval_all))) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
  labs(title = "droplet",x = "\u03B2 selection",y = "Cell types",fill = "-log10(P_val)") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 

# merging facs and droplet plots into figure 12
title = ggdraw() + draw_label("Over-dispersion vs selection regressing mean") # plot title
p = plot_grid(p_selc_disp_scale_reg,p_selc_disp_scale_reg_drop,labels = LETTERS[1:2])
plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))


#### Figure 13
## facs
df_fc_reg_facs$CT = interaction(df_fc_reg_facs$Organs,df_fc_reg_facs$Cell_type,sep = ":") # adding cell type variable

p_selc_disp_scale_log_fc_reg = ggplot(df_fc_reg_facs,aes(y = reorder(CT,beta_selc_log_fc),x = beta_selc_log_fc,fill = -log10(pval_log_fc))) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
  labs(title = "facs",x = "\u03B2 selection",y = "Cell types",fill = "-log10(P_val)") +
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


#### Figure 14
## facs
df_len_reg_facs$CT = interaction(df_len_reg_facs$Organs,df_len_reg_facs$Cell_type,sep = ":") # adding cell type variable

p_len_disp_scale_reg = ggplot(df_len_reg_facs,aes(y = reorder(CT,beta_len_all),x = beta_len_all,fill = -log10(pval_all))) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
  labs(title = "facs",x = "\u03B2 length",y = "Cell types",fill = "-log10(P_val)") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 

## droplet
df_len_reg$CT = interaction(df_len_reg$Organs,df_len_reg$Cell_type,sep = ":") # adding cell type variable

p_len_disp_scale_reg_drop = ggplot(df_len_reg,aes(y = reorder(CT,beta_len_all),x = beta_len_all,fill = -log10(pval_all))) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
  labs(title = "droplet",x = "\u03B2 length",y = "Cell types",fill = "-log10(P_val)") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 

# merging facs and droplet plots into figure 14
title = ggdraw() + draw_label("Over-dispersion vs length regressing means") # plot title
p = plot_grid(p_len_disp_scale_reg,p_len_disp_scale_reg_drop,labels = LETTERS[1:2])
plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))

#### Figure 15
## facs
p_len_disp_scale_reg_fc = ggplot(df_len_reg_facs,aes(y = reorder(CT,beta_len_fc),x = beta_len_fc,fill = -log10(pval_fc))) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
  labs(title = "facs",x = "\u03B2 length",y = "Cell types",fill = "-log10(P_val)") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 

## droplet
p_len_disp_scale_reg_fc_drop = ggplot(df_len_reg,aes(y = reorder(CT,beta_len_fc),x = beta_len_fc,fill = -log10(pval_fc))) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  scale_fill_gradient2(high = "red",mid = "white",low = "blue",midpoint = -log10(0.05)) +
  labs(title = "droplet",x = "\u03B2 length",y = "Cell types",fill = "-log10(P_val)") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10)) 

# merging facs and droplet plots into figure 15
title = ggdraw() + draw_label("Over-dispersion FC vs length regressing means") # plot title
p = plot_grid(p_len_disp_scale_reg_fc,p_len_disp_scale_reg_fc_drop,labels = LETTERS[1:2])
plot_grid(title,p,ncol=1, rel_heights=c(0.1, 1))