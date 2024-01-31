# Temp code for debugging

union_names = union(names(gene_feat), names(bad_gene_feat))
union_gene_feat = rep(0, length(union_names))
names(union_gene_feat) = union_names
union_gene_feat[names(gene_feat)] = gene_feat
union_gene_feat[names(bad_gene_feat)] = bad_gene_feat
union_mean_expr_young = rep(0, length(union_names))
names(union_mean_expr_young) = union_names
union_mean_expr_young[names(gene_mean_young)] = gene_mean_young
union_mean_expr_young[names(bad_expr_mean_young)] = bad_expr_mean_young
cor(union_mean_expr_young, union_gene_feat, use = "complete.obs", method = "spearman")


I = intersect(names(bad_gene_feat), names(gene_feat))
unique(gene_feat[I]- bad_gene_feat[I])
plot(gene_feat, gene_mean_young, pch=20, cex=0.1)
points(bad_gene_feat, bad_expr_mean_young, pch=1, cex=0.1, col="red")
points(bad_gene_feat[I], bad_expr_mean_young[I], pch=1, cex=0.1, col="green")
points(gene_feat[I], gene_mean_young[I], pch=1, cex=0.1, col="blue")


# Plot ranks (more spread)
plot(rank(gene_feat)/length(gene_feat), 
     rank(gene_mean_young)/length(gene_feat), pch=20, cex=0.1)
plot(rank(bad_gene_feat)/length(bad_gene_feat), 
     rank(bad_expr_mean_young)/length(bad_gene_feat), pch=1, cex=0.1, col="red")
plot(rank(bad_gene_feat[I])/length(I), 
     rank(bad_expr_mean_young[I])/length(I), pch=1, cex=0.1, col="green")
#points(gene_feat[I], gene_mean_young[I], pch=1, cex=0.1, col="blue")


n_feat = sum(!is.na(gene_feat))
plot((1:n_feat)/n_feat, sort(gene_feat), pch=20, cex=0.1)
n_bad_feat = sum(!is.na(bad_gene_feat))
points((1:n_bad_feat)/n_bad_feat, sort(bad_gene_feat), pch=20, cex=0.1, col="red")



n_mean = sum(!is.na(gene_mean_young))
plot((1:n_mean)/n_mean, sort(gene_mean_young), pch=20, cex=0.1)
n_bad_mean = sum(!is.na(bad_expr_mean_young))
points((1:n_bad_mean)/n_bad_mean, sort(bad_expr_mean_young), pch=20, cex=0.1, col="red")


df_4 = data.frame("Mean" = c(mean_old_rank, mean_young_rank), # take ranks
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
  labs(title = paste(samples$organs[i],cell_types_categories[k], sep = ": "), x = "Mean expression rank", y = paste0(feature.type, " rank")) + 
  theme(plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8),axis.title = element_text(size = 10),strip.text.x = element_text(size = 8,face = "bold")) +
  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
