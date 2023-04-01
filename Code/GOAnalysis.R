


genes_for_GO <- function(OD_genes_data,age = "Old")
{
  n = (ncol(OD_genes_data)-1)/3
  age_res = ifelse(age == "Old","Old+","Young+")
  n_cell_signif = rowSums(OD_genes_data == age_res,na.rm = T)
  n_cell_exclude = rowSums(OD_genes_data == "ExcludedFromTesting",na.rm = T)
  
  Y = n_cell_signif/(n-n_cell_exclude)
  P = mean(Y)
  Z = (Y-P)*sqrt((n-n_cell_exclude)/(P*(1-P)))
  names(Z) = OD_genes_data[,1]
  
  genes_for_GO = names(Z)[Z>quantile(Z,0.95)]
  
  return(genes_for_GO)
}

facs_data = as.data.frame(read_excel("Over-dispersion genes data.xlsx"))

facs_old_genes_for_GO = genes_for_GO(facs_data)

drop_data = as.data.frame(read_excel("Over-dispersion genes data.xlsx",sheet = 2))

drop_old_genes_for_GO = genes_for_GO(drop_data)

write.xlsx(data.frame("facs" = facs_old_genes_for_GO,droplet = c(drop_old_genes_for_GO,rep(NA,length(facs_old_genes_for_GO) - length(drop_old_genes_for_GO)))),
           file = "Old genes for GO.xlsx")

facs_young_genes_for_GO = genes_for_GO(facs_data,age = "Young")

drop_young_genes_for_GO = genes_for_GO(drop_data,age = "Young")

write.xlsx(data.frame("facs" = facs_young_genes_for_GO,droplet = c(drop_young_genes_for_GO,rep(NA,length(facs_young_genes_for_GO) - length(drop_young_genes_for_GO)))),
           file = "Young genes for GO.xlsx")
