library(Seurat)

droplet.files = list.files(pattern = "drop.r") # list of tissue from droplet technology
#removing large intestine and trachea because they have only one age group and also fat,skin and pancreas.
droplet.files = droplet.files[-grep(c("Fat|Large_Intestine|Skin|Pancreas|Trachea"),droplet.files)]
drop_organs = sapply(droplet.files,function(f) unlist(strsplit(f,"[.]"))[1]) 
drop_organs = unname(drop_organs) # droplet organs vector

# loading the meta data names for the droplet technology
load("meta.data.drop.Rdata")



# reading the selection score data
select = read.delim("gnomad.v2.1.1.lof_metrics.by_gene.txt")
gene_name = toupper(select$gene) # making all the gene names in to upper case letters
gene_selection = select$pLI # the selection score vector
names(gene_selection) = gene_name # naming each score with the gene he belongs to


young_ages = "3m" # The ages of the young mice 
old_ages_1 = c("21m","24m") # The ages of the old mice
old_ages_2 = c("18m","24m") # secondery old mice age incase no mice in previous "old_ages"

DF_cor_drop = c()
for(i in 1:length(drop_organs)){
  SC = readRDS(file = paste(drop_organs[i],"drop","rds",sep = ".")) # current tissue seurat object
  counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the current tissue
  young.ind = c(SC@meta.data$age %in% young_ages) # index for cells that came from 3 month old mouses
  old.ind = c(SC@meta.data$age %in% old_ages_1) # index for cells that came from old mouses
  if(sum(old.ind) == 0){
    old.ind = c(SC@meta.data$age %in% old_ages_2)
  }
  if(i==6){
    old.ind = SC@meta.data$age %in% c("18m","21m")
  }
  
  n_cell = SC@assays$RNA@counts@Dim[2] # Number of cells
  n_genes = SC@assays$RNA@counts@Dim[1] # Number of genes
  SC_gene_name = toupper(rownames(SC)) # Gene names (in uppercase)
  rownames(counts.mat) = SC_gene_name 
  
  cell_types = SC@meta.data$cell.ontology.class # Cell types vector
  cell_types_categories = meta.data.drop[[i]]$cell_ontology_class # Cell type names
  n_cell_types = max(cell_types) # Number of cell types
  earase = vector()
  
  for(ct_ind in 0:n_cell_types){
    #filtering cell types with less then 100 cells or less the 20 cells in each the age groups
    if(sum(cell_types==ct_ind)<100|sum(cell_types==ct_ind & young.ind)<20|sum(cell_types==ct_ind & old.ind)<20){
      earase = c(earase,ct_ind+1)
      next()
    }
  }
  
  gene_selc = gene_selection[gene_name %in% (SC_gene_name)] # filtering the selection score to genes that are found in the current tissue
  cur_gene_name = gene_name[gene_name %in% (SC_gene_name)] # the names of the filtered genes
  
  cells_ind = c(1:(n_cell_types+1))[-earase]
  if(length(cells_ind) == 0) cells_ind = (1:(n_cell_types+1))
  
  for(k in cells_ind){
    gene_selc = gene_selection[gene_name %in% (SC_gene_name)] # filtering the selection score to genes that are found in the current tissue
    
    # genes mean vector
    gene_mean = rowMeans(counts.mat[,cell_types==k-1])
    gene_mean = gene_mean[cur_gene_name]
    
    # Mean expression and selection correlation (old and young together)
    mean_selc_cor = cor(gene_mean,gene_selc,use = "complete.obs",method = "spearman")
    pval_all = cor.test(gene_mean,gene_selc,use = "complete.obs",method = "spearman")$p.value
    
    # genes mean vector-old
    gene_mean_old = rowMeans(counts.mat[,(cell_types==k-1)&(old.ind)])
    gene_mean_old = gene_mean_old[cur_gene_name]
    
    # genes mean vector-young
    gene_mean_young = rowMeans(counts.mat[,(cell_types==k-1)&(young.ind)])
    gene_mean_young = gene_mean_young[cur_gene_name]
    
    # Filtering genes with low expression for old (less than 10 counts)
    genes_ind = rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&(!young.ind)]) > 10
    names(genes_ind) = toupper(names(genes_ind))
    genes_ind_old = genes_ind[cur_gene_name]
    
    # Filtering genes with low expression for young (less than 10 counts)
    genes_ind = rowSums(SC@assays$RNA@counts[,(cell_types==k-1)&(young.ind)]) > 10
    names(genes_ind) = toupper(names(genes_ind))
    genes_ind_young = genes_ind[cur_gene_name]
    
    mean_young = gene_mean_young[genes_ind_young] # filtered young mean expression vector
    mean_old = gene_mean_old[genes_ind_young]# filtered young mean expression vector
    
    # Mean expression fold-change and selection correlation
    fc_cor = cor((mean_old/mean_young),gene_selc[genes_ind_young],use = "complete.obs",method = "spearman")
    pval_fc = cor.test((mean_old/mean_young),gene_selc[genes_ind_young],use = "complete.obs",method = "spearman")$p.value
    
    gene_selc_old = gene_selc[genes_ind_old]
    gene_mean_old = gene_mean_old[genes_ind_old]
    gene_selc_young = gene_selc[genes_ind_young]
    gene_mean_young = gene_mean_young[genes_ind_young]
    
    
    # mean-selection spearman correlation for both age groups
    old_cor_spearman = cor(gene_selc[genes_ind_old],mean_old, use = "complete.obs",method = "spearman")
    young_cor_spearman = cor(gene_selc[genes_ind_young],mean_young, use = "complete.obs",method = "spearman")
    # mean-selection spearman correlation P values for both age groups
    p_val_old = cor.test(gene_selc[genes_ind_old],mean_old, use = "complete.obs",method = "spearman")$p.value
    p_val_young = cor.test(gene_selc[genes_ind_young],mean_young, use = "complete.obs",method = "spearman")$p.value
    
    
    # data frame containig young and old mean-selection correlation for all cell-types
    DF_cor_drop = rbind(DF_cor_drop,data.frame("Organs" = drop_organs[i],"Cell_type" = cell_types_categories[k],mean_selc_cor,pval_all,
                                               "selc_mean_old_cor_spearman" = old_cor_spearman,"selc_mean_young_spearman" = young_cor_spearman,
                                               p_val_old,p_val_young,fc_cor,pval_fc))
    
  }
}