SC = readRDS(file = paste(organs[1],"rds",sep = "."))
n_genes = SC@assays$RNA@counts@Dim[1] # Number of genes
nam = toupper(rownames(SC)) # Genes names
od_base = rep(NA,n_genes) # empty vector
names(od_base) = nam  

OD_old_data = c()
OD_young_data = c()
OD_log_data = c()
for(i in 1:length(facs.files)){
  SC = readRDS(file = paste(organs[i],"rds",sep = ".")) # current tissue seurat object
  
  young.ind = c(SC@meta.data$age %in% c("3m")) # index for cells that came from young mouses
  n_cell = SC@assays$RNA@counts@Dim[2]
  SC_gene_name = toupper(rownames(SC)) # genes names in upper case letters
  
  cell_types = SC@meta.data$cell.ontology.class # the cell-types vector 
  cell_types_categories = meta.data[[i]]$cell_ontology_class # the names of the different cell-types
  n_cell_types = max(cell_types) # number of cell types
  
  earase = vector()
  # filltering the cell-types
  for(ct_ind in 0:n_cell_types){
    #filtering cell types with less then 100 cells or less the 20 cells in each the age groups
    if(sum(cell_types==ct_ind)<100 | sum((cell_types==ct_ind)&(young.ind)) < 20 | sum((cell_types==ct_ind)&(!young.ind)) < 20){
      earase = c(earase,ct_ind+1)
      next()
    }
    
  }
  
  
  # the vector of filltered cell-types
  cells_ind = c(1:(n_cell_types+1))[-earase]
  if(length(cells_ind) == 0) cells_ind = (1:(n_cell_types+1))
  
  for(k in cells_ind){
    test_file = paste("/tmp/DVT/DVT",organs[i],cell_types_categories[ct_ind + 1],"same-mean.RData") # BASiCS over-dispersion young vs old test file
    if(file.exists(test_file)){
      load(file = test_file) # Loading the BASiCS over-dispersion young vs old test
      df = test@Results$Disp@Table # the test result table
      df$GeneName = toupper(df$GeneName) 
      df = df[df$ResultDiffDisp != "ExcludedFromTesting",] # removing genes with different mean between young and old
      
      od_old = od_base
      od_old[df$GeneName] = df$Disp1 # Current cell type old over-dispersion for all genes
      
      od_young = od_base
      od_young[df$GeneName] = df$Disp2 # Current cell type young over-dispersion for all genes
      
      OD_old_data = cbind(OD_old_data,od_old)
      OD_young_data = cbind(OD_young_data,od_young)
      
      OD_log_data = rbind(OD_log_data,data.frame("Organs" = organs[i],"Cell_type" = cell_types_categories[k],
                                                 "Genes" = df$GeneName,"Log_OD" = df$DispLog2FC,"OD" = df$DispFC,"significant" = df$ResultDiffDisp))
    }
  }
}

# Finding the most over dispersed genes
genes = unique(OD_log_data$Genes) 
genes_to_keep = c()
for(j in 1:length(genes)){
  temp_data = OD_log_data[OD_log_data$Genes == genes[j],]
  
  if(sum(temp_data$significant %in% c("Old+","Young+")) > 20){ # If current gene is significant in more then 20 cell types
    if(sum(temp_data$Log_OD > log2(1.5) | temp_data$Log_OD < log2(0.75)) > 20){ # If current gene over-dispersion for old is more then 1.5 then the young (or the opposite) in more then 20 cell type
      genes_to_keep=c(genes_to_keep,genes[j])
    }
  }
}
OD_DF = OD_log_data[OD_log_data$Genes %in% genes_to_keep,]


# Getting the Log over-dispersion (old/young) matrix (genes by cell types)
genes = unique(OD_DF$Genes) # remainig genes
ct = unique(interaction(OD_DF$Organs,OD_DF$Cell_type,sep = ": ")) # cell types
ct = factor(ct,levels = ct)
OD_DF$ct = interaction(OD_DF$Organs,OD_DF$Cell_type,sep = ": ") # Organs-cell types interaction
OD_DF$ct = factor(OD_DF$ct,levels = ct)
OD_log_mat = matrix(NA,length(genes),length(ct))
rownames(OD_log_mat) = genes
colnames(OD_log_mat) = ct
for(j in 1:length(ct)){
  temp = OD_DF[OD_DF$ct == ct[j],]
  OD_log_mat[temp$Genes,j] = temp$Log_OD
}

OD_log_mat[is.na(OD_log_mat)] = 0 # replacing NA with 0

hp_facs = heatplot(t(OD_log_mat*100),returnSampleTree = T) # creating the heatmap

# Getting the genes clusters from the heat map 
hclust = as.hclust(hp_facs) # hclust object
clusts = cutree(hclust,k=3) # finding the 3 clusters from the hclust 
l = labels(hp_facs)
oc = clusts[l] # Vector of the cluters. 