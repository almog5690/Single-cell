library(Seurat)
library(rhdf5)
library(BASiCS)

for(i in 1:length(drop_organs)){
  SC = readRDS(file = paste(drop_organs[i],"drop","rds",sep = ".")) # Reading seurat object
  counts.mat = as.matrix(SC@assays$RNA@counts) # The counts matrix
  young.ind = c(SC@meta.data$age=="3m") # Index for cells that came from young mouses
  old.ind = c(SC@meta.data$age %in% c("21m","24m")) # Index for cells that came from old mouses
  if(sum(old.ind) == 0){
    old.ind = c(SC@meta.data$age %in% c("18m","24m"))
  }
  cell_types = SC@meta.data$cell.ontology.class  
  cell_types_categories = meta.data.drop[[i]]$cell_ontology_class # Cell types names
  n_cell_types = max(cell_types) # Number of cell types
  
  for(ct_ind in 0:n_cell_types){
    #filtering cell types with less then 100 cells or less then 20 cells for each age group
    if(sum(cell_types==ct_ind)<100 | sum((cell_types==ct_ind)&(young.ind)) < 20 | sum((cell_types==ct_ind)&(old.ind)) < 20){
      next()
    } 
    
    expressed_genes = vector()
    
    for(gene in 1:n_genes){
      old = counts.mat[gene,(cell_types==ct_ind & old.ind)] # Current gene old mice count vector
      young = counts.mat[gene,(cell_types==ct_ind & young.ind)] # Current gene young mice count vector
      if(sum(old!=0)<10 | sum(young!=0)<10){
        # Filtering genes not expressed in at least 10 cells in both young and old populations
        next()
      } else {
        expressed_genes = c(expressed_genes,gene)
      }
    }
    
    batch = SC@meta.data$mouse.id # Using mouse id as batch
    
    # We preform the BASiCS only if both age group have more then 1 mouse  
    if(!(length(table(batch[cell_types == ct_ind & young.ind]))==1|length(table(batch[cell_types == ct_ind & old.ind]))==1)){ 
      
      # BASiCS data for old and young mice  
      old_bs = newBASiCS_Data(Counts = counts.mat[expressed_genes,cell_types == ct_ind & old.ind],BatchInfo = batch[cell_types == ct_ind & old.ind]) 
      young_bs = newBASiCS_Data(Counts = counts.mat[expressed_genes,cell_types == ct_ind & young.ind],BatchInfo = batch[cell_types == ct_ind & young.ind]) 
      
      # Creating BASiCS chain for old and young mice
      chain_old = BASiCS_MCMC(Data = old_bs,N = 20000,Thin = 20,Burn = 10000,Regression = T,
                              WithSpikes = F,PrintProgress = FALSE, StoreChains = TRUE, StoreDir = "/tmp/chains",
                              RunName = paste(drop_organs[i],cell_types_categories[ct_ind+1],"old")) 
      chain_young = BASiCS_MCMC(Data = young_bs,N = 20000,Thin = 20,Burn = 10000,Regression = T,
                                WithSpikes = F,PrintProgress = FALSE, StoreChains = TRUE, StoreDir = "/tmp/chains",
                                RunName = paste(drop_organs[i],cell_types_categories[ct_ind+1],"young"))  
      
      # Differential over-dispersion test - only on genes without significant difference in mean expression
      test = BASiCS_TestDE(Chain1 = chain_old,Chain2 = chain_young,GroupLabel1 = "Old",
                           GroupLabel2 = "Young",OffSet = T,PlotOffset = F,Plot = F,
                           EpsilonM = 0, EpsilonD = log2(1.5),
                           EpsilonR = log2(1.5)/log2(exp(1)),EFDR_M = 0.10, EFDR_D = 0.10) 
      
      save(test,file = paste("/tmp/DVT/DVT",drop_organs[i],cell_types_categories[k],"drop 3-24 same-mean.RData")) 
      
    }
  }
}