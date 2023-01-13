library(Seurat)
library(rhdf5)
library(BASiCS)

BASiCS_analysis <- function(data.type){
  samples <- get_tissue_file_names(data.type)
  meta.data = get_meta_data(data.type)
  
  
  
  for(i in 1:length(samples$organs)){  # for(i in 1:2){
    read.file <- paste0(processed.data.dir, samples$organs[i], ".", processed.files.str[data.type], ".rds")
    print(paste0("Read file ", i, " out of ", length(samples$organs), ": ", basename(read.file)))
    SC = readRDS(file = read.file) # Current tissue seurat object
    counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the current tissue
    list2env(tissue_to_age_inds(data.type, samples$organs[i], groups, SC@meta.data), env=environment()) # set specific ages for all age groups in all datasets
    
    SC_gene_name = toupper(rownames(SC)) # Gene names (in uppercase)
    if(length(SC_gene_name) != dim(counts.mat)[1]) { # For Rats, names are already in the matrix
      SC_gene_name = toupper(rownames(counts.mat)) 
    } 
    rownames(counts.mat) = SC_gene_name # make sure upper 
    
    if(data.type == "CR.Rat"){ # Convert to dummy variables 
      cell_types = as.numeric(SC@meta.data$cell_types) # Cell types vector
      cell_types_categories = levels(SC$cell_types)
    } else
    {
      cell_types = SC@meta.data$cell.ontology.class # Cell types vector
      cell_types_categories = meta.data[[i]]$cell_ontology_class # Cell type names. Missing variable meta.data.drop
    }
    
    for(CT_ind in cells_ind){ # loop on cell types . 
      print(paste0("Analyze cell type: ", cell_types_categories[k]))
      
      # Filter genes with less then 10 reads for either age group
      old_sum = rowSums(counts.mat[,(cell_types==ct_ind & old.ind)]) 
      young_sum = rowSums(counts.mat[,(cell_types==ct_ind & young.ind)])
      expressed_genes = which(old_sum > 10 & young_sum > 10)
      
      if(data.type == "CR.Rat"){
        batch = SC$orig.ident # using mouse id as batch
      } else {
        batch = SC@meta.data$mouse.id # Using mouse id as batch
      }
      
      # BASiCS directory and file names
      if(data.type == "CR.Rat"){
        basics.chains.dir =  paste0(main.data.dir, 'Data/CR.Rat/chains/')
        basics.DVT.dir =  paste0(main.data.dir, 'Data/CR.Rat/DVT/')
        DVT.file.name = paste0(basics.DVT.dir,paste(samples$organs[i],cell_types_categories[ct_ind + 1],"same-mean.RData"))

      }
      if(data.type == "TM.droplet"){
        basics.chains.dir =  paste0(main.data.dir, 'Data/TabulaMuris/chains/droplet')
        basics.DVT.dir =  paste0(main.data.dir, 'Data/TabulaMuris/DVT/droplet')
        DVT.file.name = paste0(basics.DVT.dir,paste(samples$organs[i],cell_types_categories[ct_ind + 1],"drop 3-24 same-mean.RData"))
      }
      if(data.type == "TM.facs"){
        basics.chains.dir =  paste0(main.data.dir, 'Data/TabulaMuris/chains/facs')
        basics.DVT.dir =  paste0(main.data.dir, 'Data/TabulaMuris/DVT/facs')
        DVT.file.name = paste0(basics.DVT.dir,paste(samples$organs[i],cell_types_categories[ct_ind + 1],"same-mean.RData"))
      }
      
      # We preform the BASiCS only if both age group have more then 1 mouse  
      if(!(length(table(batch[cell_types == ct_ind & young.ind]))==1|length(table(batch[cell_types == ct_ind & old.ind]))==1)){ 
        
        # BASiCS data for old and young mice  
        old_bs = newBASiCS_Data(Counts = counts.mat[expressed_genes,cell_types == ct_ind & old.ind],BatchInfo = batch[cell_types == ct_ind & old.ind]) 
        young_bs = newBASiCS_Data(Counts = counts.mat[expressed_genes,cell_types == ct_ind & young.ind],BatchInfo = batch[cell_types == ct_ind & young.ind]) 
        
        # Creating BASiCS chain for old and young mice
        chain_old = BASiCS_MCMC(Data = old_bs,N = 20000,Thin = 20,Burn = 10000,Regression = T,
                                WithSpikes = F,PrintProgress = FALSE, StoreChains = TRUE, StoreDir = basics.chains.dir,
                                RunName = paste(samples$organs[i],cell_types_categories[k],data.type,"old")) 
        chain_young = BASiCS_MCMC(Data = young_bs,N = 20000,Thin = 20,Burn = 10000,Regression = T,
                                  WithSpikes = F,PrintProgress = FALSE, StoreChains = TRUE, StoreDir = basics.chains.dir,
                                  RunName = paste(samples$organs[i],cell_types_categories[k],data.type,"young"))  
        
        # Differential over-dispersion test - only on genes without significant difference in mean expression
        test = BASiCS_TestDE(Chain1 = chain_old,Chain2 = chain_young,GroupLabel1 = "Old",
                             GroupLabel2 = "Young",OffSet = T,PlotOffset = F,Plot = F,
                             EpsilonM = 0, EpsilonD = log2(1.5),
                             EpsilonR = log2(1.5)/log2(exp(1)),EFDR_M = 0.10, EFDR_D = 0.10) 
        
        save(test,file = DVT.file.name) 
        
      }
    }
  }
}

