library(Seurat)
library(rhdf5)
library(BASiCS)


# Run BASiCS for one tissue: loop on all cell types within the tissue 
BASiCS_analysis_tissue <- function(data.type, organ, rerun.flag = TRUE){
  groups = dataset_to_age_groups(data.type)
  samples <- get_tissue_file_names(data.type)
  organ.ind = which(samples$organs == organ)
  meta.data = get_meta_data(data.type)
  
  print("Start running BASICS!")
  read.file <- paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")
  #  print(paste0("Read file ", i, " out of ", length(samples$organs), ": ", basename(read.file)))
  SC = readRDS(file = read.file) # Current tissue Seurat object
  counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the current tissue
  list2env(tissue_to_age_inds(data.type, organ, groups, SC@meta.data, SC), env=environment()) # set specific ages for all age groups in all datasets
  
  SC_gene_name = toupper(rownames(SC)) # Gene names (in uppercase)
  if(length(SC_gene_name) != dim(counts.mat)[1]) { # For Rats, names are already in the matrix
    SC_gene_name = toupper(rownames(counts.mat)) 
  } 
  rownames(counts.mat) = SC_gene_name # make sure upper 
  
  # Dataset-specific part:   
  if(data.type == "CR.Rat"){ # Convert to dummy variables 
    cell_types = as.numeric(SC@meta.data$cell_types) # Cell types vector
    n_cell_types = max(as.numeric(SC@meta.data$cell_types)) # number of cell types. NOT WORKING for Rat
    cell_types_categories = levels(SC$cell_types)
  } 
  if(data.type %in% c("TM.facs", "TM.droplet")) # Tabula Muris
  {
    cell_types = SC@meta.data$cell.ontology.class # Cell types vector
    n_cell_types = max(SC@meta.data$cell.ontology.class) # number of cell types
    cell_types_categories = meta.data[[organ.ind]]$cell_ontology_class # Cell type names. Missing variable meta.data.drop
  }
  if(data.type == "Age.Anno"){ # Human data 
    cell_types = SC$CT
    n_cell_types = length(table(cell_types))
    cell_types_categories = names(table(SC$CT))
  }
  if(data.type == "MCA"){ # MCA mice
    cell_types = SC$CT # Cell types vector
    n_cell_types = length(table(cell_types)) # number of cell types
    cell_types_categories = names(table(SC$CT))
  }
  
  print("Filter cell types:")
  erase = vector()
  # filtering the cell-types
  for(ct_ind in 0:n_cell_types){ # why start from zero? 
    # filtering cell types with less then 100 cells or less the 20 cells in each the age groups
    if(data.type %in% c("Age.Anno", "MCA"))
      ct_name = cell_types_categories[ct_ind]
    else
      ct_name = ct_ind
      
    # Missing young or old ind !!!! 
    if(sum(cell_types==ct_name, na.rm=TRUE)<100 | 
       sum((cell_types==ct_name)&(young.ind), na.rm=TRUE) < 20 | 
       sum((cell_types==ct_name)&(!young.ind), na.rm=TRUE) < 20){
      erase = c(erase,  ct_ind+1) # ct_name+1)  # always keep in numbers 
      next()
    }
  }
  # the vector of filtered cell-types
  
  cells_ind = c(1:(n_cell_types+1))[-erase]
  if(length(cells_ind) == 0) cells_ind = (1:(n_cell_types+1))
  print("cells ind:")  
  print(cells_ind)
  
  if(data.type == "MCA")  # if(data.type == "Age.Anno")
    cells_ind = cells_ind - 1
  for(ct_ind in cells_ind){ # loop on cell types within tissue 
    print(paste0("Analyze cell type: ", cell_types_categories[ct_ind]))
    if(data.type == "Age.Anno")  # if(data.type == "Age.Anno")
      ct_name = cell_types_categories[ct_ind]
    else if(data.type == "MCA")  # if(data.type == "Age.Anno")
      ct_name = cell_types_categories[ct_ind]
    else
      ct_name = ct_ind-1  # Why minus one? 
    print(paste0("cell type category: ", cell_types_categories[ct_ind]))
    
#    print("All categories:")
#    print(cell_types_categories)
#    print("Total length categories: ")
#    print(length(cell_types_categories))
#    print("Current index: ")
#    print(ct_ind)
    
#    if(cell_types_categories[ct_ind] != "T cell")  # TEMP FOR DEBUG!!
#    {
#      print("Doesn't match cell type:")
#      print(cell_types_categories[ct_ind])
#      next
#    }
      
#    print(paste0("Run! ", cell_types_categories[ct_ind]))
#    print(ct_name)
    
    # Filter genes with less then 10 reads for either age group
    old_sum = rowSums(counts.mat[,(cell_types==ct_name & old.ind)]) 
    young_sum = rowSums(counts.mat[,(cell_types==ct_name & young.ind)])
    expressed_genes = which(old_sum > 10 & young_sum > 10)
    
#    print("Length expressed genes for tissue:")
#    print(length(expressed_genes))
#    print("Length old, length young:")
#    print(c( sum(old_sum>10), sum(young_sum>10) ))
    
    # Dataset-specific code for individual cell type       
    if(data.type == "CR.Rat")
      batch = SC$orig.ident # using mouse id as batch
    if(data.type %in% c("TM.facs", "TM.droplet"))
      batch = SC@meta.data$mouse.id # Using mouse id as batch
    if(data.type == "Age.Anno")
      batch = SC$orig.ident # using individual id as batch
    if(data.type == "Blood_SC")
      batch = SC$orig.ident  # ??? 
    if(data.type == "MCA"){  # extract for each cell the individual identity 
        batch = names(SC$orig.ident) # Get cells IDs
        for(i in 1:length(SC$orig.ident))
          batch[i] = str_split(batch[i], "\\.")[[1]][1]
      }   
    
    DVT = get_DVT_file_name(data.type, organ, cell_types_categories[ct_ind])
    print("Cell type statistics: TISSUE DIMS: COUNTS(old):")    
    print(dim(counts.mat[expressed_genes,cell_types == ct_name & old.ind]))
    print(paste0("BATCH LENs (young, old): ", length(batch[(cell_types == ct_name) & young.ind]), " , ",  
                                              length(batch[(cell_types == ct_name) & old.ind])))
    print(paste0("TABLE BATCH LENs (young, old): ", length(table(batch[(cell_types == ct_name) & young.ind])), " , ",  
                 length(table(batch[(cell_types == ct_name) & old.ind]))))
#    print("TABLE BATCH LENs (young, old): ")
#    print(length(table(batch[(cell_types == ct_name) & young.ind])))
#    print(length(table(batch[(cell_types == ct_name) & old.ind])))

    
    if((rerun.flag <= 0) & file.exists(DVT$DVT.file.name))   # Check if output file exists: !!! 
    {
      print(paste0("Skipping file = aready did this! ", ct_name))
      next
    }
      # We preform the BASiCS only if both age groups have more than 1 mouse (doesn't matter how many cells!)
    # NOTE! This filtering is only for BASICS, the over-dispersion - not for the mean analysis!
    if(!(length(table(batch[(cell_types == ct_name) & young.ind]))<=1|length(table(batch[(cell_types == ct_name) & old.ind]))<=1))
    { 
      if(rerun.flag == -1)
      {
        print(paste0("This cell type should run!!!! ", ct_name, " in tissue: ", organ))
        next
      }
        
      print(paste0("Running Cell-type! ", ct_name, " in tissue: ", organ))        # Skip cells with no BATCHES
      # BASiCS data for old and young individuals
      old_bs = newBASiCS_Data(Counts = counts.mat[expressed_genes,cell_types == ct_name & old.ind],
                              BatchInfo = batch[cell_types == ct_name & old.ind]) 
      young_bs = newBASiCS_Data(Counts = counts.mat[expressed_genes,cell_types == ct_name & young.ind],
                                BatchInfo = batch[cell_types == ct_name & young.ind]) 
      
      
      print(paste0("Chains directory: ",  DVT$basics.chains.dir))
      # Creating BASiCS chain for old and young mice
      chain_old = BASiCS_MCMC(Data = old_bs,N = 20000,Thin = 20,Burn = 10000,Regression = T,
                              WithSpikes = F,PrintProgress = FALSE, StoreChains = FALSE, StoreDir = DVT$basics.chains.dir)
      chain_young = BASiCS_MCMC(Data = young_bs,N = 20000,Thin = 20,Burn = 10000,Regression = T,
                                WithSpikes = F,PrintProgress = FALSE, StoreChains = FALSE, StoreDir = DVT$basics.chains.dir)
      print(paste0("Start BASICS TEST: : "))
      # Differential over-dispersion test - only on genes without significant difference in mean expression
      test = BASiCS_TestDE(Chain1 = chain_old,Chain2 = chain_young,GroupLabel1 = "Old",
                           GroupLabel2 = "Young",OffSet = T,PlotOffset = F,Plot = F,
                           EpsilonM = 0, EpsilonD = log2(1.5),
                           EpsilonR = log2(1.5)/log2(exp(1)), EFDR_M = 0.10, EFDR_D = 0.10) 
      save(test, file = DVT$DVT.file.name) 
    } else
      cat(paste0("Skipping, filtering this cell type!! (not enough mice for each group) ", ct_name, "\n---------\n---------\n\n")) # end if enough cells  
  }  # end loop on cell types
}  


# Run BASiCS command for a single cell type 
Cell_type_BASiCS = function(data.type, organ, cell_types, ct_name, counts.mat, old.ind, young.ind, batch, rerun.flag = TRUE){
  print("Starting cell type BASICS")
  cell_types[is.na(cell_types)] = ""  # get rid of NA cell types 
  print("Filtering cell type BASICS")
  
  # filtering the cell-types
  if(sum(cell_types==ct_name, na.rm=TRUE)<100 | 
     sum((cell_types==ct_name)&(young.ind), na.rm=TRUE) < 20 | 
     sum((cell_types==ct_name)&(!young.ind), na.rm=TRUE) < 20){
    print("The cell type don't have enough cells!")
    return()
  }
  
  print("Filter genes:")
  
  # Filter genes with less then 10 reads for either age group
  old_sum = rowSums(counts.mat[,(cell_types==ct_name & old.ind)], na.rm = TRUE) 
  young_sum = rowSums(counts.mat[,(cell_types==ct_name & young.ind)], na.rm = TRUE)
  expressed_genes = which(old_sum > 10 & young_sum > 10)
  
  DVT = get_DVT_file_name(data.type, organ, ct_name)

  if((rerun.flag <= 0) & file.exists(DVT$DVT.file.name))   # Check if output file exists: !!! 
  {
    print(paste0("Skipping file: aready did this! ", ct_name, " ; ", organ, " ; ", DVT@DVT.file.name))
    return(c())
  }
  
  
  if(!(length(table(batch[cell_types == ct_name & young.ind]))==1|length(table(batch[cell_types == ct_name & old.ind]))==1)){ 
#    print("Start BASICS CELL TYPE ANALYSIS")
#    cc = counts.mat[expressed_genes,cell_types == ct_name & old.ind]
#    print("OK COUNTS, dim cc:")
#    print(dim(cc))
#    bb = batch[(cell_types == ct_name) & old.ind]
#    print("OK BATCH, len bb:")
#    print(length(bb))
#    print(dim(counts.mat))
#    print(length(expressed_genes))
#    print(length(ct_name))
#    print(length(old.ind))
    
#    print("Length expressed genes for cell type:")
#    print(length(expressed_genes))
#    print("Legth old, length young:")
#    print(c( sum(old_sum>10), sum(young_sum>10) ))
    
    if(rerun.flag == -1)
    {
      print(paste0("This cell type should run!!!! ", ct_name, " in tissue: ", organ))
      return(c())
    }
    
    # BASiCS data for old and young mice  
    print("Run newBASiCS_Data cell-type Old and young")
    save(counts.mat, expressed_genes, cell_types, ct_name, old.ind, batch, file = "tmp_working_BASICS.Rdata")
    old_bs = newBASiCS_Data(Counts = counts.mat[expressed_genes, (cell_types == ct_name) & old.ind],
                            BatchInfo = batch[(cell_types == ct_name) & old.ind]) 
#    print("Start YOUNG")
    young_bs = newBASiCS_Data(Counts = counts.mat[expressed_genes, cell_types == ct_name & young.ind],
                              BatchInfo = batch[cell_types == ct_name & young.ind]) 
    
    print(paste0("Chains directory: ",  DVT$basics.chains.dir))
    if(!dir.exists(DVT$basics.chains.dir))
      dir.create(DVT$basics.chains.dir, recursive = TRUE)
    # Creating BASiCS chain for old and young mice
    chain_old = BASiCS_MCMC(Data = old_bs,N = 20000,Thin = 20,Burn = 10000,Regression = T,
                            WithSpikes = F,PrintProgress = FALSE, StoreChains = FALSE, StoreDir = DVT$basics.chains.dir)
    chain_young = BASiCS_MCMC(Data = young_bs,N = 20000,Thin = 20,Burn = 10000,Regression = T,
                              WithSpikes = F,PrintProgress = FALSE, StoreChains = FALSE, StoreDir = DVT$basics.chains.dir)
    
    print(paste0("Start BASICS TEST: : "))
    # Differential over-dispersion test - only on genes without significant difference in mean expression
    test = BASiCS_TestDE(Chain1 = chain_old,Chain2 = chain_young,GroupLabel1 = "Old",
                         GroupLabel2 = "Young",OffSet = T,PlotOffset = F,Plot = F,
                         EpsilonM = 0, EpsilonD = log2(1.5),
                         EpsilonR = log2(1.5)/log2(exp(1)), EFDR_M = 0.10, EFDR_D = 0.10) 
    print("Save and return:")
    save(test,file = DVT$DVT.file.name) 
    return(test)
    
  } # end if 
  
}


# Run BASiCS for one data type (loop on tissues)  
BASiCS_analysis <- function(data.type){
  # Set data directories: 
  set_data_dirs(data.type)
  samples <- get_tissue_file_names(data.type)
  meta.data = get_meta_data(data.type)
  groups = dataset_to_age_groups(data.type)
  
  
  for(i in 1:length(samples$organs)){  # loop on all tissues
    print(paste0("Read file ", i, " out of ", length(samples$organs), ": ", basename(read.file)))
    read.file <- paste0(processed.data.dir, samples$organs[i], ".", processed.files.str[data.type], ".rds")
    BASiCS_analysis_tissue(data.type, samples$organs[i])
  }  # end loop on tissues
}

