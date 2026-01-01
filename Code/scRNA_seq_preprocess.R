library(Seurat)
library(rhdf5)

#main.dir = "C:/Code/Github/Single-cell"  # Change to your local path. This path should be used everywhere
# setwd(main.dir)


# Read meta data for each data type
get_meta_data <- function(data.type, force.rerun = FALSE)
{
  set_data_dirs(data.type)
  meta.data.file <- paste0(processed.data.dir, 'meta.data.', processed.files.str[data.type], '.RData')
  if(file.exists(meta.data.file) & (force.rerun == FALSE)) {
    print(paste0("Metadata Exists! ", meta.data.file))
    load(meta.data.file)
    return(meta.data)
  }

  if(data.type == "TM.droplet")
  {
    meta.data = list() # getting organs metadata 
    for(i in 1:length(raw.files)){
      print(c("Opening file:", paste(drop_organs[i],"droplet.h5ad", sep = "_")))
      meta = h5read(paste0("raw files/", raw.files[i]),"uns")
#      meta = h5read(paste0("raw files/",droplet.files[i]),"uns")
      meta = meta[grep("categories",names(meta))]
      names(meta) = sapply(names(meta),function(name) substr(name,1,nchar(name) - 11))
      meta.data[[i]] = meta
    }
  }
  if(data.type == "TM.facs") # organs meta data 
  {
    # getting organs metadata 
    meta.data = list() 
    for(i in 1:length(raw.files)){
      meta = h5read(raw.files[i], "uns")
      meta = meta[grep("categories",names(meta))]
      names(meta) = sapply(names(meta),function(name) substr(name,1,nchar(name) - 11))
      meta.data[[i]] = meta
    }
  }
  if(data.type == "Blood_SC"){
    meta.data = list(list("cell_ontology_class" = c("ABC", "CD14-MC", "CD14-MC-PPBP", "CD16-MC", "CD4-Naive",        
                  "CD4-Naive-NECTIN1", "CD4-Tem", "CD4-Tm", "CD4-Treg", "CD8-CTL",          
                  "CD8-Naive", "CD8-Tem", "cDC", "Intermed-MC", "Megakaryocytes",   
                  "Memory-B-CRIP1", "Memory-B-IL32", "Naive-B", "Naive-B-NRGN", "NK-FCER1G",        
                  "NK-GZMH", "NK-GZMK", "NK-IL7R", "NK-S100A2", "NK-SELL",          
                  "pDC", "Plasma", "RBC")))
  }
  
  if(data.type == "CR.Rat"){ # Rats: (meta data file ready)
     
  }
  if(data.type == "MCA")  # mice new data 
  {
    meta.data = list() # getting organs metadata 
  }
  
  save(meta.data, file = meta.data.file)  
  return(meta.data)
}


preprocess <- function(data.type){
  # Set data directories: 
  set_data_dirs(data.type)
  
  if(data.type == "TM.droplet"){
    raw.files = list.files(path = raw.data.dir, pattern = "drop.r", full.names = T)
    raw.files = raw.files[-c(6,21,2,15,16)]  # Removing Large intestine, Fat, Skin, Pancreas and Trachea because they have only one age group.
    organs = sapply(raw.files,function(f) unlist(strsplit(f,"[.]"))[1]) # List of droplet organs
    organs = unname(organs)
  }
  if(data.type == "TM.facs"){
    raw.files = list.files(path = raw.data.dir, pattern = "facs", full.names = T) # facs raw files list
    
    organs  = sapply(raw.files,function(f) strsplit(f,split ="-|[.]")[[1]][[8]]) # List of facs organs
    organs = unname(organs)
  }
  
  # Pre-processing the raw data + downsampling
  for(i in 1:length(raw.files)){
    # Reading the meta data categories
    meta = h5read(raw.files[i],"uns") 
    meta = meta[grep("categories",names(meta))]
    names(meta) = sapply(names(meta),function(name) substr(name,1,nchar(name) - 11))
    
    SingleCell_tissue = ReadH5AD(raw.files[i],assay = "RNA",verbose = F) # Reading tissue file in to Seurat object
    names(Idents(SingleCell_tissue)) = Cells(SingleCell_tissue) # Setting the ident to cell types
    
    if(data.type == "TM.facs"){
      # Down sampling
      NCount_avg_old = mean(SingleCell_tissue$nCount_RNA[SingleCell_tissue$age!=0]) # Old mice mean count
      NCount_avg_young = mean(SingleCell_tissue$nCount_RNA[SingleCell_tissue$age==0]) # Young mice mean count
      p = NCount_avg_young/NCount_avg_old # Young-old mean count ratio
      
      if(p<1){ # If p<1 we preform down sampling
        count = SingleCell_tissue[["RNA"]]@counts # Counts matrix 
        for(cell in which(SingleCell_tissue$age>0)){
          # For each old mice cell we randomly choose (with binomial distribution) the count to keep base on p
          count[,cell] = rbinom(nrow(count),count[,cell],p)
        }
      } else{
        count = SingleCell_tissue[["RNA"]]@counts
      }
      SingleCell_tissue = SetAssayData(SingleCell_tissue,slot = "counts",new.data = count) # Creating a new Seurat object with the new counts matrix
      
    }
        
    # Normalization and feature selection.
    SingleCell_tissue = NormalizeData(SingleCell_tissue) # Size factor normalization to 10,000 counts and log transformation
    SingleCell_tissue = FindVariableFeatures(SingleCell_tissue,verbose = F) # Finding the 2000 highly variable features 
    
    SingleCell_tissue = ScaleData(SingleCell_tissue,verbose = F,do.center = F) # Scaling the data
    
    SingleCell_tissue@active.ident = factor(SingleCell_tissue@meta.data$cell.ontology.class,labels = meta$cell_ontology_class) # Setting the active ident to cell types 
    SingleCell_tissue = AddMetaData(SingleCell_tissue,metadata = ifelse(SingleCell_tissue@meta.data$age%in%c("18m","21m","24m","30m"),"Old","Young"),col.name = "Age") #Adding Old/Young vector to meta data 
    SingleCell_tissue = AddMetaData(SingleCell_tissue,metadata = SingleCell_tissue@active.ident,col.name = "Cells") # Adding cell type names to meta data
    SingleCell_tissue$age = factor(SingleCell_tissue$age,labels = meta$age) # naming the ages vector in meta data
    
    # Running PCA,TSNE and UMAP
    SingleCell_tissue = RunPCA(SingleCell_tissue,verbose = F) 
    SingleCell_tissue = RunTSNE(SingleCell_tissue) 
    SingleCell_tissue = RunUMAP(SingleCell_tissue,dims = 1:10)  
    
    saveRDS(paste0(processed.data.dir,SingleCell_tissue,file = paste(organs[i],"rds",sep = "."))) # Saving the Seurat object 
  } # end loop on facs files 
  
} # end function preprocess



# Filter cells from an expression matrix (not used yet. Should be part of analysis)
filter_cells <- function(cell_types, young.ind, old.ind, filter.params)
{
  unique_cell_types = unique(cell_types)
  n_cell_types = length(unique_cell_types)  #  max(cell_types) # Number of cell types
  erase = vector()
  
  for(j in 1:n_cell_types){ # filtering cell types with less then 100 cells or less the 20 cells in each the age groups
    ct_ind = unique_cell_types[j]
    if(sum(cell_types==ct_ind, na.rm=TRUE) < filter.params$min.cells.total |
       sum(cell_types==ct_ind & young.ind, na.rm=TRUE) < filter.params$min.cells.per.age |
       sum(cell_types==ct_ind & old.ind, na.rm=TRUE) < filter.params$min.cells.per.age){
      erase = c(erase, j)
      next()
    }
  }
  cells_ind = unique_cell_types[-erase] #  cells_ind = c(1:(n_cell_types+1))[-erase]
  if(length(cells_ind) == 0) cells_ind = unique_cell_types # (1:(n_cell_types+1))
  return(cells_ind) # indices of cells that we keep 
}




# Get fold-change values from young and old expression values
#
get_fold_change <- function(young_exp, old_exp, fc.method, 
                            filter.params, stat, SC=c(), data.type="TM.facs", CT_name="")
{
  
  if(fc.method == "log_old_minus_young")
    fc_exp = log((old_exp+filter.params$pseudo.count)/(young_exp+filter.params$pseudo.count))
  if(fc.method == "seurat")
  {
    age.group.str = dataset_to_age_groups(data.type)
    if(stat != "mean")
    {
      print("Error! cannot compute ", stat, " fold-change using Seurat")
      return(NA)
    }
    Idents(SC) <- "Cells"
    CT_SC = subset(SC, idents = CT_name) # current cell type Seurat object
    
    # Clustering using FindMarkers
    age_clusters = FindMarkers(CT_SC, ident.1 = Cells(CT_SC)[CT_SC$age %in% age.group.str$old_ages_1],
                               ident.2 = Cells(CT_SC)[CT_SC$age %in% age.group.str$young_ages],
                               test.use = "wilcox", assay = "RNA", slot = "data",
                               pseudocount.use = 0.1, verbose = T, min.pct = 0.5, logfc.threshold = log2(1.25))
    age_clusters$bh_p_val = p.adjust(age_clusters$p_val, method = "BH") # Adjusted p_values
    fc_exp = age_clusters$avg_log2FC # FC estimations
    names(fc_exp) = toupper(row.names(age_clusters)) # gene names
  }
  if(fc.method == "basics")  # works for mean and overdispersion
  {
    load(BASiCS.files$test) # load BASICS object 
    if(stat == "mean")
    {
      df.stat = test@Results$Mean@Table
      fc_exp = df.stat$MeanFC  # works for mean
    }
    if(stat %in% c("od", "overdispersion"))
    {
      df.stat = test@Results$Disp@Table
      fc_exp = df.stat$DispFC  # works for mean
    }      
    names(fc_exp) = toupper(df.stat$GeneName) # Uppercasing gene names
  }
  
  return(fc_exp)
}


# Mean filtering function
# Input: 
#
# Output: 
# vector of mean expression of genes 
gene_filtering = function(gene_mean_old, gene_mean_young, gene_mean_all = NULL,
                          expression.thresh = 0.2, FC_filter = F, SC, CT_name, 
                          data.type){
  
  age.group.str = dataset_to_age_groups(data.type)
  if(FC_filter){ # FC-filtering using FindMarkers function
    Idents(SC) <- "Cells"
    CT_SC = subset(SC, idents = CT_name) # current cell type Seurat object
    
    # if(all(is.na( Cells(CT_SC)[CT_SC$age %in% Old.indicator]))) next()
    # Clustering using FindMarkers
    age_clusters = FindMarkers(CT_SC, ident.1 = Cells(CT_SC)[CT_SC$age %in% age.group.str$old_ages_1],
                               ident.2 = Cells(CT_SC)[CT_SC$age %in% age.group.str$young_ages],
                               test.use = "wilcox",assay = "RNA",slot = "data",
                               pseudocount.use = 0.1,verbose = T, min.pct = 0.5,logfc.threshold = log2(1.25))
    age_clusters$bh_p_val = p.adjust(age_clusters$p_val,method = "BH") # Adjusted p_values
    gene_mean_FC = age_clusters$avg_log2FC # FC estimations
    names(gene_mean_FC) = toupper(row.names(age_clusters)) # gene names
    
    return(gene_mean_FC[age_clusters$bh_p_val <= 0.1])
    
  } else { # Mean filtering (young/old or all)
    ages_mean = (gene_mean_old + gene_mean_young)/2 # genes averages average  
    
    if(is.null(gene_mean_all)){ # young/old filtering
      return(c("Young" = gene_mean_young[which(ages_mean > expression.thresh)], 
               "Old" = gene_mean_old[which(ages_mean > expression.thresh)]))
      
    } else { # all filtering
      return(gene_mean_all[which(ages_mean > expression.thresh)])
    }
    
  }
}




