library(Seurat)
library(rhdf5)

#main.dir = "C:/Code/Github/Single-cell"  # Change to your local path. This path should be used everywhere
# setwd(main.dir)


get_meta_data <- function(data.type, force.rerun = FALSE)
{
  meta.data.file <- paste0(processed.data.dir, 'meta.data.', processed.files.str[data.type], '.RData')
  if(file.exists(meta.data.file) & (force.rerun == FALSE)) {
       load(meta.data.file)
       return(meta.data)
  }

  if(data.type == "TM.droplet")
  {
    meta.data = list() # getting organs metadata 
    for(i in 1:length(droplet.files)){
      print(c("Opening file:", paste(drop_organs[i],"droplet.h5ad",sep = "_")))
      meta = h5read(paste0("raw files/",droplet.files[i]),"uns")
      meta = meta[grep("categories",names(meta))]
      names(meta) = sapply(names(meta),function(name) substr(name,1,nchar(name) - 11))
      meta.data[[i]] = meta
    }
  }
  if(data.type == "TM.facs") # Fill in .. 
  {
    # getting organs metadata 
    meta.data = list() 
    for(i in 1:length(facs.files)){
      meta = h5read(facs.files[i],"uns")
      meta = meta[grep("categories",names(meta))]
      names(meta) = sapply(names(meta),function(name) substr(name,1,nchar(name) - 11))
      meta.data[[i]] = meta
    }
  }
  if(data.type == "Blood_SC"){
    meta.data = list("cell.ontology.class" = c("ABC", "CD14-MC", "CD14-MC-PPBP", "CD16-MC", "CD4-Naive",        
                  "CD4-Naive-NECTIN1", "CD4-Tem", "CD4-Tm", "CD4-Treg", "CD8-CTL",          
                  "CD8-Naive", "CD8-Tem", "cDC", "Intermed-MC", "Megakaryocytes",   
                  "Memory-B-CRIP1", "Memory-B-IL32", "Naive-B", "Naive-B-NRGN", "NK-FCER1G",        
                  "NK-GZMH", "NK-GZMK", "NK-IL7R", "NK-S100A2", "NK-SELL",          
                  "pDC", "Plasma", "RBC"))
  }
  save(meta.data, file = meta.data.file)  
  return(meta.data)
}


preprocess <- function(data.type){
  # Set data directories: 
  set_data_dirs(data.type)
  
  if(data.type == "TM.droplet"){
    raw.files = list.files(path = raw.data.dir, pattern = "drop.r",full.names = T)
    raw.files = raw.files[-c(6,21,2,15,16)]  # Removing Large intestine, Fat, Skin, Pancreas and Trachea because they have only one age group.
    organs = sapply(raw.files,function(f) unlist(strsplit(f,"[.]"))[1]) # List of droplet organs
    organs = unname(organs)
  }
  if(data.type == "TM.facs"){
    raw.files = list.files(path = raw.data.dir, pattern = "facs",full.names = T) # facs raw files list
    
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

