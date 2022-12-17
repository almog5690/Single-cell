library(Seurat)
library(rhdf5)

#main.dir = "C:/Code/Github/Single-cell"  # Change to your local path. This path should be used everywhere
# setwd(main.dir)

preprocess_droplet <- function()
{
  droplet.files = list.files(path = raw.data.dir, pattern = "drop.r")  # Raw data should be in the data sub-directory
  droplet.files = droplet.files[-c(6,21,2,15,16)] # Removing Large intestine, Fat, Skin, Pancreas and Trachea because they have only one age group.
  drop_organs = sapply(droplet.files,function(f) unlist(strsplit(f,"[.]"))[1]) # List of droplet organs
  drop_organs = unname(drop_organs)
  
  meta.data.drop = list() # getting organs metadata 
  for(i in 1:length(droplet.files)){
    print(c("Opening file:", paste(drop_organs[i],"droplet.h5ad",sep = "_")))
    meta = h5read(paste0("raw files/",droplet.files[i]),"uns")
    meta = meta[grep("categories",names(meta))]
    names(meta) = sapply(names(meta),function(name) substr(name,1,nchar(name) - 11))
    meta.data.drop[[i]] = meta
  }
  
  for(i in 1:length(drop_organs)){
    # Reading the meta data categories
    meta = h5read(paste0("/tmp/raw files/",drop.files[i]),"uns")
    meta = meta[grep("categories",names(meta))]
    names(meta) = sapply(names(meta),function(name) substr(name,1,nchar(name) - 11))
    
    SC = ReadH5AD(paste0("/tmp/raw files/",drop.files[i]),assay = "RNA",verbose = F) # Reading the tissue file in to Seurat object
    
    # Normalization and feature selection.
    SC = NormalizeData(SC) # Size factor normalization to 10,000 counts and log transformation
    SC = FindVariableFeatures(SC,verbose = F) # Finding the 2000 highly variable features
    
    SC = ScaleData(SC) # Scaling the data
    
    SC@active.ident = factor(SC@meta.data$cell.ontology.class,labels = meta$cell_ontology_class) # Setting the active ident to cell types 
    SC = AddMetaData(SC,metadata = ifelse(SC@meta.data$age%in%c("18m","21m","24m","30m"),"Old","Young"),col.name = "Age") # Adding Old/Young vector to meta data 
    SC = AddMetaData(SC,metadata = SC@active.ident,col.name = "Cells") # Adding cell type names to meta data
    SC$age = factor(SC$age,labels = meta$age) # Naming the ages vector in meta data
    
    # Running PCA, TSNE and UMAP
    SC = RunPCA(SC,verbose = F) 
    SC = RunTSNE(SC) 
    SC = RunUMAP(SC,dims = 1:10) 
    
    saveRDS(SC, file = paste(drop_organs[i],"drop","rds",sep = ".")) # Saving the seurat object
  }   # end loop on tissues 
  
  # save the drop_organs variable and possibly additional variables to an RData/RDS file. This will be loaded next time without running the preprocessing again.
  
}  # end function preprocess_droplet


# Preprocess the facs dataset 
preprocess_facs <- function()
{
  facs.files = list.files(path = raw.data.dir, pattern = "facs",full.names = T) # facs raw files list
  
  organs  = sapply(facs.files,function(f) strsplit(f,split ="-|[.]")[[1]][[8]]) # List of facs organs
  organs = unname(organs)
  
  # getting organs metadata 
  meta.data = list() 
  for(i in 1:length(facs.files)){
    meta = h5read(facs.files[i],"uns")
    meta = meta[grep("categories",names(meta))]
    names(meta) = sapply(names(meta),function(name) substr(name,1,nchar(name) - 11))
    meta.data[[i]] = meta
  }
  
  # Pre-processing the raw data + downsampling
  for(i in 1:length(facs.files)){
    # Reading the meta data categories
    meta = h5read(facs.files[i],"uns") 
    meta = meta[grep("categories",names(meta))]
    names(meta) = sapply(names(meta),function(name) substr(name,1,nchar(name) - 11))
    
    SingleCell_tissue = ReadH5AD(facs.files[i],assay = "RNA",verbose = F) # Reading tissue file in to Seurat object
    names(Idents(SingleCell_tissue)) = Cells(SingleCell_tissue) # Setting the ident to cell types
    
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
    SC = SetAssayData(SingleCell_tissue,slot = "counts",new.data = count) # Creating a new Seurat object with the new counts matrix
    
    # Normalization and feature selection.
    SC = NormalizeData(SC) # Size factor normalization to 10,000 counts and log transformation
    SC = FindVariableFeatures(SC,verbose = F) # Finding the 2000 highly variable features 
    
    SC = ScaleData(SC,verbose = F,do.center = F) # Scaling the data
    
    SC@active.ident = factor(SC@meta.data$cell.ontology.class,labels = meta$cell_ontology_class) # Setting the active ident to cell types 
    SC = AddMetaData(SC,metadata = ifelse(SC@meta.data$age%in%c("18m","21m","24m","30m"),"Old","Young"),col.name = "Age") #Adding Old/Young vector to meta data 
    SC = AddMetaData(SC,metadata = SC@active.ident,col.name = "Cells") # Adding cell type names to meta data
    SC$age = factor(SC$age,labels = meta$age) # naming the ages vector in meta data
    
    # Running PCA,TSNE and UMAP
    SC = RunPCA(SC,verbose = F) 
    SC = RunTSNE(SC) 
    SC = RunUMAP(SC,dims = 1:10)  
    
    saveRDS(SC,file = paste(organs[i],"rds",sep = ".")) # Saving the Seurat object 
  } # end loop on facs files 
} # End function preprocess facs

