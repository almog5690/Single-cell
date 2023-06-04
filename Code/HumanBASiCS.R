setwd(code.dir)
library(cowplot)
library(Seurat)
library(rhdf5)
library(ggplot2)
library(BASiCS)

human.data.dir <- paste0(main.data.dir, "Data/HumanAgeAnno/Processed/")
human.celltype.dir <- paste0(main.data.dir, "Data/HumanAgeAnno/celltype/")

human_files = list.files(human.data.dir, pattern = "rds")
organs = sapply(1:length(human_files),function(x) strsplit(human_files[x],split = "[.]")[[1]][1])

add.cell.types <- TRUE
run.basics <- FALSE

# DimPlot(SC,reduction = )

# Run once to add cell type information to Human data: 
# 'anno' is taken from the variable: processed.files.str
if(add.cell.types)
  for(i in 5:length(organs)){
    print(i)
    print(organs[i])
    SC = readRDS(file = paste0(human.data.dir, paste(organs[i],"anno.rds",sep = ".")))
    cell_types = as.data.frame(read.delim(paste0(human.celltype.dir, paste(organs[i],"txt",sep = ".")),header = T))
    SC$CT = sapply(1:ncol(SC),function(x) if(names(Idents(SC))[x] %in% cell_types[,1]){
      cell_types$ident[which(cell_types[,1] == names(Idents(SC))[x])]} else{NA})
    # ct_test = sapply(1:ncol(SC),function(x) cell_types$ident[which(cell_types[,1] == names(Idents(SC))[x])])
    saveRDS(SC, file = paste0(human.data.dir, paste(organs[i],"anno.rds",sep = ".")))
  }
# which(duplicated(cell_types[,1]))

if(run.basics)
  #  for(i in 1:length(organs)){
  for(i in 5:5){
    SC = readRDS(file = paste0(human.data.dir, paste(organs[i],"anno.rds",sep = ".")))
    BASiCS_analysis_tissue("Age.Anno", organs[i]) # Replace this by separate jobs 
  }
