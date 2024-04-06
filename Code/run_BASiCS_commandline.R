# Simple example of running BASiCS. Supports multiple datasets (Tabula Muris, Human data, MCA)  
# Main script for analysis of single-cell gene-expression data
print("R Version:")
print(version)
print("Start script:")

# Set paths 
user.name = "Unix" #  "Almog"  # Or # Unix 
main.dir = "/sci/labs/orzuk/orzuk/github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
code.dir <- paste0(main.dir, 'Code/')  # src files 
setwd(code.dir)
print(paste0("Current dir: ", getwd()))
source('scRNA_seq_config.R')
source("BASiCSAnalysis.R")
print("Loaded config file!")

# Parameters ordered: data.type, organ, ct_name 
args = commandArgs(trailingOnly=TRUE)  # Read from user 
data.type = args[1] 
organ = args[2]  
rerun.flag = 0 # 1: rerun everything, 0: rerun new files only, -1:  # TEMP!!! Just list which files to run!!!! 

set_data_dirs(data.type)  # set data type variables 


if(length(args)<3)  # only tissues
{
  print(paste0("Read BASiCS Input for data=", data.type, " tissue=", organ))  
  BASiCS_analysis_tissue(data.type, organ, rerun.flag)
} else  # run cell type 
{
  ct_name = args[3] 
  print(paste0("Read BASiCS Input for data=", data.type, " tissue=", organ, " cell-type=", ct_name))  
  
  # Need to read the single cell data first 
  print(paste0("Loading data from: ", paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")))  
  SC = readRDS(file = paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")) # Current tissue Seurat object - heavy code 
  print("First parsing of data:")  
  cell_types = SC$CT
  unique_cell_types = unique(cell_types)
  counts.mat = as.matrix(SC@assays$RNA@counts)  # heavy code - convert to matrix
  groups = dataset_to_age_groups(data.type)
  
  list2env(tissue_to_age_inds(data.type, organ, groups, SC@meta.data, SC), env=environment())
#  old.ind = SC$Age == "Old"
#  young.ind = SC$Age == "Young"
  if(data.type == "MCA"){  # extract for each cell the individual identity 
    batch = names(SC$orig.ident) # Get cells IDs
    for(i in 1:length(SC$orig.ident))
      batch[i] = str_split(batch[i], "\\.")[[1]][1]
  } else
    batch = SC$orig.ident
  
  print(paste0("Running BASiCS for data=", data.type, " tissue=", organ, " cell-type=", ct_name))  
  # Remove NA cells (should be done inside function?)
  test_file = Cell_type_BASiCS(data.type, organ, cell_types, ct_name, counts.mat, old.ind, young.ind, batch)
}
