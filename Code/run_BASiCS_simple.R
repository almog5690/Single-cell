# Simple example of running BASiCS. Try Tabula Muris or Human data  
# Main script for analysis of single-cell gene-expression data
data.type = "TM.facs"  # Blood_SC
organ = "Aorta" # Blood (?) # smallest facs file, should be fastest 
ct_name = "NK-FCER1G" # "ABC" "Naive-B" ## cell type example - change to a real name!! 
data.type = "Blood_SC"
organ = "Blood"

run.tissue = FALSE
run.celltype = !(run.tissue)


if(!exists("user.name"))  # set default values 
  user.name = "Or" #  "Almog"  # Or # Unix 
if(!exists("data.type"))  # set default values 
  data.type = "TM.facs"  # Choose one type for analysis (can later loop over multiple datasets)
if(user.name == "Almog")
  main.dir = "C:/Users/User/OneDrive/Documents/Github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
if(user.name == "Or")
  main.dir = "C:/Code/Github/Single-cell/" 
# Change to your local path. This path should be used everywhere as the main source path
if(user.name == "Unix")
  main.dir = "/sci/labs/orzuk/orzuk/github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
code.dir <- paste0(main.dir, 'Code/')  # src files 
setwd(code.dir)
source('scRNA_seq_config.R')

if(run.tissue)
{
  print("Running BASiCS:")
  BASiCS_analysis_tissue(data.type, organ)
  print("Finished Running BASiCS!")
}

if(run.celltype)
{ ### For the Human-Blood:
  # Need to read the single cell data first 
  SC = readRDS(file = paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")) # Current tissue Seurat object - heavy code 
  cell_types = SC$CT
  unique_cell_types = unique(cell_types)
  counts.mat = as.matrix(SC@assays$RNA@counts)  # heavy code- convert to matrix
  old.ind = SC$Age == "Old"
  young.ind = SC$Age == "Young"
  batch = SC$orig.ident
  test_file = Cell_type_BASiCS(data.type, organ, cell_types, ct_name, counts.mat, old.ind, young.ind, batch)
}
