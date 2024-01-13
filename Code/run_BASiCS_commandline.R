# Simple example of running BASiCS. Try Tabula Muris or Human data, or MCA  
# Main script for analysis of single-cell gene-expression data

print("R Version:")
print(version)
print("Start script:")

# Parameters ordered: data.type, organ, ct_name 
args = commandArgs(trailingOnly=TRUE)  # Read from user 
data.type = args[1] 
organ = args[2]  
ct_name = args[3] 

print(paste0("Read Input for data=", data.type, " tissue=", organ, " cell-type=", ct_name))  

# Set paths 
user.name = "Unix" #  "Almog"  # Or # Unix 
main.dir = "/sci/labs/orzuk/orzuk/github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
code.dir <- paste0(main.dir, 'Code/')  # src files 
setwd(code.dir)
print(paste0("Current dir: ", getwd()))
source('scRNA_seq_config.R')
print("Loaded config file!")

# Need to read the single cell data first 
print(paste0("Loading data from: ", paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")))  
SC = readRDS(file = paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")) # Current tissue Seurat object - heavy code 
print("First parsing of data:")  
cell_types = SC$CT
unique_cell_types = unique(cell_types)
counts.mat = as.matrix(SC@assays$RNA@counts)  # heavy code - convert to matrix
old.ind = SC$Age == "Old"
young.ind = SC$Age == "Young"
batch = SC$orig.ident

print(paste0("Running BASiCS for data=", data.type, " tissue=", organ, " cell-type=", ct_name))  
# Remove NA cells (should be done inside function?)
test_file = Cell_type_BASiCS(data.type, organ, cell_types, ct_name, counts.mat, old.ind, young.ind, batch)
