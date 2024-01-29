# Simple example of running BASiCS. Try Tabula Muris or Human data  
# Main script for analysis of single-cell gene-expression data
data.type = "MCA" # "TM.facs"  # Blood_SC
organ = "Liver" # skin has a small file "Aorta" # Blood (?) # smallest facs file, should be fastest 
ct_name = "T cell" # "Acinar cell" # "Fibroblast" # "NK-FCER1G" # "ABC" "Naive-B" ## cell type example - change to a real name!! 
# data.type = "Blood_SC"
# organ = "Blood"
# ct_name = "NK-FCER1G" # "ABC" "Naive-B" ## cell type example - change to a real name!! 
# data.type = "CR.Rat"
# organ = "Aorta"
# ct_name = "NK-FCER1G" # "ABC" "Naive-B" ## cell type example - change to a real name!! 



run.tissue = TRUE
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
  
  list2env(tissue_to_age_inds(data.type, organ, groups, SC@meta.data), env=environment())
  
  #old.ind = SC$Age == "Old"
  #young.ind = SC$Age == "Young"
  
  if(data.type == "MCA"){  # extract for each cell the individual identity 
    batch = names(SC$orig.ident) # Get cells IDs
    for(i in 1:length(SC$orig.ident))
      batch[i] = str_split(batch[i], "\\.")[[1]][1]
  } else
    batch = SC$orig.ident
  
  # Remove NA cells (should be done inside function?)
  #  na.inds = which(is.na(counts.mat), arr.ind=TRUE)
  #  na.cols = union(na.inds[,2])
  test_file = Cell_type_BASiCS(data.type, organ, cell_types, ct_name, counts.mat, old.ind, young.ind, batch)
}

prepare.scripts = TRUE
if(prepare.scripts)  # Make command line scripts for running on unix cluster 
{
  # get all organ lists! 
  processed.files = list.files(path = processed.data.dir, pattern = paste0(processed.files.str[data.type], ".rds"), full.names = T)
  organs = rep("", length(processed.files))
  for(tissue.ctr in 1:length(processed.files))
    organs[tissue.ctr] = str_split(tail(str_split(processed.files[tissue.ctr], "/")[[1]], n=1), "\\.")[[1]][1]
  
  if(!run.celltype)
    sink(paste0("scripts/", data.type, "/run_BASiCS_all_tissues_types_", data.type, ".sh"))
  else    
    sink(paste0("scripts/", data.type, "/run_BASiCS_all_cell_types_", data.type, ".sh"))
  
  for(tissue.ctr in 1:length(processed.files))
  {
    organ = organs[tissue.ctr]
    print(paste0("Write scripts organ: ", organ))
    
    if(!run.celltype)  # Prepare script per tissue: 
    {
      script.file.name = paste0("run_BASiCS_", organ, ".sh") 
      out.file.name = paste0("BASiCS_log_", organ, ".out") 
      cat(paste0("sbatch -o '", out.file.name, "' ", script.file.name, "\n")) # basics_log_", cur.ct, "' run_BASiCS_", ctr, ".sh"))
      # Create the script 
      fileConn<-file(paste0('scripts/', data.type, "/", script.file.name))  # "module load R4", not in script (should be loaded before!)
      # Set job running environment parameters
      writeLines(c("#!/bin/bash", "",  
                   "#SBATCH --time=168:00:00", 
                   "#SBATCH --ntasks=4", 
                   "#SBATCH --mem=48G", 
                   "module load R4", "", 
                   "#!/usr/bin/env Rscript", "", 
                   "module load R4", "", 
                   paste0("Rscript --vanilla ../../run_BASiCS_commandline.R ", data.type, " ", organ)), fileConn)
      close(fileConn)
#      ctr = ctr + 1
    }  else  # : one script per cell type
    {
      SC = readRDS(file = paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")) # Current tissue Seurat object - heavy code 
      ctr = 1
      for( cur.ct in unique(SC$CT))
      {
        script.file.name = paste0("run_BASiCS_", organ, '_', cur.ct, ".sh") 
        out.file.name = paste0("BASiCS_log_", organ, '_', cur.ct, ".out") 
        cat(paste0("sbatch -o '", out.file.name, "' ", script.file.name, "\n")) # basics_log_", cur.ct, "' run_BASiCS_", ctr, ".sh"))
        # Create the script 
        fileConn<-file(paste0('scripts/', data.type, "/", script.file.name))  # "module load R4", not in script (should be loaded before!)
        # Set job running environment parameters
        writeLines(c("#!/bin/bash", "",  
                     "#SBATCH --time=168:00:00", 
                     "#SBATCH --ntasks=4", 
                     "#SBATCH --mem=48G", 
                     "module load R4", "", 
                     "#!/usr/bin/env Rscript", "", 
                     "module load R4", "", 
                     paste0("Rscript --vanilla ../run_BASiCS_commandline.R ", data.type, " ", organ, " ", cur.ct)), fileConn)
        close(fileConn)
        ctr = ctr + 1
      }  # end loop on cell types
    } # end else
  }  # end loop on tissues
  #  if(!run.celltype)
  #  {
  #    sink(paste0("scripts/", data.type, "/run_BASiCS_all_tissues_types_", data.type, ".sh"))
  #    for( organ in organs)
  #    {
  #      script.file.name = paste0("run_BASiCS_", organ, ".sh") 
  #      out.file.name = paste0("BASiCS_log_", organ,".out") 
  #      cat(paste0("sbatch -o '", out.file.name, "' ", script.file.name, "\n")) # basics_log_", cur.ct, "' run_BASiCS_", ctr, ".sh"))
  #    }
  #  } else
  #  {
  #    sink(paste0("scripts/", data.type, "/run_BASiCS_all_cell_types_", data.type, ".sh"))
  #    for( cur.ct in unique(SC$CT))
  #    {
  #      script.file.name = paste0("run_BASiCS_", organ, '_', cur.ct, ".sh") 
  #      out.file.name = paste0("BASiCS_log_", organ, '_', cur.ct, ".out") 
  #      cat(paste0("sbatch -o '", out.file.name, "' ", script.file.name, "\n")) # basics_log_", cur.ct, "' run_BASiCS_", ctr, ".sh"))
  #    }
  #  }
  sink()  
}  

