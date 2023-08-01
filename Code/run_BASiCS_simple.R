# Simple example of running BASiCS
data.type = "TM.facs"  # Blood_SC
organ = "Aorta" # Blood (?) # smallest facs file, should be fastest 
if(!exists(user.name))  # set default values 
  user.name = "Unix" #  "Almog"  # Or # Unix 
if(!exists(data.type))  # set default values 
  data.type = "TM.facs"  # Choose one type for analysis (can later loop over multiple datasets)
source("scRNA_seq_config.R")
print("Running BASiCS:")
BASiCS_analysis_tissue(data.type, organ)
print("Finished Running BASiCS!")
