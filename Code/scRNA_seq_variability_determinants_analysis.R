library(Seurat)
library(rhdf5)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)
library(ggpubr)


# Data types
data.types = c("TM.facs", "TM.droplet", "CR.Rat") 
data.dirs = c("TabulaMuris", "TabulaMuris", "RatCR")   # Names of directories 
names(data.dirs) = data.types 

processed.files.str <- c("facs", "drop", "rat") # match these (then rds)
names(processed.files.str) = data.types

data.type = "TM.droplet"  # Choose one type for analysis (can later loop over multiple datasets)


# Main script for analysis of single-cell gene-expression data
user.name = "Or" #  "Almog"  # Or 
if(user.name == "Almog")
{
  main.dir = "C:/Users/User/OneDrive/Documents/Github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
  main.data.dir = "C:/Users/almog/Singel cell"  # Change to directory with all the data
} else # Or 
{
  # Directories 
  main.dir = "C:/Code/Github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
  main.data.dir = "G:/.shortcut-targets-by-id/1Kq8HX8zZy_lm6GvFIGkmL5c42sEU3VBx/SingleCell/" 
  #   paste0(main.dir, 'Data/', data.dirs[data.type], '/') # Can change so code and data may be at different locations. Also change per dataset
}
code.dir <- paste0(main.dir, 'Code/')  # src files 

gene.data.dir <- paste0(main.data.dir, 'GeneLevelData/')  

#data.dir <- paste0(main.data.dir, 'Data/', data.dirs[data.type], '/')
#raw.data.dir <- paste0(data.dir, 'Raw/')  # For raw scRNA-seq gene expression files (one per tissue). Format: h5ad (may differ for different datasets) 
#processed.data.dir <- paste0(data.dir, 'Processed/')  # For processed scRNA-seq gene expression files (one per tissue), 
#                                                       # after running scRNA_seq_preprocess.R. Format: *.rds files. 
#                                                       # This directory will also contain one meta.data file for each dataset. 
#analysis.results.dir <- paste0(data.dir, 'Analysis/')  # For analysis results (one per dataset for each analysis, e.g. mean, overdispersion ..), 
#analysis.figures.dir <- paste0(analysis.results.dir, 'Figures/')  # For analysis figures   
#basics.dir <- paste0(analysis.results.dir, 'BASiCS/')  # For BASiCS output 

setwd(code.dir)


# Set analysis types we want to do: 
y.test <- c("mge", "od", "od", "od", "Delta.mge", "Delta.od", "Delta.od") # The dependent variable 
x.test <- c("", "", "", "mge", "mge", "mge", "mge") # The covariates 
old.yound.delta <- c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) # Test also difference between young and old? 
cor.test <- c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)  # Compute also correlation coefficient?


# Set filtering parameters for filtering (may depend on the data type?): 
filter.params = c()
filter.params$min.count <- 10 
filter.params$min.cells.total <- 100
filter.params$min.cells.per.age <- 20


# Include all needed source files 
source("scRNA_seq_utilities.R")
source("scRNA_seq_preprocess.R")
source("MeanExpressionAnalysis.R") 
source("OverdispersionAnalysis.R")
source("BASiCSAnalysis.R")
source("MeanFigures.R")
source("OverdispersionFigures.R")


# Features to examine (currently, only the first two are implemented)
feature.types = c("gene.len", "selection", "TATA", "mRNA.half.life", "GC", "CpG")  # These aren't ready yet: "GC", "CpG%",  "gene.age", )

# What analysis to do: 
preprocess = FALSE
mean.analysis = TRUE
mean.figures = TRUE
var.analysis = TRUE
var.figures = TRUE

if(preprocess)
{
  preprocess_droplet()
  preprocess_facs(data.type)
}

if(mean.analysis)
{
  DF_cor.drop = mean_expression_analysis("TM.droplet", feature.types =feature.types) #  c("selection", "gene.len"))  # Should be both droplet and facs in the same function
  DF_cor.facs = mean_expression_analysis("TM.facs", feature.types = c("selection", "gene.len"))  # Should be both droplet and facs in the same function
  DF_cor.rat = mean_expression_analysis("CR.Rat", feature.types = c("selection", "gene.len"))  
}
  
if(var.analysis)
{
  DF_cor = var_analysis(data.type)
}

# Plot figures: 
if(mean.figures)
{
  draw_mean_figures(c("TM.facs", "TM.droplet"), 2) # Choose specific figure (no need data.type. Figure combines droplet+facs)
}

if(mean.figures)
{
  DF_cor = dra_var_figures(data.type)
}


