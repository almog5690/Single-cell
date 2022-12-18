library(Seurat)
library(rhdf5)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)


# Main script for analysis of single-cell gene-expression data
user.name = "Or" #  "Almog"  # Or 
if(user.name == "Almog")
{
  main.dir = "C:/???/Github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
  data.dir = "???"  # Change to directory with all the data
} else # Or 
{
  # Directories 
  main.dir = "C:/Code/Github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
  data.dir = paste0(main.dir, 'Data/TabulaMuris/') # Can change so code and data may be at different locations. Also change per dataset
}

code.dir <- paste0(main.dir, 'Code/')  # src files 

raw.data.dir <- paste0(data.dir, 'Raw/')  # For raw scRNA-seq gene expression files (one per tissue). Format: h5ad (may differ for different datasets) 
processed.data.dir <- paste0(data.dir, '/Processed/')  # For processed scRNA-seq gene expression files (one per tissue), 
                                                       # after running scRNA_seq_preprocess.R. Format: *.rds files. 
                                                       # This directory will also contain one meta.data file for each dataset. 
gene.data.dir <- paste0(main.dir, 'GeneLevelData/')  # For gene features (selection, length ...)

analysis.results.dir <- paste0(data.dir, '/Analysis/')  # For analysis results (one per dataset for each analysis, e.g. mean, overdispersion ..), 
analysis.figure.dir <- paste0(analysis.results.dir, '/Figures/')  # For analysis figures   
basics.dir <- paste0(analysis.results.dir, '/BASiCS/')  # For BASiCS output 


setwd(code.dir)


# Include all needed source files 
source("scRNA_seq_utilities.R")
source("scRNA_seq_preprocess.R")
source("MeanExpressionAnalysis.R") 
source("OverdispersionAnalysis.R")
source("BASiCSAnalysis.R")
source("MeanFigures.R")
source("OverdispersionFigures.R")

# Data types
data.types = c("TM.facs", "TM.droplet", "CR.Rat") 

# Features to examine (currently, only the first two are implemented)
feature.types = c("gene.len", "selection", "GC%", "CpG%", "TATA", "gene.age", "mRNA.half.life")


# What analysis to do: 
preprocess = FALSE
mean.analysis = TRUE
mean.figures = TRUE
var.analysis = TRUE
var.figures = TRUE


data.type = "TM.droplet"

if(preprocess)
{
  preprocess_droplet()
  preprocess_facs(data.type)
}

if(mean.analysis)
{
  DF_cor = mean_expression_droplet(data.type)  # Should be both droplet and facs in the same function
  DF_cor = mean_expression_facs()  # Should be both droplet and facs in the same function
}
  
if(var.analysis)
{
  DF_cor = var.analysis(data.type)
}

# Plot figures: 
if(mean.figures)
{
  DF_cor = draw_mean_figures(data.type, 3) # Choose specific figure (no need data.type. Figure combines droplet+facs)
}

if(mean.figures)
{
  DF_cor = dra_var_figures(data.type)
}


