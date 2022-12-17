library(Seurat)
library(rhdf5)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)


# Main script for analysis of single-cell gene-expression data

# Directories 
main.dir = "C:/Code/Github/Single-cell"  # Change to your local path. This path should be used everywhere
code.dir = paste0(main.dir, '/Code')  # src files 
raw.data.dir = paste0(main.dir, '/Data')  # For raw scRNA-seq gene expression files 
processed.data.dir = paste0(main.dir, '/Data/Processed')  # For processed scRNA-seq gene expression files 
gene.data.dir = paste0(main.dir, '/GeneLevelData')  # For gene features 
setwd(code.dir)


# Include all needed source files 
source("scRNA_seq_utilities.R")
source("MeanExpressionAnalysis.R") 
source("OverdispersionAnalysis.R")
source("BASiCSAnalysis.R")
source("MeanFigures.R")
source("OverdispersionFigures.R")

# Data types
data.types = c("TM.facs", "TM.droplet", "CR.RAT") 


# What analysis to do: 
preprocess = FALSE
mean.analysis = TRUE
mean.figures = TRUE
var.analysis = TRUE
var.figures = TRUE


if(preprocess)
{
  preprocess.scRNAseq(data.type)
}

if(mean.analysis)
{
  DF_cor = mean.analysis(data.type)
}
  
if(var.analysis)
{
  DF_cor = var.analysis(data.type)
}

# Plot figures: 
if(mean.figures)
{
  DF_cor = draw_mean_figures(data.type)
}

if(mean.figures)
{
  DF_cor = dra_var_figures(data.type)
}


