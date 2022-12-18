library(Seurat)
library(rhdf5)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)


# Main script for analysis of single-cell gene-expression data

# Directories 
main.dir = "C:/Code/Github/Single-cell/"  # Change to your local path. This path should be used everywhere
data.dir = paste0(main.dir, 'Data/TabulaMuris/') # Can change so code and data may be at different locations. Also change per dataset
code.dir <- paste0(main.dir, 'Code/')  # src files 
raw.data.dir <- paste0(data.dir, 'Raw/')  # For raw scRNA-seq gene expression files 
processed.data.dir <- paste0(data.dir, '/Processed/')  # For processed scRNA-seq gene expression files 
gene.data.dir <- paste0(main.dir, 'GeneLevelData/')  # For gene features 
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
  DF_cor = draw_mean_figures(data.type)
}

if(mean.figures)
{
  DF_cor = dra_var_figures(data.type)
}


