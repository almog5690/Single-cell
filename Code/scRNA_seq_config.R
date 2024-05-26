library(Seurat)
library(rhdf5)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)
library(ggpubr)


if(!exists('user.name'))  # set default values 
  user.name = "Or" #  "Almog"  # Or # Unix 
if(!exists('data.type'))  # set default values 
  data.type = "TM.droplet"  # Choose one type for analysis (can later loop over multiple datasets)

# Data types
data.types = c("TM.facs", "TM.droplet", "CR.Rat", "Age.Anno", "Blood_SC", "MCA") 
data.dirs = c("TabulaMuris", "TabulaMuris", "RatCR", "HumanAgeAnno", "HumanBlood", "MCA")   # Names of directories 
names(data.dirs) = data.types 
processed.files.str <- c("facs", "drop", "rat", "anno", "SC", "MCA") # match these (then rds)
names(processed.files.str) = data.types




# Main script for analysis of single-cell gene-expression data
if(user.name == "Almog")
{
  main.dir = "C:/Users/User/OneDrive/Documents/Github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
  main.data.dir = "C:/Users/User/Google Drive/SingleCell/"  # Change to directory with all the data
} else # Or 
if(user.name == "Or")
{
#  G:\.shortcut-targets-by-id\1Kq8HX8zZy_lm6GvFIGkmL5c42sEU3VBx\SingleCell   
#  G:\.shortcut-targets-by-id\1Kq8HX8zZy_lm6GvFIGkmL5c42sEU3VBx\SingleCell\GeneLevelData  
  # Directories 
  main.dir = "C:/Code/Github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
  #  main.data.dir = "G:/.shortcut-targets-by-id/1Kq8HX8zZy_lm6GvFIGkmL5c42sEU3VBx/SingleCell/" 
#  main.data.dir = "G:/My Drive/Students/AlmogYair/SingleCell/"
  main.data.dir = "G:/.shortcut-targets-by-id/1Kq8HX8zZy_lm6GvFIGkmL5c42sEU3VBx/SingleCell/"  # give shortcut (problem: may change!)
  #   paste0(main.dir, 'Data/', data.dirs[data.type], '/') # Can change so code and data may be at different locations. Also change per dataset
}
if(user.name == "Unix")
{
  main.dir = "/sci/labs/orzuk/orzuk/github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
  main.data.dir = "/sci/labs/orzuk/orzuk/projects/SingleCell/"  # Change to directory with all the data
}
code.dir <- paste0(main.dir, 'Code/')  # src files 
gene.data.dir <- paste0(main.data.dir, 'GeneLevelData/')  
res.dir <- paste0(main.data.dir, "Results")


#data.dir <- paste0(main.data.dir, 'Data/', data.dirs[data.type], '/')
#raw.data.dir <- paste0(data.dir, 'Raw/')  # For raw scRNA-seq gene expression files (one per tissue). Format: h5ad (may differ for different datasets) 
#processed.data.dir <- paste0(data.dir, 'Processed/')  # For processed scRNA-seq gene expression files (one per tissue), 
#                                                       # after running scRNA_seq_preprocess.R. Format: *.rds files. 
#                                                       # This directory will also contain one meta.data file for each dataset. 
#analysis.results.dir <- paste0(data.dir, 'Analysis/')  # For analysis results (one per dataset for each analysis, e.g. mean, overdispersion ..), 
#analysis.figures.dir <- paste0(analysis.results.dir, 'Figures/')  # For analysis figures   
#basics.dir <- paste0(analysis.results.dir, 'BASiCS/')  # For BASiCS output 

setwd(code.dir)
source("scRNA_seq_utilities.R")

set_data_dirs(data.type)  # set all directories for this data type


# Set analysis types we want to do: 
y.test <- c("mge", "od", "od", "od", "Delta.mge", "Delta.od", "Delta.od") # The dependent variable 
x.test <- c("", "", "", "mge", "mge", "mge", "mge") # The covariates 
old.yound.delta <- c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) # Test also difference between young and old? 
cor.test <- c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)  # Compute also correlation coefficient?
age.groups <- c("all", "young", "old")


# Set filtering parameters for filtering (may depend on the data type?): 
filter.params = c()
filter.params$min.count <- 10 
filter.params$min.cells.total <- 100
filter.params$min.cells.per.age <- 20
filter.params$min.genes.for.reg <- 50
filter.params$pseudo.count <- 1 # for regularization
filter.params$expression.thresh = 0.2


# Include all needed source files 
source("scRNA_seq_utilities.R")
source("scRNA_seq_preprocess.R")
source("ExpressionRegressionAnalysis.R") 
source("OverdispersionAnalysis.R")
source("BASiCSAnalysis.R")
source("MeanFigures.R")
source("OverdispersionFigures.R")

# Features to examine (currently, only the first two are implemented)
feature.types = c("gene.len", "selection", "TATA", "mRNA.half.life", "GC", "CpG")  # These aren't ready yet: "GC", "CpG%",  "gene.age", )

# Actual regression runs for the paper:  (5 models + additional five for fold-change?)
paper.expression.stat.y <- c("mean", "mean", "overdispersion", "overdispersion") # , "overdispersion")  # , "overdispersion")
paper.expression.stat.x <- c("",     "", "mean", "mean") # , "mean") # , "age")  # age is special, should use a different code in expression regression analysis!
paper.features <- list("gene.len", "selection", "gene.len", "selection") # , c("selection", "gene.len")) # last one of two feature covariates is for supp. info.

