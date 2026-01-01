library(Seurat)
library(rhdf5)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)
library(ggpubr)


# =============================================================================
# ENVIRONMENT CONFIGURATION - Change this single variable to switch environments
# =============================================================================
# Options:
#   "or.windows"  - Or Zuk, Windows local
#   "or.wsl"      - Or Zuk, WSL (Windows Subsystem for Linux)
#   "or.cluster"  - Or Zuk, Unix cluster
#   "almog.windows" - Almog Yair, Windows local
# =============================================================================
if(!exists('run.env'))
  run.env <- "or.windows"  # <-- CHANGE THIS TO SWITCH ENVIRONMENT

if(!exists('data.type'))
  data.type <- "TM.droplet"

# Data types (shared across all environments)
data.types <- c("TM.facs", "TM.droplet", "CR.Rat", "Age.Anno", "Blood_SC", "MCA")
data.dirs <- c("TabulaMuris", "TabulaMuris", "RatCR", "HumanAgeAnno", "HumanBlood", "MCA")
names(data.dirs) <- data.types
processed.files.str <- c("facs", "drop", "rat", "anno", "SC", "MCA")
names(processed.files.str) <- data.types

# Set paths based on environment
if(run.env == "or.windows") {
  main.dir <- "C:/Code/Github/Single-cell/"
  main.data.dir <- "G:/.shortcut-targets-by-id/1Kq8HX8zZy_lm6GvFIGkmL5c42sEU3VBx/SingleCell/"
} else if(run.env == "or.wsl") {
  main.dir <- "/mnt/c/Code/Github/Single-cell/"
  main.data.dir <- "/mnt/c/Code/Github/Single-cell/"  # Data in same repo for WSL
} else if(run.env == "or.cluster") {
  main.dir <- "/sci/labs/orzuk/orzuk/github/Single-cell/"
  main.data.dir <- "/sci/labs/orzuk/orzuk/projects/SingleCell/"
} else if(run.env == "almog.windows") {
  main.dir <- "C:/Users/User/OneDrive/Documents/Github/Single-cell/"
  main.data.dir <- "C:/Users/User/Google Drive/SingleCell/"
} else {
  stop(paste("Unknown run.env:", run.env,
             "\nValid options: or.windows, or.wsl, or.cluster, almog.windows"))
}

# Derived paths (computed from main.dir and main.data.dir)
code.dir <- paste0(main.dir, 'Code/')
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


# Include all needed source files (utilities already sourced above for set_data_dirs)
source("scRNA_seq_preprocess.R")
source("ExpressionRegressionAnalysis.R")  # Handles both mean and overdispersion analysis
source("BASiCSAnalysis.R")
source("MeanFigures.R")
source("OverdispersionFigures.R")

# Features to examine (currently, only the first two are implemented)
feature.types = c("gene.len", "selection", "TATA", "mRNA.half.life", "GC", "CpG")  # These aren't ready yet: "GC", "CpG%",  "gene.age", )

# Actual regression runs for the paper:  (5 models + additional five for fold-change?)
paper.expression.stat.y <- c("mean", "mean", "overdispersion", "overdispersion") # , "overdispersion")  # , "overdispersion")
paper.expression.stat.x <- c("", "", "mean", "mean") # , "mean") # , "age")  # age is special, should use a different code in expression regression analysis!
paper.features <- list("gene.len", "selection", "gene.len", "selection") # , c("selection", "gene.len")) # last one of two feature covariates is for supp. info.

# First: is correlation positive or negative, second: is it stronger in young or in old. 
# For example: neg old means that we expect negative correlations, and that the correlations will be larger in absolute value in old
paper.hypothesis <- c("neg,old", "pos,young", "pos(?),old(?)", "neg,young")  
paper.fc.hypothesis <- c("neg", "neg", "neg", "pos")  # is the correlation of old-young vs. feature positive or negative
