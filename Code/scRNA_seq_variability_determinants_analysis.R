library(Seurat)
library(rhdf5)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(colorRamps)
library(ggpointdensity)
library(ggpubr)

# Data types
data.types = c("TM.facs", "TM.droplet", "CR.Rat", "Age.Anno") 
data.dirs = c("TabulaMuris", "TabulaMuris", "RatCR", "HumanAgeAnno")   # Names of directories 
names(data.dirs) = data.types 
processed.files.str <- c("facs", "drop", "rat", "anno") # match these (then rds)
names(processed.files.str) = data.types
data.type = "TM.droplet"  # Choose one type for analysis (can later loop over multiple datasets)


# Main script for analysis of single-cell gene-expression data
user.name = "Or" #  "Almog"  # Or 
if(user.name == "Almog")
{
  main.dir = "C:/Users/User/OneDrive/Documents/Github/Single-cell/"  # Change to your local path. This path should be used everywhere as the main source path
  main.data.dir = "C:/Users/User/Google Drive/SingleCell/"  # Change to directory with all the data
} else # Or 
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
age.groups <- c("all", "young", "old")


# Set filtering parameters for filtering (may depend on the data type?): 
filter.params = c()
filter.params$min.count <- 10 
filter.params$min.cells.total <- 100
filter.params$min.cells.per.age <- 20
filter.params$min.genes.for.reg <- 50
filter.params$pseudo.count <- 1 # for regularization


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

# What analysis to do: 
for.paper = TRUE  # which analysis to do 
preprocess = FALSE
reg.analysis = TRUE
reg.figures = TRUE
# var.analysis = TRUE
# var.figures = TRUE

if(preprocess)
{
  preprocess_droplet()
  preprocess_facs(data.type)
}


# Possible Different figure types: 1.young+old log-pval barplot, 2. young vs. old new-scatter, 3. fc log-pval barplot (same as 1?)
# (4: Figure 5 of paper: BASicS results, overdispersion)

# Actual regression runs for the paper: 
paper.expression.stat.y <- c("mean", "mean", "overdispersion", "overdispersion", "overdispersion")
paper.expression.stat.x <- c("",     "",    "mean", "mean", "mean")
paper.features <- list("gene.len", "selection", "gene.len", "selection", c("selection", "gene.len")) # last one of two feature covariates is for supp. info.


if(for.paper)
{
  ctr <- 1
  n.paper.reg <- length(paper.expression.stat.y)
  paper.DF <- vector("list", length(data.types))
  if(reg.analysis)  # regression with different covariates and response expression variables 
    for(data.type in data.types)
    {
      paper.DF[[ctr]] <- vector("list", n.paper.reg)
      for(i in 1:n.paper.reg)  # run different regression types 
        paper.DF[[ctr]][[i]] = expression_regression_analysis(data.type, feature.types = paper.features[[i]], 
                                                     expression.stat.x = paper.expression.stat.x[i], expression.stat.y = paper.expression.stat.y[i], 
                                                     force.rerun = FALSE) #  c("selection", "gene.len"))  # Should be both droplet and facs in the same function
      ctr <- ctr + 1
    }
  
  if(reg.figures) # plot figures for paper 
  {
    for (fig.num in c(1,2,11,111,99))     # for each option, save the 3 figures (some are used in main or supp. info of the paper)
      for(data.type in data.types)
        draw_expr_reg_figures(data.type, expr.stat.y = "mean", expr.stat.x = c(),  # each data type separately 
                            fig.num, feature.types = feature.types) # Choose specific figure (no need data.type. Figure combines droplet+facs)
    
  }
  
} else # here do custom analysis 
{
  
  
  
  if(reg.analysis)  # regression with different covariates and response expression variables 
  {
    #  feature.types = c("gene.len", "selection")  # simplify to get a regression analysis similar to correlation analysis!!!
    for(data.type in data.types)
    {
      DF_cor_mean = expression_regression_analysis(data.type, feature.types = feature.types, force.rerun = FALSE) #  c("selection", "gene.len"))  # Should be both droplet and facs in the same function
      #    DF_cor_overdispersion = expression_regression_analysis(data.type, feature.types = feature.types, 
      #                                                           expression.stat.y = c("overdispersion"), expression.stat.x = c("mean"), # set y and covariates  
      #                                                           force.rerun = FALSE) #  c("selection", "gene.len"))  # Should be both droplet and facs in the same function
    }
  }
  
  # Plot figures: 
  if(reg.figures)
  {
    for (fig.num in c(1,2,11,111,99))
      draw_expr_reg_figures(c("TM.facs", "TM.droplet"), expr.stat.y = "mean", expr.stat.x = c(), 
                            fig.num, feature.types = feature.types) # Choose specific figure (no need data.type. Figure combines droplet+facs)
    #  draw_mean_figures(c("CR.Rat"), 1, feature.types = c("selection", "gene.len")) 
    for (fig.num in c(1,2,11,111))
      draw_expr_reg_figures(c("CR.Rat"), expr.stat.y = "mean", expr.stat.x = c(),
                            fig.num, feature.types = feature.types, # c("selection", "gene.len"), 
                            tissue = "Liver", cell_type = "7") # need different cell-type example for Rat !! 
    draw_expr_reg_figures(c("TM.facs", "TM.droplet", "CR.Rat"), expr.stat.y = "mean", expr.stat.x = c(), 66, feature.types = feature.types, n.features=5)
    draw_expr_reg_figures(c("TM.facs", "TM.droplet", "CR.Rat"), expr.stat.y = "overdispersion", expr.stat.x = c("mean"), 666, feature.types = feature.types)
  }
}  # end if paper analysis 


