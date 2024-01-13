#!/usr/bin/env Rscript
# Start with script with hard-coded arguments. Then later run this with command line arguments 

# New: read from the command line !!! 
args = commandArgs(trailingOnly=TRUE)

source("scRNA_seq_config.R")

# What analysis to do: 
for.paper = FALSE  # which analysis to do 
preprocess = FALSE
reg.analysis = TRUE
reg.figures = TRUE
# var.analysis = TRUE
# var.figures = TRUE


if(for.paper == FALSE)  # Run one regression : choose the response, the expression covariates, and the gene features 
{
  data.types = c("TM.facs")  # Data type: Tabule-Muris facs 
  feature.types = c("alpha.missense")  # gene feature: alpha-missnese gene conservation 
  expression.stat.y = c("mge")  # response: over-dispersion
  expression.stat.x = c("") # expression covariate: mean gene expression 
}


if(preprocess)
{
  preprocess_droplet()
  preprocess_facs(data.type)
}


# Possible Different figure types: 1.young+old log-pval barplot, 2. young vs. old new-scatter, 3. fc log-pval barplot (same as 1?)
# (4: Figure 5 of paper: BASicS results, overdispersion)


#data.type= "TM.Droplet"
#TRY.TM.DROPLET = expression_regression_analysis("TM.droplet", feature.types = paper.features[[1]], 
#                               expression.stat.x = paper.expression.stat.x[i], expression.stat.y = paper.expression.stat.y[1], 
#                               force.rerun = FALSE)
#
#data.type = "Blood_SC"
#TRY.HUMAN.BLOOD = expression_regression_analysis("Blood_SC", feature.types = paper.features[[1]], 
#                                                expression.stat.x = paper.expression.stat.x[1], expression.stat.y = paper.expression.stat.y[1], 
#                                                force.rerun = FALSE)
#
#draw_expr_reg_figures(c("Blood_SC"), expr.stat.y = "mean", expr.stat.x = c(), 
#                      2, feature.types = "gene.len") # Choose specific figure (no need data.type. Figure combines droplet+facs)
#
#draw_expr_reg_figures(c("TM.Droplet"), expr.stat.y = "mean", expr.stat.x = c(), 
#                      2, feature.types = "gene.len") # Choose specific figure (no need data.type. Figure combines droplet+facs)



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
  
} else # here do ONE custom analysis 
{
  if(reg.analysis)  # regression with different covariates and response expression variables 
  {
    #  feature.types = c("gene.len", "selection")  # simplify to get a regression analysis similar to correlation analysis!!!
    for(data.type in data.types)
    {
      DF_cor_mean = expression_regression_analysis(data.type, feature.types = feature.types, 
                                                   expression.stat.x = expression.stat.x, 
                                                   expression.stat.y = expression.stat.y, 
                                                   force.rerun = FALSE) 
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


