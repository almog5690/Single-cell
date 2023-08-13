source("scRNA_seq_config.R")

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



TRY.TM.DROPLET = expression_regression_analysis(data.type, feature.types = paper.features[[i]], 
                               expression.stat.x = paper.expression.stat.x[i], expression.stat.y = paper.expression.stat.y[i], 
                               force.rerun = TRUE)

data.type = "Blood_SC"
TRY.HUMAN.BLOOD = expression_regression_analysis(data.type, feature.types = paper.features[[i]], 
                                                expression.stat.x = paper.expression.stat.x[i], expression.stat.y = paper.expression.stat.y[i], 
                                                force.rerun = TRUE)

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


