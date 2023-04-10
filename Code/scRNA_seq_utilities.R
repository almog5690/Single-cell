# Utilities needed for single cell analysis and visualization 
library(readxl)


# Set directories for raw and processed data, gene features data, analysis results, figures .. for a given data type 
set_data_dirs <- function(data.type)
{
  # Set globally: 
  data.dir <<- paste0(main.data.dir, 'Data/', data.dirs[data.type], '/')
  raw.data.dir <<- paste0(data.dir, 'Raw/')  # For raw scRNA-seq gene expression files (one per tissue). Format: h5ad (may differ for different datasets) 
  processed.data.dir <<- paste0(data.dir, 'Processed/')  # For processed scRNA-seq gene expression files (one per tissue), 
  # after running scRNA_seq_preprocess.R. Format: *.rds files. 
  # This directory will also contain one meta.data file for each dataset. 
  analysis.results.dir <<- paste0(data.dir, 'Analysis/')  # For analysis results (one per dataset for each analysis, e.g. mean, overdispersion ..), 
  analysis.figures.dir <<- paste0(analysis.results.dir, 'Figures/')  # For analysis figures   
  basics.dir <<- paste0(data.dir, 'BASiCS/')  # For BASiCS output 
  
  # Return all in a list? or just keep them as updated global variables  
}

# Get tissues. Use the processed files if they exist (otherwise need to go back to raw files?)
get_tissue_file_names <- function(data.type)
{
  # loading the meta data names for the droplet technology
  if(data.type == "TM.droplet")
  {
    load(paste0(processed.data.dir, "/meta.data.drop.RData"))
    file.patt <- "drop.r"
    file.names = list.files(path = processed.data.dir, pattern = "drop.rds") # list of tissue from droplet technology
    #removing large intestine and trachea because they have only one age group and also fat,skin and pancreas.
    file.names = file.names[-grep(c("Fat|Large_Intestine|Skin|Pancreas|Trachea"), file.names)]
    organs = unname(sapply(file.names,function(f) unlist(strsplit(f,"[.]"))[1]))
  }
  if(data.type == "TM.facs")
  {
    load(paste0(processed.data.dir, "/meta.data.facs.RData"))
    file.patt <- "facs"
    file.names = list.files(path = processed.data.dir, pattern = "*facs.rds") # facs processed files list
    organs  = unname(sapply(file.names,function(f) unlist(strsplit(f,"[.]"))[1])) # facs organs names
    #    organs[5] = "Brain_Non-Myeloid" # why needed?
  }
  if(data.type == "CR.Rat")  # to fill 
  {
    organs = c("Aorta","BAT","BM","Brain","Muscle","Skin","WAT","Kidney","Liver")
    file.names = paste(organs,"rds", sep = ".")
  }
  return(list(file.names = file.names, organs = organs))
}


# Filter cells from an expression matrix (not used yet. Should be part of analysis)
filter_cells <- function(cell_types, young.ind, old.ind, filter.params)
{
  n_cell_types = max(cell_types) # Number of cell types
  erase = vector()
  
  for(ct_ind in 0:n_cell_types){ # filtering cell types with less then 100 cells or less the 20 cells in each the age groups
    if(sum(cell_types==ct_ind) < filter.params$min.cells.total |
       sum(cell_types==ct_ind & young.ind) < filter.params$min.cells.per.age |
       sum(cell_types==ct_ind & old.ind) < filter.params$min.cells.per.age){
      erase = c(erase, ct_ind+1)
      next()
    }
  }
  cells_ind = c(1:(n_cell_types+1))[-erase]
  if(length(cells_ind) == 0) cells_ind = (1:(n_cell_types+1))
  return(cells_ind) # indices of cells that we keep 
}


# Read feature files (not used yet. NEXT ONE!!! )
# Output will be a list with names being the genes, and values being numerical values.
# Gene names are according to Ensembl (???)
# organism - allow multiple ones in the future (currently only mice supported)
read_gene_features  <- function(feature.names, organism = "mice", force.rerun = FALSE)
{
  gene.features.outfile <- paste0(main.data.dir, 'Data/', 'gene.features.', 
                                  paste0( feature.names, collapse="_"), '.RData')
  if(file.exists(gene.features.outfile) & (force.rerun==FALSE))
  {
    load(gene.features.outfile)
    return(gv)
  }
  
  n.features <- length(feature.names)
  gv <- vector("list", n.features)
  names(gv) <- feature.names
  for(feature.name in feature.names)
  {
    if(feature.name == "selection")
    {
      select <- read.delim(paste0(gene.data.dir, "gnomad.v2.1.1.lof_metrics.by_gene.txt"))
      gene.names = toupper(select$gene) # making all the gene names in to upper case letters
      gene.values = select$pLI # the selection score vector
      names(gene.values) = gene.names # naming each score with the gene it belongs to
    }
    if(feature.name %in% c("gene.len", "gene.length", "genes.len"))
    {
      # reading the transcript length data
      length_data = read.delim(paste0(gene.data.dir, "mart_export.txt"))
      gene.values = length_data[,7] # gene length vector
      names(gene.values) = toupper(length_data[,6]) # Add gene names. Note: there are many values per gene 
    }
    if(feature.name %in% c("TATA", "GC", "CpG"))
    {
      if(!exists('tata.data')) # read only once 
        tata.data <- read.delim(paste0(gene.data.dir, "mm10_tata_annotation_v3.3.tsv"))
      if(feature.name == "TATA")# TATA BOX
        gene.values <- as.integer(unlist(lapply(tata.data$TATA.Box.TBP..Promoter.Homer.Distance.From.Peak.sequence.strand.conservation., nchar)) > 0)
      if(feature.name == "GC")
        gene.values <- tata.data$GC.
      if(feature.name == "CpG")
        gene.values <- tata.data$CpG.
      names(gene.values) <- toupper(tata.data$Gene.Name)
    }
    if(feature.name == "mRNA.half.life")
    {
      half.life.data <-  read_excel(paste0(gene.data.dir, "13059_2022_2811_MOESM3_ESM.xlsx"), sheet = "mouse", skip =2)
      gene.values <- half.life.data$Half_life_PC1
      names(gene.values) <- toupper(half.life.data$Gene.name)
    }
    if(feature.name == "gene.age")  # extract gene age (T.B.D.)
    {
      
    }
    
    gene.values = aggregate(x = gene.values, by = list(names(gene.values)), FUN = mean)  # perform unique
    gv[[feature.name]] <- gene.values$x
    names(gv[[feature.name]]) <- gene.values$Group.1
  }  
  save(gv, file=gene.features.outfile)   # Save dataframe to file ( Save also to excel? ) 
  return(gv)
}


# Get union of items in list 
list_to_common_dataframe <- function(l)
{
  n <- length(l)
  row.names <- c()  
  for( i in 1:n)
    row.names <- union(row.names, names(l[[i]]))
  n.row <- length(row.names)
  df_common <- data.frame(matrix(NA, ncol = n, nrow = n.row))
  colnames(df_common) <- names(l)
  rownames(df_common) <- row.names
  for( i in 1:n)
    df_common[names(l[[i]]), names(l)[i]] <- l[[i]]
  return(df_common)
}

# Compute features of expression (mean, overdispersion, variance ... )
# for a given tissue/cell type
# organ - what tissues
# cell.types - which ones. Default: all types in a given tissue
# expression.stats - which expression statistics to extract 
# age.groups - compute for each group separately
# Need both organ and Seurat output (it doesn't contain the organ/tissue)
# Take list of data frames for all cell types  
extract_expression_statistics <- function(data.type, organ, cell.types=c(), expression.stats = c("mean", "overdispersion"), 
                                          age.groups = c("young", "old", "all", "fc"), 
                                          SeuratOutput=c(), BASiCSOutput = c(), force.rerun = FALSE)
{
  set_data_dirs(data.type)
  expression.statistics.outfile <- paste0(analysis.results.dir, 'expression.stats.', data.type, '.', organ, '.',
                                          paste0( expression.stats, collapse="_"), '.RData')
  samples <- get_tissue_file_names(data.type)
  organ.ind = which(samples$organs == organ)
  meta.data = get_meta_data(data.type)
  groups = dataset_to_age_groups(data.type)
  n.stats = length(expression.stats)
  
  if(file.exists(expression.statistics.outfile) & (force.rerun==FALSE))
  {
#    print("Loading file!")
#    print(expression.statistics.outfile)
    load(expression.statistics.outfile)
    return(list(DF.expr.stats=DF.expr.stats, cell_types_categories=cell_types_categories, cells_ind=cells_ind))
  }
  DF.expr.stats <- list() # list of data-frames, one for each cell type 
  
  # get statistics names 
  stats.col.names <- c()
  for(stat in expression.stats)
    for(age.group in age.groups)  # loop on age groups
      stats.col.names <- c(stats.col.names, paste0(stat, "_", age.group))
  stats.col.names <- c(stats.col.names, "filter") # names of genes filtered 
  if("overdispersion" %in% expression.stats) # For overdispersion read DVT files if they exist/run basics  
  {
#    stats.col.names <- c(stats.col.names, "overdispersion_fc")
    stats.col.names <- c(stats.col.names, "overdispersion_diff")
  }
#  print("Col names:")
#  print(stats.col.names)

#  print("Load tissue, file name:")
#  print( paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")  )
  if(length(SeuratOutput)==0) # empty, on first time the loop runs 
    SeuratOutput = readRDS(file = paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")) # Current tissue Seurat object
  list2env(tissue_to_age_inds(data.type, organ, groups, SeuratOutput@meta.data), env=environment()) # set specific ages for all age groups in all datasets
  counts.mat = as.matrix(SeuratOutput@assays$RNA@data) # the data matrix for the current tissue
  n.genes <- dim(counts.mat)[1]  # same number of genes for all cell types (could have many NAs)
  n.cells <- dim(counts.mat)[2]
  n.groups <- length(age.groups)
  
  if(data.type == "CR.Rat"){ # Convert to dummy variables 
    cell_types = as.numeric(SeuratOutput@meta.data$cell_types) # Cell types vector
    cell_types_categories = levels(SeuratOutput$cell_types)
  } else
  {
    cell_types = SeuratOutput@meta.data$cell.ontology.class # Cell types vector
    cell_types_categories = meta.data[[organ.ind]]$cell_ontology_class # Cell type names. Missing variable meta.data.drop
  }
  cells_ind = filter_cells(cell_types, young.ind, old.ind, filter.params)  #  unique(cell_types)+1 # NO FILTERING NOW !!!
  n.cell.types <- length(cells_ind) # number of cell types in tissue
  
  for(cell.type in cells_ind)  # First load data if not loaded already 
  {
    DF.expr.stats[[cell.type]] <- as.data.frame(matrix(NA ,nrow = n.genes, ncol = length(stats.col.names))) # n.stats*n.groups+1) # Set negatives 
    colnames(DF.expr.stats[[cell.type]]) <- stats.col.names
    rownames(DF.expr.stats[[cell.type]]) <- toupper(rownames(SeuratOutput@assays$RNA@data))

    if(("overdispersion" %in% expression.stats) & (length(BASiCSOutput) == 0)) # For overdispersion read DVT files if they exist/run basics  
    {
      BASiCS.files <- dataset_to_BASiCS_file_names(data.type, organ, cell_types_categories[cell.type]) # Load BASiCS results
      print(paste0("Load BASiCS file: ", (BASiCS.files[[1]])))
      
      if(!file.exists(BASiCS.files$test)) # try renaming (for Rats)
        BASiCS.files$test <- str_replace(BASiCS.files$test, "DVT ", "DVT test ")
      # Don't let corrupt files kill the run !! 
      read.flag <- tryCatch( {
        load(BASiCS.files$test)         # ,  temp_env <- new.env())  # one variable called 'test' #      env_list <- as.list(temp_env)
        },
        error = function(cond){return (FALSE)})
      if(read.flag == FALSE)
      {
        print("Error! couldn't load!")
        next  # move to next cell type
      }
      
      # Extract overdispersion features !!! 
      df.od = test@Results$Disp@Table
      df.od$GeneName = toupper(df.od$GeneName) # Uppercasing gene names
      df.od = df.od[!df.od$ResultDiffDisp %in% c("ExcludedFromTesting"),] # Filtering genes significant difference in mean expression
    }
    
    for(age.group in age.groups)  # loop on age groups
    {
      cur.cell.ind = switch(age.group, # Indices of cells in each group
                       "all" = all.ind, 
                       "young" = young.ind, 
                       "old" = old.ind, 
                       "fc" = young.ind)  & (cell_types==cell.type-1) # take only cell type 
      # Filter genes with low expression for old/young (less than 10 counts). Keep indices of genes FILTER BY ALL!!!
      if(age.group == "fc")  #union young or old 
        cur.gene.ind = rowSums(SeuratOutput@assays$RNA@counts[, young.ind]) > filter.params$min.count |
        (rowSums(SeuratOutput@assays$RNA@counts[, old.ind]) > filter.params$min.count )
      else
        cur.gene.ind = rowSums(SeuratOutput@assays$RNA@counts[, all.ind]) > filter.params$min.count # cur.cell.ind
      
      names(cur.gene.ind) = toupper(names(cur.gene.ind))

      # Next, extract different statistics 
      for(stat in expression.stats)
      {
        if(stat == "mean")
        {
          if(age.group == "fc")  # union young or old 
          {
            bulk.data.types <- c()  # add bulk data types with fold change information
            if(data.type %in% bulk.data.types)
              DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_fc")] <- c() # read excel TBD
            else
              DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_fc")] <- log( # add regularization
                (DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_old")] + filter.params$pseudo.count) / 
                  (DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_young")] + filter.params$pseudo.count)) # Take only filtered cells. why? 
          } else
            # Repeat for all age groups and features 
            DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_", age.group)] <- rowMeans(counts.mat[cur.gene.ind, cur.cell.ind])  # Take only filtered cells. why? 
        }
        if(stat == "overdispersion")
        {
          od.str = switch(age.group, 
                          "all" = "DispOverall", 
                          "old" = "Disp1", 
                          "young" = "Disp2", 
                          "fc" = "DispFC")
          DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_", age.group)] <- df.od[cur.gene.ind, od.str]
        }
      }
    } # end loop on age groups    
    
    # Add filtering information (global, not specific for age group or expression statistic. May need to change!! )
    DF.expr.stats[[cell.type]][, "filter"] <- cur.gene.ind  # Take only filtered cells
    # Add special expression features of overdispersion 
    if(stat %in% "overdispersion")
      DF.expr.stats[[cell.type]][df.od$GeneName, paste0(stat, "_diff")] <- df.od[, "ResultDiffDisp"]  # Take only filtered cells
  } # end loop on cell types
  
  save(DF.expr.stats, cell_types_categories, cells_ind, file=expression.statistics.outfile) # Save results !!! 
  return(list(DF.expr.stats=DF.expr.stats, cell_types_categories=cell_types_categories, cells_ind=cells_ind))
}


# Computing density for plotting
get_density <- function(x, y, ...) { # function for figures 4 and 5
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# Divide to old and young in each dataset
dataset_to_age_groups <- function(data.type) { # set specific ages for all age groups in all datasets
  old_ages_2 = c()
  if(data.type == "TM.droplet")
  {
    young_ages = "3m" # The ages of the young mice 
    old_ages_1 = c("21m","24m") # The ages of the old mice
    old_ages_2 = c("18m","24m") # secondary old mice age in case no mice in previous "old_ages"
    test_str = "drop 3-24 same-mean"
    young_str = "3m"
    old_str = "21-24m"
  }  
  if(data.type == "TM.facs") 
  {
    young_ages = "3m" # The ages of the young mice 
    old_ages_1 = c("18m", "21m", "24m")
    test_str = "same-mean"
    young_str = "young"
    old_str = "old"
  }
  if(data.type == "CR.Rat")   
    
  {
    young_ages = "Y" # The ages of the young mice 
    old_ages_1 = "O"
    test_str = "rats"
    young_str = "young"
    old_str = "old"
  }
  return(list(young_ages = young_ages, old_ages_1 = old_ages_1, old_ages_2 = old_ages_2, 
              test_str = test_str, young_str = young_str, old_str = old_str))
}


# Get indices of samples 
tissue_to_age_inds <- function(data.type, organ, age.groups, meta.data) { # set specific ages for all age groups in all datasets
  young.ind = c(meta.data$age %in% age.groups$young_ages) # index for cells that came from 3 month old mouses
  old.ind = c(meta.data$age %in% age.groups$old_ages_1) # index for cells that came from old mouses
  if(sum(old.ind) == 0){ # Empty
    old.ind = c(meta.data$age %in% age.groups$old_ages_2)
  }
  if((data.type == "TM.droplet") & (organ=="Lung")){ # special tissue for droplet (which?). Should move this line
    old.ind = meta.data$age %in% c("18m","21m")
  }
  n_cell <- length(meta.data$age)
  all.ind = rep(TRUE, n_cell)
  return(list(young.ind=young.ind, old.ind=old.ind, all.ind=all.ind))
}


# Getting file names for reading BASiCS results 
# cell.type - string representing cell 
dataset_to_BASiCS_file_names <- function(data.type, tissue, cell.type)
{
#  print("Start get BASiCS names")
  groups = dataset_to_age_groups(data.type)
  # Read cell type categories 
  RData_finish = ifelse(data.type == "TM.droplet", "drop 3-24 same-mean.RData","same-mean.RData")
  if(tissue == tolower(tissue)) # if lower, set all upper !
    tissue.cell.type <- paste(toupper(tissue), cell.type)
  else
    tissue.cell.type <- paste(tissue, cell.type)
  if(data.type %in% c("TM.droplet", "TM.facs")) # two technologies in the same dataset
    test_file = paste0(basics.dir, paste(paste0("DVT/",data.type,"/DVT"), tissue.cell.type, RData_finish))   # Differential over-dispersion test results file
  else
    test_file = paste0(basics.dir, paste("DVT/DVT", tissue.cell.type, RData_finish))   # Differential over-dispersion test results file
  
  old_file = paste0(basics.dir, "chains/chain_", paste(tissue.cell.type, groups$old_str), ".Rds")     # Old Markov chain file name
  young_file = paste0(basics.dir, "chains/chain_", paste(tissue.cell.type, groups$young_str), ".Rds")     # Young Markov chain file name
#  print("End get BASiCS names")
  
  return(list(test=test_file, old=old_file, young=young_file))
}


get_DVT_file_name <- function(data.type,tissue,cell_type)
{
  # BASiCS directory and file names
  if(data.type == "CR.Rat"){
    basics.chains.dir =  paste0(main.data.dir, 'Data/CR.Rat/chains/')
    basics.DVT.dir =  paste0(main.data.dir, 'Data/CR.Rat/DVT/')
    DVT.file.name = paste0(basics.DVT.dir,paste("DVT",tissue,cell_type,"same-mean.RData"))
    
  }
  if(data.type == "TM.droplet"){
    basics.chains.dir =  paste0(main.data.dir, 'Data/TabulaMuris/chains/TM.droplet/')
    basics.DVT.dir =  paste0(main.data.dir, 'Data/TabulaMuris/DVT/TM.droplet/')
    DVT.file.name = paste0(basics.DVT.dir,paste("DVT test",tissue,cell_type,"3-24m drop.RData"))
  }
  if(data.type == "TM.facs"){
    basics.chains.dir =  paste0(main.data.dir, 'Data/TabulaMuris/chains/TM.facs/')
    basics.DVT.dir =  paste0(main.data.dir, 'Data/TabulaMuris/DVT/TM.facs/')
    DVT.file.name = paste0(basics.DVT.dir,paste("DVT",tissue,cell_type,"same-mean.RData"))
  }
  
  return(DVT.file.name)
}


tdif <- function(r13,r23,r12,nsize)
{
  
  tdif = (r13 - r23)*sqrt((nsize-3)*(1+r12)/(2*(1-r13^2-r23^2-r12^2 + 2*r13*r23*r12)))
  pdif = 2*(1-pt(abs(tdif),nsize-3)) 
  
  return(c(tdif,pdif))
}

Interaction_reg <- fuunction(Disp_old,Disp_young,Mean_old,Mean_young,gene_names,gene_feature,
                             Organ,Cell_type,feature_name)
{
  
  Disp_both = c(Disp_old,Disp_young) # Over-dispersion for old and young together
  Mean_both = c(Mean_old,Mean_young)# Mean expression for old and young together
  feature_both = rep(gene_feature[gene_names],2) # doubling the current gene feature
  Age = rep(c("Old","Young"),each = length(Disp_old)) # Age vector
  
  # linear regression with interaction
  lm_interaction = lm(log2(Disp_both) ~ feature_both + Age + log2(Mean_both) + Age*feature_both)
  
  s = summary(lm_interaction) # linear regression summary
  
  # sd of the feature coefficient for youngs
  sd_sum = sqrt(vcov(lm_interaction)[2,2]+vcov(lm_interaction)[5,5] + 2*vcov(lm_interaction)[2,5])
  
  beta_sum = sum(s$coefficients[c(2,5),1]) # young feature coefficient
  t_sum = beta_sum/sd_sum # young feature T-statistic
  p_inter = pt(t_sum,df = s$df[2]) # young feature coefficient p_value
  
  
  disp_reg_data_inter = data.frame("Organs" = Organ,"Cell_type" = Cell_type,
                                   lm_interaction$coefficients[2],lm_interaction$coefficients[5],beta_sum,
                                   summary(lm_interaction)$coef[2,4],summary(lm_interaction)$coef[5,4],p_inter)
                                                             
  names(disp_reg_data_inter)[3:8] = c(paste(feature_name,"beta",c("Old","Inter","Young"),sep = "_"),
                                 paste(feature_name,"p_val",c("Old","Inter","Young"),sep = "_"))
  
  return(disp_reg_data_inter)
}

