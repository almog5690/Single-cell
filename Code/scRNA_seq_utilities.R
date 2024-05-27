# Utilities needed for single cell analysis and visualization 
library(readxl)
library(GGally)
library(dplyr)


# Set directories for raw and processed data, gene features data, analysis results, figures .. for a given data type 
set_data_dirs <- function(data.type)
{
  dirs.index <- which(data.types == data.type)
#  dirs.index <- case_when(data.type == "TM.facs" ~ 1,data.type == "TM.droplet" ~ 2,
#                         data.type == "CR.Rat" ~ 3,data.type == "Age.Anno" ~ 4, 
#                         data.type == "Blood_SC" ~ 5, data.type == "MCA" ~ 6)
  # Set globally: 
  data.dir <<- paste0(main.data.dir, 'Data/', data.dirs[dirs.index], '/')
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
  set_data_dirs(data.type)
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
    file.names = paste(organs, "rds", sep = ".")
  }
  if(data.type == "Blood_SC"){
    organs = "Blood"
    file.names = "Blood.SC.rds"
  }
  if(data.type == "MCA"){
    file.patt <- "MCA"
#    load(paste0(processed.data.dir, "/meta.data.MCA.RData"))
    file.names = list.files(path = processed.data.dir, pattern = "*MCA.rds") # facs processed files list
    organs = unname(sapply(file.names,function(f) unlist(strsplit(f,"[.]"))[1]))
    
    # List of organs with OLD and YOUNG: 
    organs_old_young = c("Bladder", "Brain", "Heart", "Intestine", "Kidney", "Liver", "Lung", "Pancreas", 
          "Prostate", "Spleen", "Stomach", "Testis", "Uterus")
  }
  return(list(file.names = file.names, organs = organs))
}


# Filter cells from an expression matrix (not used yet. Should be part of analysis)
filter_cells <- function(cell_types, young.ind, old.ind, filter.params)
{
  unique_cell_types = unique(cell_types)
  n_cell_types = length(unique_cell_types)  #  max(cell_types) # Number of cell types
  erase = vector()
  
  for(j in 1:n_cell_types){ # filtering cell types with less then 100 cells or less the 20 cells in each the age groups
      ct_ind = unique_cell_types[j]
      if(sum(cell_types==ct_ind, na.rm=TRUE) < filter.params$min.cells.total |
         sum(cell_types==ct_ind & young.ind, na.rm=TRUE) < filter.params$min.cells.per.age |
         sum(cell_types==ct_ind & old.ind, na.rm=TRUE) < filter.params$min.cells.per.age){
        erase = c(erase, j)
        next()
      }
  }
  cells_ind = unique_cell_types[-erase] #  cells_ind = c(1:(n_cell_types+1))[-erase]
  if(length(cells_ind) == 0) cells_ind = unique_cell_types # (1:(n_cell_types+1))
  return(cells_ind) # indices of cells that we keep 
}




# Read feature files (selection, length, GC-content etc. )
# Output will be a list with names being the genes, and values being numerical values.
# Gene names are according to Ensembl (???)
# Input: 
# feature.names - list of features to read 
# organism - allow multiple ones in the future (currently only mice supported)
# force.rerun - default: false, otherwise always re-run to read and parse files 
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
    if(feature.name == "alpha.missense")  # new conservation metric! compare it to gnomad!!
    {
      alpha <- read.delim(paste0(gene.data.dir, "alpha_missense_science.adg7492_data_s4.txt"))
      gene.names = toupper(alpha$gene) # making
      gene.values = alpha$mean_am_pathogenicity
      names(gene.values) = gene.names
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
      # extract age of each gene 
    }
    
    # Features from the impc database
    if(feature.name %in% c("achilles", "impc", "haploinsufficiency"))
    {
        impc.data <- read.csv(paste0(gene.data.dir, "impc_essential_genes_full_dataset.csv"), header=TRUE, sep = ',')
        if(feature.name == "achilles") # human cell line essentiality
          gene.values <- impc.data$achilles_mean_gene_effect
        if(feature.name == "impc") # mouse essentiality
        {
          gene.values.str <- impc.data$impc_via_category
          # Now parse: vital vs. non-vital
          viability.str <-  c("Homozygous-Viable", "Homozygous-Lethal", "Homozygous-Subviable", 
                              "Homozygous-Subviable,Homozygous-Lethal", "Hemizygous-Viable",                      
                              "Homozygous-Subviable,Homozygous-Viable", "Hemizygous-Lethal")
          viability.val <- c(1, 0, 0.5, 0.25, 0.75, 0.75, 0.25)
          mapping_table <- data.frame(strings = viability.str, values = viability.val)
          values_table <- setNames(mapping_table[,2], mapping_table[,1])
          gene.values <- sapply(gene.values.str, function(x) ifelse(x %in% names(values_table), values_table[x], NA))
        }
        if(organism == "mice")
          names(gene.values) <- toupper(impc.data$mouse_symbol)  # take mouse gene names (not human)   
        else  # human
          names(gene.values) <- toupper(impc.data$human_symbol)  # take human gene names (not mouse)   
    }
    gene.values = aggregate(x = gene.values, by = list(names(gene.values)), FUN = mean)  # perform unique
    gv[[feature.name]] <- gene.values$x
    names(gv[[feature.name]]) <- gene.values$Group.1
  }  # end loop on different features   
  save(gv, file=gene.features.outfile)   # Save dataframe to file in R format for faster reading next time  
  
  # new: create also a unified data-frame with all gene names (many missing values: NA)
  all.names <- c() # names(gv[[1]])
  for(feature.name in feature.names)
  {
    all.names <- union(all.names, names(gv[[feature.name]]))
  }
  gv.df <- data.frame(matrix(NA, nrow = length(all.names), ncol = length(feature.names)))
  names(gv.df) <- feature.names
  rownames(gv.df) <- all.names
  for(feature.name in feature.names)
  {
    gv.df[names(gv[[feature.name]]), feature.name] <- gv[[feature.name]]
  }
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
# cell.types - which ones. Default: all cell-types in a given tissue
# expression.stats - which expression statistics to extract 
# age.groups - compute for each group separately
# Need both organ and Seurat output (it doesn't contain the organ/tissue)
# Take list of data frames for all cell types  
# fc.method = How to compute fold-change: Default: Compute ourselves just take log of old/young ratio
#                                         Seurat: (find.markers) - only for mean (currently used for analysis)           
#                                         BASiCS: Can work for both mean and overdispersion (not implemented)
# Output: 
# DF.expr.stats - data frame with mean, overdispersion ... 
# cell_types_categories - 
# cells_ind - ??
extract_expression_statistics <- function(data.type, organ, cell.types=c(), expression.stats = c("mean", "overdispersion"), 
                                          age.groups = c("young", "old", "all", "fc"), 
                                          SeuratOutput=c(), BASiCSOutput = c(), fc.method = "log_old_minus_young", force.rerun = FALSE)
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
    print("Loading file!")
    print(expression.statistics.outfile)
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
    # SeuratOutput = readRDS(file = "D:/Human-Blood/Blood.SC.rds") # HARD-CODED Current tissue Seurat object
  
  print("Finished reading Seurat object")
  print(names(SeuratOutput@meta.data))
  
  list2env(tissue_to_age_inds(data.type, organ, groups, SeuratOutput@meta.data, SeuratOutput), env=environment()) # set specific ages for all age groups in all datasets
  counts.mat = as.matrix(SeuratOutput@assays$RNA@data) # the data matrix for the current tissue
  print("Converted Seurat object to matrix")
  n.genes <- dim(counts.mat)[1]  # same number of genes for all cell types (could have many NAs)
  n.cells <- dim(counts.mat)[2]
  n.groups <- length(age.groups)
  
  if(data.type == "CR.Rat"){ # Convert to dummy variables 
    cell_types = as.numeric(SeuratOutput@meta.data$cell_types) # Cell types vector
    cell_types_categories = levels(SeuratOutput$cell_types)
  } else if (data.type == "Blood_SC")
  {
    cell_types = SeuratOutput@meta.data$CT # Cell types vector
    cell_types_categories = meta.data[[organ.ind]]$cell_ontology_class # Cell type names. Missing variable meta.data.drop
  } else if (data.type == "MCA") # Other data types, including MCA (?)
  {
    cell_types = SeuratOutput@meta.data$CT # Cell types vector
    cell_types_categories = unique(SeuratOutput@meta.data$CT) # Cell type names. Missing variable meta.data.drop
  } else
  {
    cell_types = SeuratOutput@meta.data$cell.ontology.class # Cell types vector
    cell_types_categories = meta.data[[organ.ind]]$cell_ontology_class # Cell type names. Missing variable meta.data.drop
  }
  cells_ind = filter_cells(cell_types, young.ind, old.ind, filter.params)  #  unique(cell_types)+1 # NO FILTERING NOW !!!
  n.cell.types <- length(cells_ind) # number of cell types in tissue

      
  for(cell.type in 1:length(cells_ind))  # First load data if not loaded already 
  {
    print(paste0("Parse data for cell type: ", cell_types_categories[cell.type]))
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
      
      # Extract overdispersion features
      df.od = test@Results$Disp@Table
      df.od$GeneName = toupper(df.od$GeneName) # Uppercasing gene names
      df.od = df.od[!df.od$ResultDiffDisp %in% c("ExcludedFromTesting"),] # Filtering genes significant difference in mean expression
    }  # end if overdispersion (read from BASiCS)
    
    print("Go over age groups")
    for(age.group in age.groups)  # loop on age groups
    {
      cur.cell.ind = switch(age.group, # Indices of cells in each group
                       "all" = all.ind, 
                       "young" = young.ind, 
                       "old" = old.ind, 
                       "fc" = young.ind)  & (cell_types==cells_ind[cell.type]) # take only cell type . Assume that cell typea are numbers starts with zero 
      cur.cell.ind[is.na(cur.cell.ind)] = FALSE
      # Filter genes with low expression for old/young (less than 10 counts). Keep indices of genes FILTER BY ALL!!!
      if(age.group == "fc")  # union young or old 
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
          if(age.group == "fc")  # union young or old. fc = log(old) - log(young) 
          {
            # TODO: ENABLE READING FROM BASICS OUTPUT INSTEAD OF LOG-DIFFERENEC
            if(fc.method == "log_old_minus_young")
              DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_fc")] <- log( # add regularization
                (DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_old")] + filter.params$pseudo.count) / 
                  (DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_young")] + filter.params$pseudo.count)) # Take only filtered cells. why? 
            if(fc.method == "seurat") # works only for mean .. 
            {
              DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_fc")] <- 9999 # XXX ?? ADD CODE 
            }
            if(fc.method == "basics") # can work for both mean and od, currently used only for od. 
            {
              load(BASiCS.files$test) # load BASICS object 
              df.mean = test@Results$Mean@Table
              df.mean$GeneName = toupper(df.mean$GeneName) # Uppercasing gene names
              DF.expr.stats[[cell.type]][cur.gene.ind , paste0(stat, "_fc")] = df.mean$MeanFC # ?? ADD
            }
              
          } else             # Repeat for all age groups and features 
          {
            if(sum(cur.cell.ind) == 1)
              DF.expr.stats[[cell.type]][cur.gene.ind , paste0(stat, "_", age.group)] <- 
                counts.mat[cur.gene.ind, cur.cell.ind]
            else
              DF.expr.stats[[cell.type]][cur.gene.ind , paste0(stat, "_", age.group)] <- 
                rowMeans(counts.mat[cur.gene.ind, cur.cell.ind])  # Take only filtered cells. why? 
          }
        }
        
        if(stat == "overdispersion")
        {
          od.str = switch(age.group, 
                          "all" = "DispOverall", 
                          "old" = "Disp1", 
                          "young" = "Disp2", 
                          "fc" = "DispFC")
          # Find common gene indices, and intersect. Copy in the right order!!! 
          common.od.gene.names = intersect(toupper(names(cur.gene.ind)), toupper(df.od$GeneName))
          rownames(df.od) <- df.od$GeneName
          DF.expr.stats[[cell.type]][common.od.gene.names  , paste0(stat, "_", age.group)] <- df.od[common.od.gene.names, od.str]  # Copy overdispersion information 
          
#          DF.expr.stats[[cell.type]][cur.gene.ind  , paste0(stat, "_", age.group)] <- df.od[cur.gene.ind, od.str]  # Copy overdispersion information 
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
  print("Finished reading and saving data matrix")
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
  young_str = "young"
  old_str = "old"
  
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
  }
  if(data.type == "CR.Rat")   
    
  {
    young_ages = "Y" # The ages of the young mice 
    old_ages_1 = "O"
    test_str = "rats"
  }
  if(data.type == "Blood_SC"){
    young_ages = "Young" # The ages of the young mice 
    old_ages_1 = "Old"
    test_str = "BloodHuman"
  }
  
  if(data.type == "MCA"){  # Need to fill! 
    young_ages = "Young"
    old_ages_1 = "Old"
    test_str = "MCA"
  }
  
  return(list(young_ages = young_ages, old_ages_1 = old_ages_1, old_ages_2 = old_ages_2, 
              test_str = test_str, young_str = young_str, old_str = old_str))
}


# Get indices of samples (old, young ..)
# Input: 
# data.type - which dataset
# organ - tissue name
# age.groups - set groups for young and old
# meta.data - structure for meta data information
tissue_to_age_inds <- function(data.type, organ, age.groups, meta.data, SC=c()) { # set specific ages for all age groups in all datasets
  if(data.type == "Age.Anno")  # ignore groups and meta data. Remove very young. Set middle group as 'young'
  {
    young.ind = grepl("mid",Idents(SC))
    old.ind = grepl("old",Idents(SC))
    n_cell <- length(Idents(SC))
  } else if(data.type == "Blood_SC")
  {
    young.ind = c(meta.data$Age %in% age.groups$young_ages) # index for cells that came from 3 month old mouses
    old.ind = c(meta.data$Age %in% age.groups$old_ages_1) # index for cells that came from old mouses
    if(sum(old.ind) == 0){ # Empty
      old.ind = c(meta.data$Age %in% age.groups$old_ages_2)
    }
    n_cell <- length(meta.data$age)
  } else if(data.type == "MCA")
  {
    old.ind = SC$Age == "Old"
    young.ind = SC$Age == "Young"
    n_cell <- length(SC$CT)
  }  else  # All other datasets 
  {
    young.ind = c(meta.data$age %in% age.groups$young_ages) # index for cells that came from 3 month old mouses
    old.ind = c(meta.data$age %in% age.groups$old_ages_1) # index for cells that came from old mouses
    if(sum(old.ind) == 0){ # Empty
      old.ind = c(meta.data$age %in% age.groups$old_ages_2)
    }
    if((data.type == "TM.droplet") & (organ=="Lung")){ # special tissue for droplet (which?). Should move this line
      old.ind = meta.data$age %in% c("18m","21m")
    }
    n_cell <- length(meta.data$age)
  } # if human
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


# Determine file name for BASiCS data files 
get_DVT_file_name <- function(data.type, tissue, cell_type)
{
  dirs.index <- which(data.types == data.type)
#  dirs.index <- case_when(data.type == "TM.facs" ~ 1,data.type == "TM.droplet" ~ 2,
#                          data.type == "CR.Rat" ~ 3,data.type == "Age.Anno" ~ 4, 
#                          data.type == "Blood_SC" ~ 5, data.type == "MCA" ~ 6)
  
  
  if(data.type %in% c("TM.facs", "TM.droplet"))  # special case 
  {
    basics.dir <- paste0(main.data.dir, 'Data/TabulaMuris/BASiCS/')  # This is already in the config !! 
    basics.chains.dir =  paste0(basics.dir, 'chains/', data.type, "/")
    basics.DVT.dir =  paste0(basics.dir, 'DVT/', data.type, "/")
    if(data.type == "TM.droplet")
      DVT.file.name = paste0(basics.DVT.dir, paste("DVT test",tissue,cell_type,"3-24m drop.RData"))
    else
      DVT.file.name = paste0(basics.DVT.dir, paste("DVT",tissue,cell_type,"same-mean.RData"))
  } else  # all other datasets 
  {
    basics.dir <- paste0(main.data.dir, 'Data/', data.dirs[dirs.index], '/BASiCS/')
    basics.chains.dir =  paste0(basics.dir, 'chains/')
    basics.DVT.dir =  paste0(basics.dir, 'DVT/')
    DVT.file.name = paste0(basics.DVT.dir, paste("DVT",tissue,cell_type,"same-mean.RData"))
  }  
  
  return(list(DVT.file.name=DVT.file.name, basics.chains.dir = basics.chains.dir, basics.DVT.dir = basics.DVT.dir))
}


#######################################################################
# Function for testing the difference between two Spearman correlations 
# Input: 
# r13, r23, r12 - three pairwise correlations
# nsize - sample size
# spearman_dif_t = function(r13, r23, r12, nsize){
tdif <- function(r13,r23,r12,nsize)
{
  tdif = (r13 - r23)*sqrt((nsize-3)*(1+r12)/(2*(1-r13^2-r23^2-r12^2 + 2*r13*r23*r12)))
  pdif = 2*(1-pt(abs(tdif),nsize-3)) 
  return(c(tdif,pdif))
}


# Function for testing the difference between two regression coefficients (TBD!)
beta_diff_t = function(beta_1, X, Y)
{
  beta.diff <- 0
  return(True)  
}
  
 
# Interactions  
Interaction_reg <- function(Disp_old,Disp_young,Mean_old,Mean_young,gene_names,gene_feature,
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


# Post-processing utilities for tables ...
post_process_reg_table <- function(paper.DF, output.df.file.name = "", use.beta = TRUE)
{
  stat_str = ifelse(use.beta, "beta", "cor")
    
    n.data.types = length(data.types)
  n.rows = 2*n.paper.reg # young,old and fc for each reg. type  # + 2 for regressing on age (need to add separately)
  
  model.df = data.frame(matrix(ncol = 3, nrow = n.rows))
  colnames(model.df) <- c("Model", "Covariate", "Hypothesis")
  for(i in 1:n.paper.reg)  # run different regression types 
  {
    model.df[i*2-1, "Model"] <- paste0(paper.expression.stat.y[[i]],  " ~ ", paper.features[[i]])
    model.df[i*2, "Model"] <- paste0("d(", paper.expression.stat.y[[i]],  ") ~ ", paper.features[[i]])
    
    model.df[i*2-1, "Covariate"] <- paper.expression.stat.x[[i]]
    if(paper.expression.stat.x[[i]] != "")
      model.df[i*2, "Covariate"] <- paste0("d(", paper.expression.stat.x[[i]], ")")
    
    model.df[i*2-1, "Hypothesis"] <- paper.hypothesis[[i]]
    model.df[i*2, "Hypothesis"] <- paper.fc.hypothesis[[i]]
  }
    
  column_names <- paste(rep(data.types, each=3),  rep(c("pos/neg", "y>o/o>y", "pval"), 3), sep=":")
  output.df <- data.frame(matrix(ncol = length(column_names), nrow = n.rows))
  colnames(output.df) <- column_names
  
  data.ctr = 1
  for(data.type in data.types)  # loop on datasets
  {
    print(paste0("Post-process Regression Results data: ", data.type))
    n.cell.types = nrow(paper.DF[[data.ctr]][[1]])
    for(i in 1:n.paper.reg)  # run different regression types 
    {
      num.young.bigger <- sum(abs(paper.DF[[data.ctr]][[i]] [paste0(paper.features[[i]], "_young_", stat_str)]) > 
            abs(paper.DF[[data.ctr]][[i]] [paste0(paper.features[[i]], "_old_", stat_str)]))
      num.old.bigger <- n.cell.types - num.young.bigger
      num.pos <- sum(paper.DF[[data.ctr]][[i]] [paste0(paper.features[[i]], "_young_", stat_str)] > 0) + 
        sum(paper.DF[[data.ctr]][[i]] [paste0(paper.features[[i]], "_old_", stat_str)] > 0)
      num.neg <- sum(paper.DF[[data.ctr]][[i]] [paste0(paper.features[[i]], "_young_", stat_str)] <= 0) + 
        sum(paper.DF[[data.ctr]][[i]] [paste0(paper.features[[i]], "_old_", stat_str)] <= 0)
      
      output.df[i*2-1, data.ctr*3-2] = ifelse(num.pos > num.neg, 
                                              paste0("pos: ", as.character(num.pos), '/', as.character(num.neg)), 
                                              paste0("neg: ", as.character(num.pos), '/', as.character(num.neg)))
      
      output.df[i*2-1, data.ctr*3-1] = ifelse(num.young.bigger > num.old.bigger,
                                              paste0("young: ", num.young.bigger, "/", num.old.bigger), 
                                              paste0("old: ", num.young.bigger, "/", num.old.bigger)) # F young/old counts
      output.df[i*2-1, data.ctr*3] = 2*min(pbinom(num.young.bigger, n.cell.types, 0.5), 
                                           pbinom(num.old.bigger, n.cell.types, 0.5))  # G young/old pvalue. Binomial two sided test
      
      num.fc.pos <- sum(paper.DF[[data.ctr]][[i]][paste0(paper.features[[i]], "_fc_", stat_str)] > 0)
      num.fc.neg <- n.cell.types - num.fc.pos
      
      tmp_str = paste0(num.fc.pos, "/", num.fc.neg)
      output.df[i*2, data.ctr*3-2] = ifelse(num.fc.pos > num.fc.neg, paste0("pos: ", tmp_str), paste0("neg: ", tmp_str)) # fold change
      output.df[i*2, data.ctr*3-1] = ""
      #      output.df[i*2, data.ctr*3-1] = paste0(num.fc.pos, "/", num.fc.neg)  # fold-change pos/neg
      output.df[i*2, data.ctr*3] = t.test(paper.DF[[data.ctr]][[i]][paste0(paper.features[[i]], "_fc_", stat_str)])$p.val # t-test for all fc coefficients: is the mean significantly different from zero? 
    }
    data.ctr <- data.ctr + 1
  }
  
  if (!file.exists(res.dir)) 
    dir.create(res.dir)
  
  output.df <- cbind(model.df, output.df) # add few columns with model names 
  write.csv(output.df, output.df.file.name) # Return and save to file 
  return(output.df)
  
}

# Mean filtering function
gene.filtering = function(gene_mean_old, gene_mean_young, gene_mean_all = NULL,
                          expression.thresh = 0.2, FC_filter = F, SC, CT_name, 
                          Old.indicator, Young.indicator){
  if(FC_filter){ # FC-filtering using FindMarkers function
    Idents(SC) <- "Cells"
    CT_SC = subset(SC,idents = CT_name) # current cell type Seurat object
    
    # if(all(is.na( Cells(CT_SC)[CT_SC$age %in% Old.indicator]))) next()
    # Clustering using FindMarkers
    age_clusters = FindMarkers(CT_SC,ident.1 = Cells(CT_SC)[CT_SC$age %in% Old.indicator],ident.2 = Cells(CT_SC)[CT_SC$age %in% Young.indicator],
                               test.use = "wilcox",assay = "RNA",slot = "data",pseudocount.use = 0.1,verbose = T,min.pct = 0.5,logfc.threshold = log2(1.25))
    age_clusters$bh_p_val = p.adjust(age_clusters$p_val,method = "BH") # Adjusted p_values
    gene_mean_FC = age_clusters$avg_log2FC # FC estimations
    names(gene_mean_FC) = toupper(row.names(age_clusters)) # gene names
    
    return(gene_mean_FC[age_clusters$bh_p_val <= 0.1])
    
  } else { # Mean filtering (young/old or all)
    ages_mean = (gene_mean_old + gene_mean_young)/2 # genes averages average  
    
    if(is.null(gene_mean_all)){ # young/old filtering
      return(c("Young" = gene_mean_young[which(ages_mean > expression.thresh)], 
               "Old" = gene_mean_old[which(ages_mean > expression.thresh)]))
      
    } else { # all filteering
      return(gene_mean_all[which(ages_mean > expression.thresh)])
    }
    
  }
}
