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
  basics.dir <<- paste0(analysis.results.dir, 'BASiCS/')  # For BASiCS output 
  
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
read_gene_features  <- function(feature.names)
{
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
  return(gv)
}



# Compute features of expression (mean, overdispersion, variance ... )
# for a given tissue/cell type
# cell.types - which ones. Default: all types in a given tissue
# expression.stats - which expression statistics to extract 
# Need both organ and Seurat output (it doesn't contain the organ/tissue)
# Take list of data frames for all cell types  
extract_expression_statistics <- function(data.type, expression.stats = c("mean", "overdispersion"), 
                                          organ, SeuratOutput=c(), cell.types=c(), force.rerun = FALSE)
{
  
  set_data_dirs(data.type)
  expression.statistics.outfile <- paste0(analysis.results.dir, 'mean.analysis.', data.type, '.', 
                                  paste0( feature.types, collapse="_"), '.RData')
  samples <- get_tissue_file_names(data.type)
  meta.data = get_meta_data(data.type)
  if(file.exists(mean.analysis.outfile) & (force.rerun==FALSE))
  {
    load(mean.analysis.outfile)
    return(DF_cor)
  }
  
  
  # First load data if not loaded already 
  if(length(SeuratOutput)==0) # empty
    SeuratOutput = readRDS(file = paste0(processed.data.dir, organ, ".", processed.files.str[data.type], ".rds")) # Current tissue Seurat object
  counts.mat = as.matrix(SC@assays$RNA@data) # the data matrix for the current tissue
  
  # For overdispersion read DVT files if they exist/run basics  
  
  
  # Next, extract mean
  for(stat in expression.stats)
  {
    if(stat == "mean")
    {
      # Repear for all age groups and features 
      gene.mean.by.age.group[age.group] <- rowMeans(counts.mat[cur_gene_name[[feature.type]], cur.ind])  # Take only filtered cells
      
    }
    if(stat == "overdispersion")
    {
      
    }
  }
  
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
dataset_to_BASiCS_file_names <- function(data.type, tissue, cell.type)
{
  groups = dataset_to_age_groups(data.type)
  # Read cell type categories 
  RData_finish = ifelse(data.type == "TM.droplet","drop 3-24 same-mean.RData","same-mean.RData")
  tissue.cell.type <- paste(tissue,cell_type)
  test_file = paste0(basics.dir, paste(paste0("DVT/",data.type,"/DVT"), tissue.cell.type, RData_finish))   # Differential over-dispersion test results file
  old_file = paste0(basics.dir, "chains/chain_", paste(tissue.cell.type, groups$old_str), ".Rds")     # Old Markov chain file name
  young_file = paste0(basics.dir, "chains/chain_", paste(tissue.cell.type, groups$young_str), ".Rds")     # Young Markov chain file name

  return(list(test=test_file, old=old_file, young=young_file))
}
