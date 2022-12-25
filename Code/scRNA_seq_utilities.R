# Utilities needed for single cell analysis and visualization 
library(readxl)

# Get tissues. Use the processed files if they exist (otherwise need to go back to raw files?)
get_tissue_file_names <- function(data.type)
{
  # loading the meta data names for the droplet technology
  # load("meta.data.drop.Rdata")
  # load("facs.files.RData")
  # load("meta.data.facs.RData")
  
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
  }
  
  return(list(file.names = file.names, organs = organs))
}


# Filter cells from an expression matrix
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


# Read feature files
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
    if(feature.name %in% c("gene.len", "gene.length"))
    {
      # reading the transcript length data
      length_data = read.delim(paste0(gene.data.dir, "mart_export.txt"))
      gene.values = length_data[,7] # gene length vector
      names(gene.values) = toupper(length_data[,6]) # Add gene names. Note: there are many values per gene 
    }
    if(feature.name == "TATA")
    {
      tata.data <- read.delim(paste0(gene.data.dir, "mm10_tata_annotation_v3.3.tsv"))
      # TATA BOX
      gene.values <- as.integer(unlist(lapply(tata.data$TATA.Box.TBP..Promoter.Homer.Distance.From.Peak.sequence.strand.conservation., nchar)) > 0)
      gene.values2 <- tata.data$GC.
      gene.values2 <- tata.data$CpG.
      names(gene.values) <- tata.data$Gene.Name
    }
    if(feature.name == "mRNA.half.life")
    {
      half.life.data <-  read_excel(paste0(gene.data.dir, "13059_2022_2811_MOESM3_ESM.xlsx"), sheet = "mouse", skip =2)
      gene.values <- half.life.data$Half_life_PC1
      names(gene.values) <- half.life.data$Gene.name
    }
    
    gene.values = aggregate(x = gene.values, by = list(names(gene.values)), FUN = mean)  # perform unique
    gv[[feature.name]] <- gene.values$x
    names(gv[[feature.name]]) <- gene.values$Group.1
  }  
  
  return(gv)
}

# Computing density for plotting
get_density <- function(x, y, ...) { # function for figures 4 and 5
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


dataset_to_age_groups <- function(data.type) { # set specific ages for all age groups in all datasets
  old_ages_2 = c()
  if(data.type == "TM.droplet")
  {
    young_ages = "3m" # The ages of the young mice 
    old_ages_1 = c("21m","24m") # The ages of the old mice
    old_ages_2 = c("18m","24m") # secondary old mice age in case no mice in previous "old_ages"
  }  
  if(data.type == "TM.facs") 
  {
    young_ages = "3m" # The ages of the young mice 
    old_ages_1 = c("18m", "21m", "24m")
  }
  if(data.type == "CR.Rat")   
  {
    young_ages = "Y" # The ages of the young mice 
    old_ages_1 = "O"
  }
  
  return(list(young_ages = young_ages, old_ages_1 = old_ages_1, old_ages_2 = old_ages_2))
}



