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


# Read feature files
# Output will be a list with names being the genes, and values being numerical values.
# Gene names are according to Ensembl (???)
read_gene_features  <- function(feature.name)
{
  if(feature.name == "selection")
  {
    select <- read.delim(paste0(gene.data.dir, "gnomad.v2.1.1.lof_metrics.by_gene.txt"))
    gene.names = toupper(select$gene) # making all the gene names in to upper case letters
    gene.values = select$pLI # the selection score vector
    names(gene.values) = gene.names # naming each score with the gene it belongs to
  }
  if(feature.name == "gene.len")
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
  
  return(gene.values)
}

# Computing density for plotting
get_density <- function(x, y, ...) { # function for figures 4 and 5
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
