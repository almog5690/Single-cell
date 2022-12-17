# Utilities needed for single cell analysis and visualization 


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
    file.names = list.files(path = processed.data.dir, pattern = "*.rds") # facs processed files list
    file.names = file.names[!(file.names %in% list.files(path = processed.data.dir, pattern = "drop.rds"))]
    organs  = unname(sapply(file.names,function(f) unlist(strsplit(f,"[.]"))[1])) #facs organs names
#    organs[5] = "Brain_Non-Myeloid" # why needed?
  }
    
  if(data.type == "CR.Rat")  # to fill 
  {
  }
  
  return(list(file.names = files.names, organs = organs))
}


# Computing density for plotting
get_density <- function(x, y, ...) { # function for figures 4 and 5
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
