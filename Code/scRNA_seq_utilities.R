# Utilities needed for single cell analysis and visualization 

get_density <- function(x, y, ...) { # function for figures 4 and 5
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
