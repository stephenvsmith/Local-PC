###################################################################################################

# Author: Stephen Smith
# Date Created: 10/2/2019
# Description: Functions to set up everything we need to build network results

###################################################################################################

ObtainPathStarts <- function(ps1="/Users/stephensmith/",
                             ps2="/home/stephen/",
                             ps3="Dropbox/Academics/Research/") {
  
  # Establish which computer this is
  computer_name <- Sys.info()["nodename"]
  my_laptop <- computer_name == "Stephens-MacBook-Pro.local"
  
  # Set up paths for saving files
  path_start <- ifelse(my_laptop,ps1,ps2)
  ps2 <- ifelse(my_laptop,
                paste0(ps1,ps3),
                paste0(ps2,ps3))
  
  return(c(path_start,ps2))
  
}

LoadLibraries <- function() {
  # Load appropriate libraries
  require(bnlearn)
  if (!my_laptop){
    setwd("/home/stephen/Desktop/Research/package")
    install("localpc")
  }
  require(localpc)
  
  return()
}


GenerateDataGrid <- function(){
  # Grid for generating data
  data.grid = data.frame(network = "asia",
                         data.type = "continuous",
                         n.dat = 1,
                         n.obs = 1000,
                         c.ratio = 0,
                         max.in.degree = Inf,
                         lb = 0.5,  # lower bound of coefficients
                         ub = 1,  # upper bound of coefficients
                         low = 0,  # lower bound of variances if continuous, of number of levels if discrete
                         high = 1,  # upper bound of variances if continuous, of number of levels if discrete
                         scale = TRUE,
                         stringsAsFactors = FALSE)
  return(data.grid)
}

# This function removes all variables that aren't used
SimplifyDAG <- function(mat,target=NULL){
  zeros <- rep(0,nrow(mat))
  rm_ind <- c()
  for (i in 1:nrow(mat)){
    if (all(mat[i,]==zeros) & all(mat[,i]==zeros)){
      if (!is.null(target)){
        if (i %in% target)
          next
        else rm_ind <- c(rm_ind,i)
      }
    }
  }
  return(rm_ind)
}
