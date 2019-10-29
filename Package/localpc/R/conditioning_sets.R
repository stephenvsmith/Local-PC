#' Conditioning Sets Function
#'
#' This function creates a list of lists for each node, with each node's list containing
#' lists that contain vectors indicating the separation sets for the two nodes.
#'
#' @param p is the number of nodes in the graph

create_conditioning_sets <- function(p){

  # Create the list that we will return with the conditioning sets
  S <- list()
  for (i in 1:p){ # loop through all nodes
    # Create a list for each node
    S[[i]] <- list()
    for (j in 1:p){
      # List i,j contains the separation sets for nodes i and j
      # Default NA - i,j are assumed to be independent
      S[[i]][[j]] <- NA
    }
  }

  return(S)
}
