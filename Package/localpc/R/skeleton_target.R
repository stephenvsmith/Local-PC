### This is the other loop for the PC Algorithm step for finding the skeleton

skeleton_target <- function(i,var_list){

  var_list[["i"]] <- i
  if (var_list$verbose) cat("i value is: ",i,"\n")

  # Determine adjacent values for possible j
  var_list[["jvals"]] <- which(var_list[["C"]][i,] != 0)

  ### loop over all adjacent j values
  for (j in var_list$jvals){
    var_list <- check_neighbor(j,var_list)
  }

  return(var_list)
}

#### Check potential neighbor for target during PC Algorithm skeleton step

check_neighbor <- function(j,var_list){

  var_list[["j"]] <- j
  if (var_list$verbose) cat("j value is: ",j,"\n")

  # (i,j) is the pair we are considering
  # adj is the remaining adjacent nodes to i not equal to j
  adj <- find_neighbors(i,var_list$true_dag)
  var_list$adj <- adj[adj!=j]

  if (length(var_list$adj)>=var_list$l){ # Makes sure there is a valid amount of adjacent nodes for each step

    var_list <- get_k_vals(var_list)
    var_list <- determine_independence(var_list)

  }

  return(var_list)
}


#' Find Neighbors
#'
#' This function finds all the neighbors for a particular node.
#' A neighbor is defined as the set of parents, children, and spouses of the node.
#'
#' @param t the node of interest
#' @param true_dag provides the adjacency matrix for the true DAG
#' Note for the DAG adjacency matrix: A_{ij}=1 if i -> j


find_neighbors <- function(t,true_dag){
  # Parents of node t
  parents <- which(true_dag[,t]==1)
  # Children of node t
  children <- which(true_dag[t,]==1)
  # Spouses share the same children
  spouses <- c()
  for (c in children){ # loop through all of t's children
    # Find all parents of current child
    potential <- which(true_dag[,c]==1)
    # Add all parents that are not t to spouses
    spouses <- c(spouses,potential[potential!=t])
  }
  return(unique(c(parents,children,spouses)))
}

#### Get potential separating sets

get_k_vals <- function(var_list){

  if (var_list$l==0){
    var_list$kvals <- NULL
  } else if (length(var_list$adj)==1){
    var_list$kvals <- matrix(var_list$adj,nrow = 1,ncol = 1)
  } else {

    # If there is more than one possible adjacent node
    # build our various ksets for conditional independence
    var_list$kvals <- combn(var_list$adj,var_list$l)

  }


  return(var_list)

}
