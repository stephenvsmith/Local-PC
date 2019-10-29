#' Building the initial graph
#'
#' This helps us to build the subgraphs for the local pc algorithm
#' We need the true DAG to provide the first and second order neighbors
#' of the set of target nodes
#'
#' @param target a vector of nodes for which we want to conduct local pc on
#' @param true_dag a matrix containing all the information about the true dag
#'

build_initial_graph2 <- function(target,true_dag){
  # Ctilde is our initial graph
  Ctilde <- matrix(0,ncol=ncol(true_dag),nrow=nrow(true_dag))
  if (!is.null(rownames(true_dag))){
    rownames(Ctilde) <- rownames(true_dag)
    colnames(Ctilde) <- colnames(true_dag)
  }
  neighbors <- target # Tracks all the neighbors we have to find neighborhoods for (second order neighbors)

  neighbors_list <- list()

  neighborhood <- list() # tracks which neighborhood each node is located in

  ### We will now find the first order neighbors of each target node
  for (t in target){ # Build a subgraph of t U Neighborhood(t)
    tmp_list <- build_subgraph_one_target(t,Ctilde,true_dag)
    Ctilde <- tmp_list[[1]]; neighbors <- union(neighbors,tmp_list[[2]])
    if (length(target) > 1) neighborhood[[t]] <- tmp_list[[2]]
  }

  ### List containing the neighbors of each of the target's neighbors
  neighbors_list <- sapply(neighbors,find_neighbors,true_dag,USE.NAMES = TRUE,simplify = FALSE)
  #neighbors_list2 <- sapply(unique(unlist(neighbors)), find_neighbors,true_dag,USE.NAMES = TRUE,simplify = FALSE)
  names(neighbors_list) <- neighbors

  ### Remaining connections between neighbors of neighbors from different neighborhoods
  if (length(target)>1){
    Ctilde <- remaining_edges(target,neighborhood,neighbors_list,Ctilde)
    Ctilde <- finalize_edges(target,neighborhood,neighbors_list,Ctilde)
  }

  return(list("Ctilde"=Ctilde,"neighbors"=neighbors_list))
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


####################################################

# build subgraph for one target

####################################################

build_subgraph_one_target <- function(t,Ctilde,true_dag) {

  fo_n <- find_neighbors(t,true_dag) # obtain first-order neighbors of the target
  neighbors <- union(t,fo_n)

  ### We now create a complete graph for this clique

  # Start with connecting the first-order neighbors
  Ctilde <- connect_fo(fo_n,Ctilde)

  # Connecting target node with first-order neighbors
  Ctilde <- connect_fo_target(t,fo_n,Ctilde)

  return(list(Ctilde,neighbors))

}


##############################################

# Connect first-order neighbors

##############################################

connect_fo <- function(fo_n,Ctilde) {

  for (f in fo_n){
    # We connect first-order neighbors to each other
    # But we prevent self-loops
    for (ff in fo_n[fo_n!=f]){
      Ctilde[f,ff] <- 1
      Ctilde[ff,f] <- 1
    }
  }

  return(Ctilde)
}

##############################################

# Connect first-order neighbors to target

##############################################

connect_fo_target <- function(t,fo_n,Ctilde) {
  # Connecting target node with first order neighbors
  for (f in fo_n){
    Ctilde[t,f] <- 1
    Ctilde[f,t] <- 1
  }

  return(Ctilde)
}

##############################################

# Connecting the neighbors of one target variable together

##############################################

connect_neighbors_1target <- function(neighbors,Ctilde){
  for (n in neighbors){
    # Prevent self-loops
    for (n2 in neighbors[neighbors!=n]){
      Ctilde[n,n2] <- 1
      Ctilde[n2,n] <- 1
    }
  }
  return(Ctilde)
}


######################################################################################################

# Connect neighbors that are currently in different neighborhoods

######################################################################################################

remaining_edges <- function(target,neighborhood,neighbors_list,Ctilde) {

  for (t in target){
    #if (t == 5) browser()
    for (n in setdiff(neighborhood[[t]],t)) { # n is a neighbor of target t
      for (tprime in setdiff(target,t)){ # tprime is a different target node
        nprime_possibilities <-
        for (nprime in neighborhood[[tprime]]){ # nprime is a neighbor of tprime or is tprime
          if (nprime %in% c(n,neighbors_list[[as.name(n)]])){ # if nprime is in the neighborhood of a neighbor of t
            for (n2 in neighborhood[[t]]){ # connecting nprime with all the neighbors of t, where n and nprime are neighbors in different target neighborhoods
              if (n2 != nprime){
                Ctilde[n2,nprime] <- 1
                Ctilde[nprime,n2] <- 1
              }
            }
          }
        }
      }
    }
  }


  return(Ctilde)
}

#############################################################################################

# If all target nodes share a neighbor,

#############################################################################################

finalize_edges <- function(target,neighborhood,neighbors_list,Ctilde) {
  n_hoods <- lapply(target,function(t) return(neighborhood[[t]]))
  for (i in 1:length(target)){
    for (j in setdiff(1:length(target),i)){
      for (ni in n_hoods[[i]]){ # neighbors of target i
        for (nj in setdiff(n_hoods[[j]],ni)){ # neighbors of target j
          for (nii in neighbors_list[[as.name(ni)]]){
            for (njj in neighbors_list[[as.name(nj)]]) {
              if (nii == njj){ # common neighbors of neighbors of targets i and j
                Ctilde[ni,nj] <- 1
                Ctilde[nj,ni] <- 1
              }
            }
          }
        }
      }
    }
  }
  return(Ctilde)
}
