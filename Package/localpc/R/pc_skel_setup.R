############################################################################################################

# Setup for the PC Algorithm (Skeleton)

############################################################################################################

pc_skel_setup <- function(dataset,true_dag,C_tilde,
                          pop,lmax,
                          verbose,tol){

  # number of nodes (Use dataset for sample, true_dag for population)
  p <- ifelse(pop,ncol(true_dag),ncol(dataset))

  # number of observations (if using dataset for sample version)
  n <- ifelse(pop,NA,nrow(dataset))

  if (is.null(C_tilde)){ # We are doing complete PC algorithm
    # Adjacency matrix for complete graph
    C_tilde <- matrix(1,nrow = p,ncol = p)
    diag(C_tilde) <- rep(0,p)
  }

  # We will use S to store the separation sets for each pair of variables
  S <- create_conditioning_sets(p)

  # Get sample correlation matrix if this is sample PC Algorithm
  cor.mat <- NULL
  if (!is.null(dataset)){
    cor.mat <- cor(dataset)
  }

  p.vals <- c() # stores the p-values for the tests

  return(list("p"=p,"C_tilde"=C_tilde,
              "S"=S,"cor.mat"=cor.mat,
              "p.vals"=p.vals,"dataset"=dataset,
              "true_dag"=true_dag,
              "pop"=pop,"lmax"=lmax,
              "verbose"=verbose,"tol"=tol,
              "n"=n))
}



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

