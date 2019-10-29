############################################################################################################

# Setup for the PC Algorithm (Skeleton)

############################################################################################################

pc_skel_setup <- function(dataset,true_dag,C_tilde,
                          neighbors,pop,lmax,
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
              "true_dag"=true_dag,"neighbors"=neighbors,
              "pop"=pop,"lmax"=lmax,
              "verbose"=verbose,"tol"=tol,
              "n"=n))
}
