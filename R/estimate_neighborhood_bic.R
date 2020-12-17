# Neighborhood Estimation Function

library(gglasso)
library(glmnet)

###############################################################
# Function returns an estimated adjacency matrix for the
# graph (undirected) using the Lasso for various values of
# lambda and choosing the best one using BIC
###############################################################
estimate_neighborhood_bic <- function(data,target,verbose=TRUE){

  p <- ncol(data)
  est_dag <- matrix(rep(0,p^2),ncol = p,nrow = p)
  neighbors1 <- target_regression_bic(data,target,est_dag,verbose)
  order1nbrs <- neighbors1$neighbors
  est_dag <- neighbors1$est_dag

  for (nbr in order1nbrs){
    est_dag <- target_regression_bic(data,nbr,est_dag,verbose)$est_dag
  }

  return(est_dag)
}

###############################################################
# Workhorse function that runs lasso and updates the estimated
# adjacency matrix graph so that it contains the estimated
# neighbors of the current target.
###############################################################
target_regression_bic <- function(data,target,dag,verbose){

  nodes <- setdiff(1:ncol(data),target)
  y <- data[,target]
  x <- data[,nodes]

  p <- length(nodes)+length(target)

  model1 <- glmnet(x,y,family = "gaussian")
  nonzero_coef <- get_neighbors_bic(model1,nrow(data))

  neighbors <- nodes[nonzero_coef]
  if (verbose){
    if (!is.null(colnames(data))){
      cat("Neighbors of",colnames(data)[target],"are",paste(colnames(data)[neighbors],collapse = ", "),"\n")
    }
  }

  results <- update_graph(dag,target,neighbors)

  return(list(
    "neighbors"=neighbors,
    "est_dag"=results))
}

###############################################################
# Helper function to update the estimated adjacency matrix
# based on the target and its neighbors
###############################################################
update_graph <- function(graph,target,neighbors){

  for (nbr in neighbors){
    graph[target,nbr] <- 1
    graph[nbr,target] <- 1
  }

  return(graph)
}

###############################################################
# Connects each node in the neighborhood with each other
# (Undirected edges)
###############################################################
neighborhood_graph <- function(p,neighbors){
  g <- matrix(0,nrow = p,ncol = p)
  for (i in neighbors){
    n2 <- setdiff(neighbors,i)
    for (j in n2){
      g[i,j] <- 1
    }
  }
  return(g)
}

###############################################################
# Returns which coefficients are not equal to 0
# Coefficients are taken from model with lowest BIC
###############################################################
get_neighbors_bic <- function(model,n){
  index <- which.min(deviance(model)+model$df*log(n))
  coeff <- model$beta[,index]
  neighbors <- which(coeff != 0)
  return(neighbors)
}

