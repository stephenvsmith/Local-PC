# Neighborhood Estimation Function

library(huge)
library(glmnet)

estimate_neighborhood <- function(data,target,method="lasso",rule="AND",rule_1se=FALSE){
  
  if (method=="lasso"){
    est_dag <- matrix(rep(0,p^2),ncol = p,nrow = p)
    neighbors1 <- target_regression(data,target,est_dag,rule_1se)
    order1nghbrs <- neighbors1$neighbors
    est_dag <- neighbors1$est_dag
    
    # if (rule=="AND"){
    #   neighbors <- intersect(neighbors1,neighbors2)
    # } else {
    #   neighbors <- union(neighbors1,neighbors2)
    # }
    
  } else if (method %in% c("mb","ct","glasso","tiger")) {
    if (is.data.frame(data)){
      data <- as.matrix(data)
    }
    estimate <- huge::huge(data,method = method,verbose = FALSE)
    best_est <- huge::huge.select(estimate,criterion = "ric",verbose = FALSE)
    sparse_graph <- best_est$refit
    relationships <- sparse_graph[target,]
    neighbors <- which(relationships != 0)
  }
  if (!exists("neighbors")){
    browser()
  }
  return(neighbors)
}

target_regression <- function(data,target,dag,rule_1se){
  
  nodes <- setdiff(1:ncol(data),target)
  y <- data[,target]
  x <- data[,nodes]
  
  p <- length(nodes)+length(target)

  num_folds <- ifelse(nrow(data)<=30,3,10)
  cat("We will use ",num_folds,"-fold cross-validation to determine the neighbors of node ",target,".\n",sep = "")
  model1 <- glmnet::cv.glmnet(x,y,family = "gaussian",nfolds = num_folds,alpha = 1)
  if (rule_1se){
    coefficients <- coef(model1, s = "lambda.1se")
  } else {
    coefficients <- as.vector(coef(model1))
  }
  coefficients <- coefficients[-1] # Remove the intercept
  nonzero_coef <- which(coefficients!=0)
  
  neighbors <- nodes[nonzero_coef]
  
  results <- update_graph(est_dag,target,nonzero_coef)
  return(list(
    "neighbors"=neighbors,
    "est_dag"=results))
}

update_graph <- function(graph,target,neighbors){
  sapply(neighbors,function(nbr){
    graph[target,nbr] <- 1
    graph[nbr,target] <- 1
  })
  browser()
  return(graph)
}

remaining_regressions <- function(data,target,nodes,rule_1se){
  nodes <- setdiff(nodes,target)
  indices <- sapply(nodes,function(n){
    n2 <- c(setdiff(nodes,n),target)
    neighbors <- target_regression(data,n,n2,rule_1se)
    return(target %in% neighbors)
  })
  potential_neighbors <- nodes[indices]
  return(potential_neighbors)
}

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

