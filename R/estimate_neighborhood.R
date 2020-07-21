# Neighborhood Estimation Function

library(huge)
library(glmnet)

estimate_neighborhood <- function(data,target,method="lasso",rule="AND",rule_1se=FALSE){

  if (method=="lasso"){
    nodes <- setdiff(1:ncol(data),target)
    neighbors1 <- target_regression(data,target,nodes,rule_1se)
    neighbors2 <- remaining_regressions(data,target,nodes,rule_1se)
    if (rule=="AND"){
      neighbors <- intersect(neighbors1,neighbors2)
    } else {
      neighbors <- union(neighbors1,neighbors2)
    }
    
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

target_regression <- function(data,target,nodes,rule_1se){

  y <- data[,target]
  x <- data[,nodes]
  if (is.data.frame(x)){
    x <- as.matrix(x)
  }
  num_folds <- ifelse(nrow(data)<=30,3,10)
  model1 <- glmnet::cv.glmnet(x,y,family = "gaussian",nfolds = num_folds,alpha = 1)
  if (rule_1se){
    coefficients <- coef(model1, s = "lambda.1se")
  } else {
    coefficients <- as.vector(coef(model1))
  }
  coefficients <- coefficients[-1] # Remove the intercept
  nonzero_coef <- which(coefficients!=0)
  
  neighbors <- nodes[nonzero_coef]
  return(neighbors)
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

