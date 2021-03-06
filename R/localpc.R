#' Local PC Algorithm (Efficient)
#'
#' @param data default is NULL, the sample version when we have data
#' @param cor.mat default is NULL, the correlation matrix for the nodes (population version)
#' @param target vector of nodes of interest
#' @param true_dag given for the population version
#'
#' @export

localpc <- function(data=NULL,true_dag=NULL,target,G=NULL,lmax=3,tol=0.01,
                    pop=TRUE,verbose = TRUE,verbose_small=TRUE,
                    orient_v=TRUE,fci_step1=FALSE){

  skel_res <- pc_skeleton(dataset = data,C_tilde = G,true_dag=true_dag,
                          pop = pop,lmax = lmax,fci_step1 = fci_step1,
                          verbose = verbose,verbose_small = verbose_small,tol = tol)
  if (orient_v){
    G_new <- pc_vstruct(G = skel_res$adjacency,S = skel_res$sep_sets,verbose=verbose)
  } else {
    G_new <- skel_res$adjacency
  }

  S <- skel_res$sep_sets
  p_vals_vec <- skel_res$p_vals

  return(list("G"=G_new,"S"=S,"p_vals"=p_vals_vec,"num_tests"=skel_res$num_tests))
}


localpc_cpp <- function(data=NULL,true_dag=NULL,target,G=NULL,lmax=3,tol=0.01,
                        verbose = TRUE,verbose_small=TRUE){
  node_names <- colnames(data)
  if (is.data.frame(data)){
    data <- as.matrix(data)
  }
  
  if (is.data.frame(true_dag)){
    true_dag <- as.matrix(true_dag)
  }
  
  # We change the target to target - 1 in order to accommodate change to C++
  return(pc_sample_efficient_cpp(true_dag,data,target-1,node_names,lmax,1-tol,verbose,verbose_small=FALSE))

}

localpc_est_nbhd_cpp <- function(data=NULL,target,G=NULL,lmax=3,tol=0.01,
                                 rule="AND",rule_1se=FALSE,cv_folds=10,
                                 verbose = TRUE){
  node_names <- colnames(data)
  if (is.data.frame(data)){
    data <- as.matrix(data)
  }
  
  # We need to estimate the target's first- and second-order neighbors
  true_dag_est <- estimate_neighborhood_bic(data,target,verbose=verbose)
  
  # We change the target to target - 1 in order to accommodate change to C++
  return(pc_sample_efficient_cpp(true_dag_est,data,target-1,node_names,lmax,tol,verbose,verbose_small=FALSE))
  
}
