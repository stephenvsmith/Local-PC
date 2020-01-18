#' Local PC Algorithm (Efficient)
#'
#' @param data default is NULL, the sample version when we have data
#' @param cor.mat default is NULL, the correlation matrix for the nodes (population version)
#' @param target vector of nodes of interest
#' @param true_dag given for the population version
#'
#' @export

localpc <- function(data=NULL,true_dag=NULL,target,G=NULL,lmax=3,tol=0.01,pop=TRUE,verbose = TRUE,verbose_small=TRUE){

  skel_res <- pc_skeleton(dataset = data,C_tilde = G,true_dag=true_dag,
                          pop = pop,lmax = lmax,verbose = verbose,tol = tol)

  G_new <- pc_vstruct(G = skel_res$adjacency,S = skel_res$sep_sets,verbose=verbose)

  S <- skel_res$sep_sets
  p_vals_vec <- skel_res$p_vals

  return(list("G"=G_new,"S"=S,"p_vals"=p_vals_vec))
}
