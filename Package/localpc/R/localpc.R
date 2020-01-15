#' Local PC Algorithm (Efficient)
#'
#' @param data default is NULL, the sample version when we have data
#' @param cor.mat default is NULL, the correlation matrix for the nodes (population version)
#' @param target vector of nodes of interest
#' @param true_dag given for the population version
#'
#' @export

localpc <- function(data=NULL,true_dag=NULL,target,G=NULL,lmax=3,tol=0.01,pop=TRUE,verbose = TRUE,verbose_small=TRUE){

  build <- build_initial_graph(target,true_dag)
  
  Ctilde <- G

  neighbors <- build$neighbors

  skel_res <- pc_skel_loc(dataset = data,C_tilde = Ctilde,
                          true_dag=true_dag,neighbors=neighbors,
                          pop = pop,lmax = lmax,verbose = verbose,tol = tol)

  G <- pc_vstruct(G = skel_res$adjacency,S = skel_res$sep_sets,verbose=verbose)

  S <- skel_res$sep_sets
  p_vals_vec <- skel_res$p_vals

  return(list("G"=G,"S"=S,"p_vals"=p_vals_vec))
}
