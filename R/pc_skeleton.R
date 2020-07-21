#' PC Algorithm for Skeleton (for Local PC)
#'
#' @param dataset default is NULL, only used when we are using sample version and have data
#' @param C_tilde default is NULL, used in the local algorithm with each clique being considered forming a complete graph
#' @param cor.mat default is NULL, the correlation matrix for the nodes (population algorithm)
#' @param lmax the maximum size of the separating sets that we will consider (default is 3)
#' @param verbose messages concerning the operation of the algorithm (default is TRUE)
#' @param tol the significance level (default is 0.005)
#'
#' @export

pc_skeleton <- function(dataset=NULL,C_tilde=NULL,true_dag=NULL,
                        pop=TRUE,lmax=3,fci_step1,
                        verbose=TRUE,verbose_small=TRUE,tol=0.01){

  # Initial setup
  var_list <- pc_skel_setup(dataset,true_dag,C_tilde,
                            pop,lmax,
                            verbose,verbose_small,tol)

  # Initialize size of subset
  var_list[["l"]] = -1
  # Initialize adjacency matrix for finished graph
  var_list[["C"]] <- C_tilde
  
  var_list$fci_step1 <- fci_step1

  ### loop over all l values (until lmax)
  while (var_list[["l"]] < var_list[["lmax"]]){

    var_list[["l"]] <- var_list[["l"]] + 1

    if (verbose) cat("l value is: ",var_list[["l"]],"\n")

    ### loop over all i values
    for (i in 1:var_list$p){
      var_list <- skeleton_target(i,var_list)
    }

  }
  return(list("adjacency"=var_list$C,"sep_sets"=var_list$S,"p_vals"=var_list$p.vals,"num_tests"=var_list$num_tests))
}
