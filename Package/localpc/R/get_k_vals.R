#### Get potential separating sets

get_k_vals <- function(var_list){

  if (var_list$l==0){
    var_list$kvals <- NULL
  } else if (length(var_list$adj)==1){
    var_list$kvals <- matrix(var_list$adj,nrow = 1,ncol = 1)
  } else {

    # If there is more than one possible adjacent node
    # build our various ksets for conditional independence
    var_list$kvals <- combn(var_list$adj,var_list$l)

  }


  return(var_list)

}
