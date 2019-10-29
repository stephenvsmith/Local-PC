#### working through the k_vals

determine_independence <- function(var_list){

  if (is.null(var_list$kvals)){

    ### No separating sets to consider for independence test (independence test)
    var_list <- indep_test_nosep(var_list)

  } else{

    ### Separating sets to consider (conditional independence test)
    var_list <- indep_test_sep(var_list)

  }

  return(var_list)
}
