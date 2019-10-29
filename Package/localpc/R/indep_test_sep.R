### Independence test results with separating sets

indep_test_sep <- function(var_list){

  # Reset the stop variable
  var_list$stop <- FALSE

  # If l is greater than 0
  # loop through the various k sets
  for (k in 1:ncol(var_list$kvals)){

    var_list$k_set <- var_list$kvals[,k]
    var_list <- get_pval(var_list)
    var_list$p.vals <- c(var_list$p.vals,var_list$pval)

    var_list <- test_results(var_list)

    # We may stop after finding a separating set
    if (var_list$stop) break

  }

  return(var_list)

}

