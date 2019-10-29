### results from independence test (no separating set)

indep_test_nosep <- function(var_list){
  # This is when l = 0
  # check correlation between node i and node j
  var_list <- get_pval(var_list)
  var_list$p.vals <- c(var_list$p.vals,var_list$pval)

  var_list <- test_results(var_list)

  return(var_list)
}
