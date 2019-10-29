#### Check potential neighbor for target during PC Algorithm skeleton step

check_neighbor <- function(j,var_list){

  var_list[["j"]] <- j
  if (var_list$verbose) cat("j value is: ",j,"\n")

  # (i,j) is the pair we are considering
  # adj is the remaining adjacent nodes to i not equal to j
  adj <- var_list$neighbors[[as.name(var_list$i)]]
  var_list$adj <- adj[adj!=j]

  if (length(var_list$adj)>=var_list$l){ # Makes sure there is a valid amount of adjacent nodes for each step

    var_list <- get_k_vals(var_list)
    var_list <- determine_independence(var_list)

  }

  return(var_list)
}
