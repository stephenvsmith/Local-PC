### This is the other loop for the PC Algorithm step for finding the skeleton

skeleton_target <- function(i,var_list){

  var_list[["i"]] <- i
  if (var_list$verbose) cat("i value is: ",i,"\n")

  # Determine adjacent values for possible j
  var_list[["jvals"]] <- which(var_list[["C"]][i,] != 0)

  ### loop over all adjacent j values
  for (j in var_list$jvals){
    var_list <- check_neighbor(j,var_list)
  }

  return(var_list)
}
