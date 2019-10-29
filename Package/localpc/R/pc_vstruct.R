#' Obtain v-structures
#'
#' After obtaining the skeleton from the first step of the PC algorithm
#' we use one of Meek's rules to orient the v-structures
#'
#' @param G the current graph after constructing the skeleton
#' @param S the list of separating sets for each node

pc_vstruct <- function(G,S,verbose=TRUE){
  p <- ncol(G)
  for (i in 1:p){
    if (all(G[i,]==0)) next # Node i has no children and should not be considered
    #if (verbose) cat("i: ",i,"\n")
    # Maintain a list of nodes adjacent to i
    i_adj <- which(G[i,] != 0)
    # We want to find all non-adjacent j first
    j_vals <- which(G[i,] == 0)
    j_vals <- j_vals[j_vals!=i]
    # Loop through each non-adjacent j
    for (j in j_vals) {
      # Node j has no children, j is parent to i, or we are repeating an analysis
      # and this j should not be considered
      if (all(G[j,]==0) | G[j,i]==1 | j<i ) next
      #if (verbose) cat("j: ",j,"\n")
      # List of nodes adjacent to j
      j_adj <- which(G[j,] != 0)
      # Final all k that are adjacent to i and j (common neighbors)
      k_vals <- intersect(j_adj,i_adj)
      # If there are no common neighbors, move to next j
      if (length(k_vals)==0){
        next
      }
      # Loop through all common neighbors
      for (k in k_vals){
        if (verbose) cat("k: ",k,"\n")
        # We create a v structure if k is not in the separating set for i and j
        if (!(k %in% S[[i]][[j]])){
          if (verbose) {
            cat("Separation Set:",paste(S[[i]][[j]],collapse = ","),"\n")
            cat("V-Structure: ",i,"->",k,"<-",j,"\n")
          }
          G[k,i] <- 0
          G[k,j] <- 0
        }
      }
    }
  }
  return(G)
}
