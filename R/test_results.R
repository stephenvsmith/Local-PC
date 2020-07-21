### Finish results if p-value is high, thus meaning we accept null hypothesis that they are conditionally independent
### H0: rho = 0

test_results <- function(var_list){

  if (var_list$pval >= var_list$tol){

    if (is.null(var_list$kvals)){

      ### Fail to reject the null hypothesis (independence test)
      if (var_list$verbose_small){
        cat("--------------------------------------\n")
        cat("The p-value is: ",var_list$pval,"\n")
        cat("Removing edge from ",var_list$i," to ",var_list$j,"\n")
        cat("Separation Set:","empty\n")
        cat("--------------------------------------\n")
      }

      # Delete the edge
      var_list$C[var_list$i,var_list$j] <- 0
      var_list$C[var_list$j,var_list$i] <- 0

      # No separating sets for independence test
      var_list$S[[var_list$i]][[var_list$j]] <- "empty"
      var_list$S[[var_list$j]][[var_list$i]] <- "empty"

      # Remove j from possible adjacent vertices to i
      var_list$jvals <- var_list$jvals[-which(var_list$jvals==var_list$j)]
    } else {

      # Also fail to reject the null hypothesis (conditional independence)

      if (var_list$verbose_small){
        cat("--------------------------------------\n")
        cat("The p-value is: ",var_list$pval,"\n")
        cat("Removing edge from ",var_list$i," to ",var_list$j,"\n")
        cat("Separation Set:",paste(var_list$k_set,collapse = ","),"\n")
        cat("--------------------------------------\n")
      }

      # Delete the edge
      var_list$C[var_list$i,var_list$j] <- 0
      var_list$C[var_list$j,var_list$i] <- 0

      # Add k as the separation set
      var_list$S[[var_list$i]][[var_list$j]] <- var_list$k_set
      var_list$S[[var_list$j]][[var_list$i]] <- var_list$k_set

      # Remove j from the possible adjacent values
      var_list$jvals <- var_list$jvals[-which(var_list$jvals==var_list$j)]
      var_list$stop <- TRUE # We don't need to check any more k sets


    }

  }


  return(var_list)

}
