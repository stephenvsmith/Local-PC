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

### results from independence test (no separating set)

indep_test_nosep <- function(var_list){
  # This is when l = 0
  # check correlation between node i and node j
  var_list <- get_pval(var_list)
  var_list$p.vals <- c(var_list$p.vals,var_list$pval)

  var_list <- test_results(var_list)

  return(var_list)
}


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

#' Get p-value for Node Correlation
#'
#' @param cor.mat the correlation matrix for all the nodes
#' @param i the potential parent node
#' @param j the potential child node
#' @param k_set the potential separating set


get_pval <- function(var_list){

  if (var_list$pop){

    ## Set up in order to use d-separation
    names <- colnames(var_list$true_dag)
    e <- bnlearn::empty.graph(names)
    true_dag <- as.matrix(var_list$true_dag)
    bnlearn::amat(e) <- true_dag

    # Determine independence
    if (var_list$l == 0){ # if there is no set being conditioned on

      var_list$pval <- as.numeric(bnlearn::dsep(e,names[var_list$i],names[var_list$j]))

    }else{
      #if ( var_list$i == 41 | var_list$j == 41) browser()
      var_list$pval <- as.numeric(bnlearn::dsep(e,names[var_list$i],names[var_list$j],names[var_list$k_set]))
    }

  } else {

    # We first deal with the case where there is no separating set
    if(var_list$l == 0){

      #cat("i: ",i,",j: ",j,", No separating set\n")
      var_list$pval <- pcalg::gaussCItest(var_list$i,var_list$j,NULL,
                                          list("C"=var_list$cor.mat,"n"=var_list$n))

    } else {
      # case where there is a separating set

      #cat("i: ",i,",j: ",j,", k: ",paste(k,collapse = ","),"\n")
      #browser()
      # x <- var_list$dataset[,var_list$i]
      # y <- var_list$dataset[,var_list$j]
      # #if (length(var_list$k)) browser()
      # z <- as.matrix(var_list$dataset[,var_list$k_set])
      # lm1 <- lm(x ~ z)
      # lm2 <- lm(y ~ z)
      # rho <- cor(lm1$residuals,lm2$residuals)
      # Z_ij <- 0.5*log((1+rho)/(1-rho))

      var_list$pval <- pcalg::gaussCItest(var_list$i,var_list$j,var_list$k_set,
                                          list("C"=var_list$cor.mat,"n"=var_list$n))
    }
  }

  var_list$num_tests <- var_list$num_tests + 1

  return(var_list)
}

