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
      Z_ij <- 0.5*log((1+var_list$cor.mat[var_list$i,var_list$j])/(1-var_list$cor.mat[var_list$i,var_list$j]))

    } else {
      # case where there is a separating set

      #cat("i: ",i,",j: ",j,", k: ",paste(k,collapse = ","),"\n")
      #browser()
      x <- var_list$dataset[,var_list$i]
      y <- var_list$dataset[,var_list$j]
      #if (length(var_list$k)) browser()
      z <- as.matrix(var_list$dataset[,var_list$k_set])
      lm1 <- lm(x ~ z)
      lm2 <- lm(y ~ z)
      rho <- cor(lm1$residuals,lm2$residuals)
      Z_ij <- 0.5*log((1+rho)/(1-rho))
    }
    var_list$pval <- 2*(1-pnorm(sqrt(var_list$n-length(var_list$k_set)-3)*abs(Z_ij)))
  }

  return(var_list)
}
