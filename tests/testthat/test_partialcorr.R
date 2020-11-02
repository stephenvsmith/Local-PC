context("Testing partial correlation functions")

test_that("Partial Correlation function is accurate",{
  data("asiadf")
  C <- cor(asiadf)
  true_result <- cor(asiadf)[1,3]
  est_result <- get_partial_correlation(cor(asiadf),0,2,double())
  expect_equal(est_result,true_result)
  
  
  true_result <- ppcor::pcor(asiadf)
  for (i in 0:7){
    for (j in 0:7){
      est_result <- get_partial_correlation(C,i,j,setdiff(0:7,c(i,j)))
      expect_equal(est_result,true_result$estimate[i+1,j+1])
    }
  }
  
})

test_that("Make sure you are obtaining the correct test results",{
  C <- cor(asiadf)
  k1 <- numeric(0)
  n <- nrow(asiadf)
  alpha <- 0.05
  for (i in 0:7){
    for (j in 0:7){
      true_res <- pcalg::gaussCItest(i+1,j+1,k1,list("C"=C,"n"=n))
      true_res2 <- pcalg::condIndFisherZ(i+1,j+1,k1, C, n, qnorm(1-alpha/2), verbose= FALSE)
      est_res <- condIndTest(C,i,j,k1,n,signif_level=alpha)
      #cat("i =",i,"| j =",j,"| k =","| pval (Est):",est_res$pval,"| pval (True): ",true_res,"| Estimate:",est_res$result,"| True Result:",true_res2,"\n")
      expect_equal(est_res$pval,true_res)
      expect_true(est_res$result == true_res2)
      
      set.seed(111)
      k2 <- sample(setdiff(0:7,c(i,j)),2,replace = FALSE)
      true_res <- pcalg::gaussCItest(i+1,j+1,k2+1,list("C"=C,"n"=n))
      true_res2 <- pcalg::condIndFisherZ(i+1,j+1,k2+1, C, n, qnorm(1-alpha/2), verbose= TRUE)
      est_res <- condIndTest(C,i,j,k2,n,signif_level=alpha)
      if (est_res$result != true_res2){
        cat("\ni =",i,"| j =",j,"| k =",paste(k2,collapse = " "),"| pval (Est):",est_res$pval,"| pval (True): ",true_res,"| Estimate:",est_res$result,"| True Result:",true_res2,"\n")
      }
      expect_equal(est_res$pval,true_res)
      expect_true(est_res$result == true_res2)
      
      k3 <- sample(setdiff(0:7,c(i,j)),3,replace = FALSE)
      true_res <- pcalg::gaussCItest(i+1,j+1,k3+1,list("C"=C,"n"=n))
      true_res2 <- pcalg::condIndFisherZ(i+1,j+1,k3+1, C, n, qnorm(1-alpha/2), verbose= TRUE)
      est_res <- condIndTest(C,i,j,k3,n,signif_level=alpha)
      if (est_res$result != true_res2){
        cat("\ni =",i,"| j =",j,"| k =",paste(k3,collapse = " "),"| pval (Est):",
            est_res$pval,"| pval (True): ",true_res,
            "| Estimate:",est_res$result,"| True Result:",true_res2,"\n")
      }
      expect_equal(est_res$pval,true_res)
      expect_true(est_res$result == true_res2)
      
      k4 <- setdiff(0:7,c(i,j))
      true_res <- pcalg::gaussCItest(i+1,j+1,k4+1,list("C"=C,"n"=n))
      true_res2 <- pcalg::condIndFisherZ(i+1,j+1,k4+1, C, n, qnorm(1-alpha/2), verbose= TRUE)
      est_res <- condIndTest(C,i,j,k4,n,signif_level=alpha)
      if (est_res$result != true_res2){
        cat("\ni =",i,"| j =",j,"| k =",paste(k4,collapse = " "),"| pval (Est):",est_res$pval,"| pval (True): ",true_res,"| Estimate:",est_res$result,"| True Result:",true_res2,"\n")
      }
      expect_equal(est_res$pval,true_res)
      expect_true(est_res$result == true_res2)
    }
  }

})