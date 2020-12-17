context("Test Lasso and BIC Neighborhood Selection")

test_that("Lasso and BIC is working",{
  data("asiadf")
  data("asiaDAG")
  data <- as.matrix(asiadf)
  colnames(data) <- colnames(asiaDAG)
  
  for (t in 1:8){
    est_local_graph <- estimate_neighborhood_bic(data,t)
    expect_equal(class(est_local_graph),c("matrix","array"),verbose=FALSE)
  }
})
