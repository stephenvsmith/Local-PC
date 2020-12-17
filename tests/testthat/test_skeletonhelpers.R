context("Testing skeleton helper functions")

test_that("Correctly obtains the edges adjacent to the current node.",{
  data("asiaDAG")
  nodes <- colnames(asiaDAG)
  p <- ncol(asiaDAG)
  expect_equal(nodes[get_current_edges(0,p,asiaDAG)+1],"tub")
  expect_equal(nodes[get_current_edges(which(nodes=="either")-1,p,asiaDAG)+1],c("xray","dysp"))
})

test_that("Combination generation function (combn) works correctly",{
  l <- 0
  vec <- 1:8
  expect_equal(is.na(combn_cpp(vec,l)[1,1]),TRUE)
  
  for (l in 1:5){
    expect_equal(all(combn_cpp(vec,l)==combn(vec,l)),TRUE)
  }
})

test_that("Obtains the correct neighborhood from the reference graph",{
  data("asiaDAG")
  nodes <- colnames(asiaDAG)
  p <- ncol(asiaDAG)
  
  expect_equal(nodes[get_neighbors_from_dag(0,p,asiaDAG)+1],"tub")
  expect_equal(sort(nodes[get_neighbors_from_dag(which(nodes=="either")-1,p,asiaDAG)+1]),c("bronc","dysp","lung","tub","xray"))
})

test_that("Separation test is accurate",{
  data("asiaDAG")
  data("asiadf")
  true_dag <- asiaDAG
  nodes <- colnames(asiaDAG)
  p <- ncol(asiaDAG)
  target <- which(nodes=="either")-1
  neighborhood <- sort(c(target,get_neighbors_from_dag(target,p,asiaDAG)))
  C <- matrix(1,ncol = p,nrow = p)
  diag(C) <- rep(0,p)
  S <- create_conditioning_sets_efficient_cpp2(neighborhood)
  R <- cor(asiadf)
  sep <- 0
  pval <- 0.5
  signif_level <- 0.05
  verbose <- FALSE
  num_tests <- 0
  n <- nrow(asiadf)
  
  i <- 0
  j <- 1
  l <- 0
  kvals <- combn(get_potential_sep(i,j,neighborhood,length(neighborhood),asiaDAG),l)
  
  # Tub should be considered independent of lung when not conditioned on anything else
  check_separation_sample_efficient(l,i,j,
                                    kvals,sep,true_dag,nodes,neighborhood,C,
                                    S,pval,num_tests,R,n,
                                    signif_level,verbose)
  expect_equal(C[1,2],0)
  expect_equal(C[2,1],0)
  
  C[1,2] <- 1
  C[2,1] <- 1
  l <- 1
  kvals <- combn_cpp(which(neighborhood==which(nodes=="either")),l)
  # Tub should not be considered independent of lung when conditioned on either
  # This is due to a v-structure in the graph
  check_separation_sample_efficient(l,i,j,
                                    kvals,sep,true_dag,nodes,neighborhood,C,
                                    S,pval,num_tests,R,n,
                                    signif_level,verbose)
  expect_equal(C[1,2],1)
  expect_equal(C[2,1],1)
  
})