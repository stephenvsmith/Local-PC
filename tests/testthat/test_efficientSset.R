context("Testing efficient separating set functions")

test_that("Get correct neighbors and check separating set update",{
  data("asiaDAG")
  nodes <- colnames(asiaDAG)
  p <- ncol(asiaDAG)
  target <- 5 # "either" is the target
  neighbors <- get_neighbors_from_dag(target,p,asiaDAG)
  neighbor_names <- nodes[neighbors+1]
  expect_equal(nodes[target+1],"either")
  expect_equal(setdiff(neighbor_names,c("tub","lung","xray","dysp","bronc")),character(0))
  
  set.seed(111)
  neighborhood_full <- sample(1:8,8,replace = FALSE)
  S <- create_conditioning_sets_efficient_cpp2(neighborhood_full)
  change_S_0_efficient(S,neighborhood_full[6],neighborhood_full[2])
  change_S_0_efficient(S,neighborhood_full[2],neighborhood_full[6])
  expect_equal(neighborhood_full[c(6,2)],c(1,3))
  expect_equal(S[["1"]][["3"]],-1)
  expect_equal(S[["3"]][["1"]],-1)
  
  neighborhood <- c(target,neighbors)
  expect_equal(length(neighborhood),6)
  expect_equal(neighborhood,c(5, 4, 1, 6, 3, 7))
  S <- create_conditioning_sets_efficient_cpp2(neighborhood)
  i <- 2
  j <- 6
  kvals <- 4:5
  change_S_efficient(S,neighborhood[i],neighborhood[j],kvals)
  change_S_efficient(S,neighborhood[j],neighborhood[i],kvals)
  expect_equal(S[[as.name(neighborhood[i])]][[as.name(neighborhood[j])]],kvals)
  expect_equal(S[[as.name(neighborhood[j])]][[as.name(neighborhood[i])]],kvals)

})
