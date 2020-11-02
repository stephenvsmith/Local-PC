context("Skeleton Specification Construction")

test_that("Separating Set data structure is a list with proper specifications",{
  S <- create_conditioning_sets_efficient_cpp2(c(23,19,10,25,18))
  expect_is(S,"list")
  expect_equal(length(S),5)
  
  expect_is(S[["23"]],"list")
  expect_equal(length(S[["23"]]),4)
  
  expect_equal(names(S),c("23","19","10","25","18"))
  expect_equal(names(S[["10"]]),c("23","19","25","18"))
})

test_that("Initial setup is correct",{
  data("asiaDAG")
  res <- pc_sample_skeleton_setup_efficient_cpp(asiaDAG,5,colnames(asiaDAG),3,FALSE)
  expect_equal(all(res$names[res$neighborhood+1] %in% c("either","bronc","tub","xray","lung","dysp")),TRUE)
  expect_equal(length(res$S),length(res$neighborhood))
  expect_equal(length(res$S[["5"]]),length(res$neighborhood)-1)
})