context("Skeleton Specification Construction")

test_that("Separating Set data structure is a list with proper specifications",{
  S <- create_conditioning_sets_efficient_cpp(c(23,19,10,25,18))
  expect_is(S,"list")
  expect_equal(length(S),5)
  
  expect_is(S[["23"]],"list")
  expect_equal(length(S[["23"]]),4)
  
  expect_equal(names(S),c("23","19","10","25","18"))
  expect_equal(names(S[["10"]]),c("23","19","25","18"))
})