test_that("Length output", {
  expect_equal(length(Proportional_sampling(matrix(runif(20000,0,10), ncol=200), membership_vector = rep(c('A','B','C','D'), each = 50), Ntotal=20, ReturnSCE = FALSE)), 2)
})

test_that("Dim Matrix", {
  expect_equal(dim(Proportional_sampling(matrix(runif(20000,0,10), ncol=200), membership_vector = rep(c('A','B','C','D'), each = 50), Ntotal=20))[2], 20)
})
