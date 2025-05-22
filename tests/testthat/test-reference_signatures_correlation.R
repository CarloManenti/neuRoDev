S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))

test_that("Dimensions of output", {
  expect_equal(dim(reference_signatures_correlation(S, refS))[1], dim(S)[2])
})

test_that("Dimensions of output", {
  expect_equal(dim(reference_signatures_correlation(S, refS))[2], dim(refS)[2]+3)
})
