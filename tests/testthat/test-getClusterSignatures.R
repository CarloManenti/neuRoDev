M <- matrix(runif(20000,0,100), ncol=200)
rownames(M) <- paste0('Gene-', seq(1, dim(M)[1]))
colnames(M) <- paste0('Cell-', seq(1, dim(M)[2]))

test_that("Returns a list", {
  expect_s4_class(getClusterSignatures(M, resolution=2),
                  'SimpleList')
})

test_that("Length of 9 of returned list", {
  expect_equal(length(getClusterSignatures(M, resolution=2)),
               9)
})
