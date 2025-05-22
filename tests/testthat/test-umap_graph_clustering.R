S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
M <- reference_signatures_correlation(S, refS)

test_that("Length output", {
  expect_equal(length(umap_graph_clustering(M)), 7)
})

test_that("Returns a graph", {
  expect_s3_class(umap_graph_clustering(M)$umap_knn_igraph, 'igraph')
})
