S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
M <- reference_signatures_correlation(S, refS)
umap <- umap_graph_clustering(M)
l <- umap$umap_out$layout
l <- cbind(l, 'Size' = rep(1, dim(l)[1]))
l <- cbind(l, 'Color' = rep(c('blue', 'green'), each=5))
l <- cbind(l, 'Group' = rep(c('A','B'), each=5))

test_that("Type output", {
  expect_s3_class(umap_pointsize(l, color_attr = rep(c('A','B'), each=5)), 'ggplot')
})
