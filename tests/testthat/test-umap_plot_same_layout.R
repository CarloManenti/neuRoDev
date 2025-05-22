S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
annotated_M <- reference_signatures_correlation(S, refS)
new_clusterS <- FC_signatures(matrix(runif(80,0,10), ncol = 4))
rownames(new_clusterS) <- paste0('Gene-', seq(1, dim(new_clusterS)[1]))
colnames(new_clusterS) <- paste0('New-', seq(1, dim(new_clusterS)[2]))
new_M <- reference_signatures_correlation(new_clusterS, refS)

test_that("Length output", {
  expect_equal(length(umap_plot_same_layout(annotated_reference = annotated_M, signatures_cor = new_M, color_attr = annotated_M$`Best.Assignment`)), 4)
})

test_that("Names output", {
  expect_identical(names(umap_plot_same_layout(annotated_reference = annotated_M, signatures_cor = new_M, color_attr = annotated_M$`Best.Assignment`)), c('Original', 'New', 'MappingQuality', 'NearestNeighborsAnnotation'))
})

test_that("Type output", {
  expect_s3_class(umap_plot_same_layout(annotated_reference = annotated_M, signatures_cor = new_M, color_attr = annotated_M$`Best.Assignment`)$Original$umap_plot, 'ggplot')
})
