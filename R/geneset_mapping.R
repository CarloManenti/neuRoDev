#' GeneSet mapping
#'
#' @param mapping_obj The result of the mapping
#' @param geneset_matrix The matrix containing the average expression of
#' a set of genesets over clusters
#' @param show_top_n The number of genesets to show in the Heatmap with
#' highest signature values
#' @param avoid_repetitions Boolean variable
#' If some genesets have the same name but different
#' cell types, only the geneset with the first cell type will be considered
#' @param cluster_columns A boolean, if TRUE columns are clustered.
#' @param column_order The order of the columns. If NULL, no specific order is
#' given, the default is kept.
#'
#' @return A List with the geneset scores, the geneset scores signatures
#' and the Heatmap
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' annotated_M <- reference_signatures_correlation(S, refS)
#' new_clusterS <- FC_signatures(matrix(runif(80,0,10), ncol = 4))
#' rownames(new_clusterS) <- paste0('Gene-', seq(1, dim(new_clusterS)[1]))
#' colnames(new_clusterS) <- paste0('New-', seq(1, dim(new_clusterS)[2]))
#' new_M <- reference_signatures_correlation(new_clusterS, refS)
#' mapping_obj <- add_to_reference(annotated_M,
#' new_M,
#' annotated_M$`Best.Assignment`)
#' anno_matrix <- matrix(runif(100,0,10), ncol = 10)
#' colnames(anno_matrix) <- colnames(S)
#' rownames(anno_matrix) <- paste0('Annotation-', seq(1, dim(anno_matrix)[1]))
#' anno_matrix <- t(t(anno_matrix)/colSums(anno_matrix))
#' geneset_mapping(mapping_obj, anno_matrix, avoid_repetitions = FALSE)
geneset_mapping <- function(mapping_obj,
                            geneset_matrix,
                            show_top_n = 3,
                            avoid_repetitions = TRUE,
                            cluster_columns = TRUE,
                            column_order = NULL) {

  if (!is.null(mapping_obj$New$umap_obj$refined_network)) {
    index_matrix <- mapping_obj$New$umap_obj$umap_out$refined_network$indexes
    distance_matrix <- mapping_obj$New$umap_obj$umap_out$refined_network$distances

    alt_index <- mapping_obj$Original$umap_obj$umap_out$refined_network$indexes
  } else {
    index_matrix <- mapping_obj$New$umap_obj$umap_out$knn$indexes
    distance_matrix <- mapping_obj$New$umap_obj$umap_out$knn$distances

    alt_index <- mapping_obj$Original$umap_obj$umap_out$knn$indexes

  }

  new_clusters <- rownames(index_matrix)[which(!rownames(index_matrix) %in% rownames(alt_index))]

  new_cluster_distances <- distance_matrix[new_clusters,]

  new_cluster_similarities <- do.call(cbind, lapply(new_clusters, function(i) {

    dist_i <- distance_matrix[i,]
    names(dist_i) <- rownames(index_matrix)[index_matrix[i,]]

    dist_i <- dist_i[which(!names(dist_i) %in% new_clusters)]

    sim_i <- 1-dist_i

    sim_i_ordered <- sim_i[colnames(geneset_matrix)]

  }))

  colnames(new_cluster_similarities) <- new_clusters

  new_cluster_similarities <- t(t(new_cluster_similarities)/colSums(new_cluster_similarities))

  scores <- geneset_matrix %*% new_cluster_similarities

  if(avoid_repetitions) {
    rownames(scores) <- unlist(lapply(strsplit(rownames(scores),
                                               '.',
                                               fixed = TRUE),
                                      function(i) {i[2]}))
    idxs <- unlist(lapply(unique(rownames(scores)),
                          function(i) {which(rownames(scores) == i)[1]}))
    scores <- scores[idxs,]
  }

  scores_sigs <- FC_signatures(scores)

  highest_idxs <- apply(scores_sigs, 2, function(i) {
    order(i, decreasing = TRUE)[seq(1, show_top_n)]
  })

  all_idxs <- unique(as.vector(highest_idxs))

  h1 <- ComplexHeatmap::Heatmap(scores_sigs[all_idxs,],
                                width = grid::unit(min(5*ncol(scores_sigs), 200), 'mm'),
                                height = grid::unit(min(5*nrow(scores_sigs[all_idxs,]), 200), 'mm'),
                                rect_gp = grid::gpar(col = 'black'), name = 'GeneSet\nScore\nSignature',
                                cluster_columns = cluster_columns,
                                column_order = column_order)

  h1 <- ComplexHeatmap::draw(h1,
             heatmap_legend_side = 'left',
             annotation_legend_side = 'left')

  return(S4Vectors::List('Score' = scores,
                         'ScoresSignatures' = scores_sigs,
                         'Heatmap' = h1))

}
