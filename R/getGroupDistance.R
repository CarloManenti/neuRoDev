#' getGroupDistance
#'
#' It computes the mean, max and minimum distances between clusters of each group
#' given a membership label, an annotated reference (like those returned after
#' correlating signatures) and a umap object
#'
#' @param annotated_reference The annotated reference, like those returned by
#' reference_signatures_correlation
#' @param membership_label A membership label that identifies a column in the
#' annotated reference
#' @param umap_obj A umap object like those returned by umap_signature_plot
#'
#' @return Three matrices and three heatmaps, one for the mean distance, one for
#' the max distance and one for the min distance
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' M <- reference_signatures_correlation(S, refS)
#' group = M$`Best.Assignment`
#' umap_obj = umap_signature_plot(M,
#' color_attr = group,
#' label_attr = group)
#' getGroupDistance(M, 'Best.Assignment', umap_obj)
getGroupDistance <- function(annotated_reference, membership_label, umap_obj) {

  membership <- annotated_reference[[membership_label]]

  if('umap_obj' %in% names(umap_obj)) {
    network <- umap_obj$umap_obj$refined_network
  } else {
    network <- umap_obj$refined_network
  }

  whole_distance_matrix <- igraph::distances(network, weights = 1-igraph::E(network)$weight)

  mean_distances_matrix <- matrix(0, ncol = length(unique(membership)), nrow = length(unique(membership)))

  colnames(mean_distances_matrix) <- unique(membership)
  rownames(mean_distances_matrix) <- unique(membership)

  max_distances_matrix <- mean_distances_matrix
  min_distances_matrix <- mean_distances_matrix

  for(i in unique(membership)) {
    cl1 <- annotated_reference$Cluster[which(membership == i)]

    for(x in unique(membership)) {

      if(mean_distances_matrix[x,i] != 0) {
        next()
      }

      cl2 <- annotated_reference$Cluster[which(membership == x)]

      distance <- whole_distance_matrix[cl1, cl2]

      if(length(cl1) == 1 && length(cl2) == 1) {
        mean_distances_matrix[i,x] <- mean(distance)
        mean_distances_matrix[x,i] <- mean(distance)
        max_distances_matrix[i,x] <- max(distance)
        max_distances_matrix[x,i] <- max(distance)
        min_distances_matrix[i,x] <- min(distance)
        min_distances_matrix[x,i] <- min(distance)
        next()
      }

      mean_distances_matrix[i,x] <- mean(distance[which(distance != 0)])
      mean_distances_matrix[x,i] <- mean(distance[which(distance != 0)])
      max_distances_matrix[i,x] <- max(distance[which(distance != 0)])
      max_distances_matrix[x,i] <- max(distance[which(distance != 0)])
      min_distances_matrix[i,x] <- min(distance[which(distance != 0)])
      min_distances_matrix[x,i] <- min(distance[which(distance != 0)])

    }

  }

  constant <- min(120/nrow(mean_distances_matrix), 5)

  h1 <- ComplexHeatmap::Heatmap(1-mean_distances_matrix,
                                name = 'Mean similarity',
                                col = grDevices::blues9,
                                rect_gp = grid::gpar(col = 'black'),
                                width = ncol(mean_distances_matrix)*grid::unit(constant, "mm"),
                                height = nrow(mean_distances_matrix)*grid::unit(constant, "mm"),
                                column_names_gp = grid::gpar(cex = max(min(1,30/nrow(mean_distances_matrix)),0.45)),
                                row_names_gp = grid::gpar(cex = max(min(1,30/nrow(mean_distances_matrix)),0.45)))

  h2 <- ComplexHeatmap::Heatmap(1-max_distances_matrix,
                                name = 'Max similarity',
                                col = grDevices::blues9,
                                rect_gp = grid::gpar(col = 'black'),
                                width = ncol(mean_distances_matrix)*grid::unit(constant, "mm"),
                                height = nrow(mean_distances_matrix)*grid::unit(constant, "mm"),
                                column_names_gp = grid::gpar(cex = max(min(1,30/nrow(mean_distances_matrix)),0.45)),
                                row_names_gp = grid::gpar(cex = max(min(1,30/nrow(mean_distances_matrix)),0.45)))

  h3 <- ComplexHeatmap::Heatmap(1-min_distances_matrix,
                                name = 'Min similarity',
                                col = grDevices::blues9,
                                rect_gp = grid::gpar(col = 'black'),
                                width = ncol(mean_distances_matrix)*grid::unit(constant, "mm"),
                                height = nrow(mean_distances_matrix)*grid::unit(constant, "mm"),
                                column_names_gp = grid::gpar(cex = max(min(1,30/nrow(mean_distances_matrix)),0.45)),
                                row_names_gp = grid::gpar(cex = max(min(1,30/nrow(mean_distances_matrix)),0.45)))

  return(S4Vectors::List('Mean_distance' = mean_distances_matrix,
                         'Max_distance' = max_distances_matrix,
                         'Min_distance' = min_distances_matrix,
                         'Mean_similarity_Heatmap' = h1,
                         'Max_similarity_Heatmap' = h2,
                         'Min_similarity_Heatmap' = h3))

}
