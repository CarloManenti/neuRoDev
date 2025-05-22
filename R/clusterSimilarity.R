#' Computes the correlation of the signatures given defined groups
#'
#' @param signatures The signatures to consider
#' @param membership The membership vector that defines the groups
#'
#' @return The correlation matrix and the related Heatmap
#' @export
#'
#' @examples
#' #' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' membership <- c(rep('A', 3), rep('B', 2), rep('C', 2), rep('D', 2))
#' clusterSimilarity(S, membership)
#'
clusterSimilarity <- function(signatures,
                              membership) {

  groups <- unique(membership)

  similarity_matrix <- matrix(1, nrow = length(groups), ncol = length(groups))

  colnames(similarity_matrix) <- groups
  rownames(similarity_matrix) <- groups

  for(g1 in groups) {

    if(is.vector(signatures[,which(membership == g1)])) {
      sig_1 <- signatures[,which(membership == g1)]
    } else {
      sig_1 <- Matrix::rowMeans(signatures[,which(membership == g1)])
    }

    for(g2 in groups) {

      if((similarity_matrix[g1,g2] == 1) && (g1 != g2)) {

        if(is.vector(signatures[,which(membership == g2)])) {
          sig_2 <- signatures[,which(membership == g2)]
        } else {
          sig_2 <- Matrix::rowMeans(signatures[,which(membership == g2)])
        }

        cor_value <- stats::cor(sig_1, sig_2)
        similarity_matrix[g1,g2] <- cor_value
        similarity_matrix[g2,g1] <- cor_value
      }

    }

  }

  constant <- min(120/nrow(similarity_matrix), 5)

  h <- ComplexHeatmap::Heatmap(similarity_matrix,
                               name = 'Correlation',
                               rect_gp = grid::gpar('black'),
                               width = ncol(similarity_matrix)*grid::unit(constant, "mm"),
                               height = nrow(similarity_matrix)*grid::unit(constant, "mm"),
                               column_names_gp = grid::gpar(cex = max(min(1,30/nrow(similarity_matrix)),0.45)),
                               row_names_gp = grid::gpar(cex = max(min(1,30/nrow(similarity_matrix)),0.45)))

  return(S4Vectors::SimpleList('Matrix' = similarity_matrix, 'Heatmap' = h))

}
