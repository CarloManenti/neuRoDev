#' Computes an igraph network from the correlation matrix
#'
#' @param signatures_cor The signatures correlations
#' @param n_neighbors The number of neighbors to consider in the sNN
#'
#' @return An igraph network
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
#' compute_snn_network(M, n_neighbors = 10)
compute_snn_network <- function(signatures_cor,
                                n_neighbors) {

  n_neighbors <- min(nrow(signatures_cor), n_neighbors)

  dist_matrix <- as.matrix(stats::dist(signatures_cor))

  nearest_neighbors <- t(apply(dist_matrix, 1, order)[seq(1,n_neighbors),])

  shared_neighbors_matrix <- matrix(0, nrow = nrow(signatures_cor),
                                    ncol = nrow(signatures_cor))

  for (i in seq(1, nrow(signatures_cor))) {
    for (j in seq(i, nrow(signatures_cor))) {

      shared_v <- intersect(nearest_neighbors[i,],
                            nearest_neighbors[j,])

      shared_count <- length(shared_v)

      if(shared_count != 0) {
        rank_i <- match(shared_v, nearest_neighbors[i,])
        rank_j <- match(shared_v, nearest_neighbors[j,])
        rank_diff <- abs(rank_i-rank_j)
        max_rank_diff <- n_neighbors-1
        if(all(rank_diff == 0)) {
          shared_count <- shared_count + 0.999
        } else if(all(rank_diff == max_rank_diff)) {
          shared_count <- shared_count - 0.999
        } else {
          shared_count <- shared_count + mean(0.999 - 1.998 * (rank_diff/max_rank_diff))
        }
      }

      shared_neighbors_matrix[i, j] <- shared_count

      shared_neighbors_matrix[j, i] <- shared_count

    }
  }

  shared_neighbors_matrix <- shared_neighbors_matrix/(n_neighbors+0.999)

  rownames(shared_neighbors_matrix) <- rownames(signatures_cor)
  colnames(shared_neighbors_matrix) <- rownames(signatures_cor)

  shared_neighbors_matrix_no_diag <- shared_neighbors_matrix
  diag(shared_neighbors_matrix_no_diag) <- 0

  snn_network <- igraph::graph_from_adjacency_matrix(shared_neighbors_matrix_no_diag,
                                                 mode = 'undirected',
                                                 weighted = TRUE)
  return(S4Vectors::List('network' = snn_network,
                         'sNN.Matrix' = shared_neighbors_matrix))

}
