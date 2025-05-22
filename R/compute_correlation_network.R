#' Computes an igraph network from the correlation matrix
#'
#' @param signatures_cor The signatures correlations
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
#' compute_correlation_network(M)
compute_correlation_network <- function(signatures_cor) {

  correlation_distance <- as.matrix(stats::dist(signatures_cor))

  correlation_distance <- correlation_distance / max(correlation_distance)

  correlation_similarity <- 1 - correlation_distance

  diag(correlation_similarity) <- 0

  correlation_network <- igraph::graph_from_adjacency_matrix(correlation_similarity,
                                                             mode = 'undirected',
                                                             weighted = TRUE)

  return(correlation_network)

}
