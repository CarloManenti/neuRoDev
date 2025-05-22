#' Compute outlier score
#'
#' @param M A numerical matrix or data frame, like those returned by
#' `reference_signatures_correlation`
#' @param n_neighbors The number of neighbors to use to build the SNN network
#' @param new_clusters If this is used to evaluate newly mapped points, the
#' names of the newly mapped points.
#'
#' @return A score vector for each point, the lower the score, the more likely
#' that point is an outlier
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
#' o <- compute_outlier_score(M)
compute_outlier_score <- function(M,
                                  n_neighbors=15,
                                  new_clusters=NULL) {

  if(any(endsWith(tolower(colnames(M)), 'value'))) {
    sub_M <- M[,seq((max(which(endsWith(tolower(colnames(M)), 'value'))))+1,
                    dim(M)[2])]
  } else {
    sub_M <- M
  }

  to_remove_idxs <- which(apply(sub_M,
                                2,
                                function(i)
                                {any(is.na(suppressWarnings(as.numeric(i))))
                                }))

  if(length(to_remove_idxs) > 0) {
    to_remove_idxs <- seq(min(to_remove_idxs), dim(sub_M)[2])
    sub_M <- sub_M[,-to_remove_idxs]
  }

  M <- as.matrix(sub_M)

  dist_matrix <- as.matrix(stats::dist(M))

  snn_network <- compute_snn_network(signatures_cor = M,
                                     n_neighbors = n_neighbors)

  snn_matrix <- snn_network$sNN.Matrix

  if(!is.null(new_clusters)) {
    snn_matrix <- snn_matrix[new_clusters,which(!colnames(snn_matrix) %in% new_clusters)]
    snn_distance <- apply(snn_matrix, 1, function(i) {mean(sort(i, decreasing = TRUE)[seq(2, n_neighbors+1)])})
  } else {
    snn_distance <- apply(snn_matrix, 1, function(i) {mean(sort(i, decreasing = TRUE)[seq(2, n_neighbors+1)])})
  }

  return(snn_distance)

}
