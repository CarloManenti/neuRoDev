#' Computes the mean rank of all genes across signatures
#'
#' @param signatures The wanted signatures
#' @param n_genes The number of top genes to be highlighted
#' @param weights Possible weights to have a weighted mean
#' @param new_cluster Possible indication of a new cluster whose ranks will be
#' averaged against the average of all other clusters (so will matter more)
#'
#' @return A SimpleList of the selected top genes and the mean rank of all genes
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' getSignatureMeans(S, n_genes = 10)
getSignatureMeans <- function(signatures,
                              n_genes,
                              weights=NULL,
                              new_cluster=NULL) {

  ranks <- apply(signatures, 2, rankSignatures)

  if(!is.null(new_cluster)) {
    new_sig <- signatures[,new_cluster]
    new_ranks <- ranks[,new_cluster]

    other_sig <- signatures[,which(colnames(signatures) != new_cluster)]
    other_ranks <- ranks[,which(colnames(ranks) != new_cluster)]

    if(!is.null(weights)) {
      means_other <- apply(other_ranks, 1, function(i) {stats::weighted.mean(i, w = weights)})
    } else {
      means_other <- Matrix::rowMeans(other_ranks)
    }

    other_new_ranks <- rankSignatures(means_other)

    means <- Matrix::rowMeans(cbind(new_ranks, other_new_ranks))

    means <- sort(means, decreasing=TRUE)

  } else {
    if(!is.null(weights)) {
      means <- apply(ranks, 1, function(i) {stats::weighted.mean(i, w = weights)})
    } else {
      means <- Matrix::rowMeans(ranks)
    }
  }

  means <- sort(means, decreasing=TRUE)

  return(S4Vectors::SimpleList('Selected_genes' = names(means)[seq(1,n_genes)], 'All_ranks' = means))

}
