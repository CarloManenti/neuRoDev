#' Computes the genes that have a signature higher than a threshold
#'
#' @inheritParams FC_pairwise
#' @param xFC The difference threshold for a gene to be included
#'
#' @return The genes that pass the threshold
#' @export
#'
#' @examples
#' set.seed(123)
#' FC_pairwise_xFCexclusive(matrix(runif(200, 0, 10), ncol=10),
#' QueryCol=3,
#' RefColSet=c(1,4,5),
#' xFC=1.1)
FC_pairwise_xFCexclusive <- function(expMat,
                                     QueryCol,
                                     RefColSet=NULL,
                                     xFC=1) {

  # Description: given an expression matrix expMat and a query column QueryCol,
  # it computes the difference between each other column of the expression
  # matrix and the query column, and then sums the differences between them.
  # If RefColSet you can define a subset of the expMat to consider in the
  # subtraction with the query column. Then, it returns only the genes that
  # have a difference higher than xFC across all comparisons

  if(is.null(colnames(expMat))) {
    colnames(expMat) <- paste0('Column-', seq_len(dim(expMat)[2]))
    QueryCol <- paste0('Column-', QueryCol)
  }

  if(is.null(rownames(expMat))) {
    rownames(expMat) <- paste0('Rows-', seq_len(dim(expMat)[1]))
  }

  if(is.null(RefColSet)) {
    RefColSet <- colnames(expMat)[-which(colnames(expMat)==QueryCol)]
  }

  temp <- do.call(cbind, lapply(RefColSet, function(i) {
    expMat[,QueryCol]-expMat[,i]
  }))

  colnames(temp) <- RefColSet
  rownames(temp) <- rownames(expMat)
  temp <- names(which(rowMeans(temp>xFC)==1))
  return(temp)
}
