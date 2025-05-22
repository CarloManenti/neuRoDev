#' Computes the signatures of an expression matrix
#'
#' @param expMat An expression matrix
#'
#' @return The signatures for each of the columns of the matrix
#' @export
#'
#' @examples
#' set.seed(123)
#' FC_signatures(matrix(runif(200,0,10), ncol = 10))
FC_signatures <- function(expMat) {

  # Description: it computes the signature of all the columns of an expression
  # matrix expMat
  if(is.null(colnames(expMat))) {
    colnames(expMat) <- paste0('Column-', seq_len(dim(expMat)[2]))
  }
  out <- expMat - Matrix::rowMeans(expMat)

  return(out)
}
