#' Computes the highly specific genes based on signatures
#'
#' @param expMat An expression matrix
#'
#' @return For each column in expMat, returns the genes that have a signature value higher than 1
#' @export
#'
#' @examples
#' set.seed(123)
#' FC_signatures_exclusiveSets(matrix(runif(200,0,10), ncol = 10))
FC_signatures_exclusiveSets <- function(expMat) {

  # Description: it computes the exclusive set signature of all the columns of
  # an expression matrix expMat

  if(is.null(colnames(expMat))) {
    colnames(expMat) <- paste0('Column-', seq_len(dim(expMat)[2]))
  }

  out <- lapply(colnames(expMat), function(i) {
    FC_pairwise_xFCexclusive(expMat, QueryCol = i)
  })

  names(out) <- colnames(expMat)
  return(out)
}
