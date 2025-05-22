#' Given an expression matrix, it returns the specificity value for each gene
#'
#' @param expression_matrix An expression matrix
#'
#' @return The specificity values for each gene
#' @export
#'
#' @examples
#' get_specificity(matrix(sample(seq(1,100), 10000, replace = TRUE),
#' ncol = 100))
get_specificity <- function(expression_matrix) {
  expression_matrix <- as.matrix(expression_matrix)
  x.t <- ncol(expression_matrix)
  x.sum <- Matrix::colSums(expression_matrix, na.rm=TRUE)
  x.p <- t(t(expression_matrix)/x.sum)
  p.gen <- Matrix::rowMeans(x.p,na.rm=TRUE)
  Si <- (apply((x.p/p.gen)*log2(x.p/p.gen),1,sum,na.rm=TRUE))/x.t
  return(Si)
}
