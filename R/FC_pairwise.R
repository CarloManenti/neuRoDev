#' Computes the summed difference of a column to all the others in a matrix (signature)
#'
#' @param expMat An expression matrix
#' @param QueryCol A query column
#' @param RefColSet A subset of the matrix to consider in the difference with QueryCol
#' @param sortedFC A boolean variable to define if the result should be sorted (decreasing=TRUE)
#'
#' @return The summed difference between QueryCol and the other columns (the signature)
#' @export
#'
#' @examples
#' set.seed(123)
#' FC_pairwise(matrix(runif(200, 0, 10), ncol=10),
#' QueryCol=3,
#' RefColSet=c(1,4,5),
#' sortedFC=TRUE)
FC_pairwise <- function(expMat,
                        QueryCol,
                        RefColSet=NULL,
                        sortedFC=FALSE) {

  # Description: given an expression matrix expMat and a query column QueryCol,
  # it computes the difference between each other column of the
  # expression matrix and the query column, and then sums the differences
  # between them. If RefColSet you can define a subset of the expMat to
  # consider in the subtraction with the query column.
  # If sortedFC=TRUE, it sorts the result (decreasing=TRUE),
  # such that genes that are the most different will appear earlier.

  if(is.null(RefColSet)) {
    temp <- expMat[,QueryCol]-expMat
    colnames(temp) <- colnames(expMat)
  } else {
    temp <- expMat[,QueryCol]-expMat[,RefColSet]
    colnames(temp) <- RefColSet
  }

  rownames(temp) <- rownames(expMat)
  FC <- Matrix::rowSums(temp)

  if(sortedFC) {
    FC <- sort(FC, decreasing = TRUE)
  }

  return(FC)
}
