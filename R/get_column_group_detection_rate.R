#' Defines the portion of members of each group that have non-zero values in each row
#'
#' @param M A matrix
#' @param group A membership vector for columns of M
#' @param doParallel A boolean variable to run the code in parallel
#'
#' @return The portion of members of each group that have non-zero value in each row
#' @export
#'
#' @examples
#' get_column_group_detection_rate(matrix(sample(seq(0,1),
#' 100,
#' replace=TRUE),
#' ncol=10),
#' group=c(rep(c(1,2,3), each = 3), 4))
get_column_group_detection_rate <- function(M,
                                            group,
                                            doParallel=FALSE) {

  # Description: given a matrix M it returns, for each group in group, the
  # portion of members of the group that have a non-zero value of each row.
  # It can be done in parallel (doParallel=TRUE).

  Ugroups <- sort(names(which(table(group)>1)))

  if(doParallel) {
    out <- do.call(cbind, parallel::mclapply(mc.cores = 8, X=Ugroups, function(i) {
      Matrix::rowMeans(methods::as(M[,group==i], "sparseMatrix")>0)}))
  }
  else {
    out <- do.call(cbind, lapply(Ugroups, function(i) {
      Matrix::rowMeans(methods::as(M[,group==i], "sparseMatrix")>0)
    }))
  }
  colnames(out) <- Ugroups
  rownames(out) <- rownames(M)
  return(out)
}
