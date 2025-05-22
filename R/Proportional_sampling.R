#' Proportional sampling
#'
#' @param M An expression matrix
#' @param membership_vector A membership vector for the columns of the matrix
#' @param ReturnSCE A boolean variable that defines if the subsampled matrix
#' is returned or if only the indexes
#' @param Ntotal The number of cells in total
#'
#' @return Either the matrix filtered or a list of indexes and labels
#' @export
#'
#' @examples
#' set.seed(123)
#' Proportional_sampling(matrix(runif(20000,0,10), ncol=200),
#' membership_vector = rep(c('A','B','C','D'),
#' each = 50),
#' Ntotal=20)
Proportional_sampling <- function(M,
                                  membership_vector,
                                  ReturnSCE=TRUE,
                                  Ntotal=2000) {

  Ntotal <- min(Ntotal, (dim(M)[2]-1))

  if(is.character(membership_vector)) {
    membership_vector <- as.numeric(as.factor(membership_vector))
  }

  Freqs <- table(membership_vector)/sum(membership_vector)
  X <- lapply(unique(membership_vector), function(i) which(membership_vector==i))
  names(X) <- unique(membership_vector)
  Ns <- round(Ntotal*vapply(X, length, integer(1))/ncol(M))

  X <- lapply(names(X), function(i) sample(X[[i]], Ns[i]))
  names(X) <- unique(membership_vector)
  Labs <- rep(names(X), vapply(X, length, integer(1)))
  x <- list(idx=as.numeric(unlist(X)), labs=Labs)

  if(ReturnSCE) {
    return(M[,x$idx])
  } else {
    return(x)
  }
}
