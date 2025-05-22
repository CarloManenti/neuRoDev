#' Weighted summary
#'
#' @param v A numerical vector
#' @param weights A numerical vector of weights
#'
#' @return The weighted summary
#' @export
#'
#' @examples
#' set.seed(123)
#' weighted_summary(seq(1,100), weights=runif(100,0,1))
weighted_summary <- function(v,
                             weights) {

  # Description: given a numerical vector and a vector of weights,
  # it returns the weighted summary of the vector

  summar <- Hmisc::wtd.quantile(v,
                         weights = weights)
  w_mean <- stats::weighted.mean(v,
                          weights = weights)
  summar <- c(summar, w_mean)
  summar <- summar[c(1,2,3,6,4,5)]
  names(summar) <- c('Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.')
  return(summar)
}
