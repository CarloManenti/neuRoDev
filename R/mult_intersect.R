#' Multiple intersection
#'
#' @param vec_list A vector list
#'
#' @return Common elements between all vectors in the list
#' @export
#'
#' @examples
#' v1 <- c('Z', 'A', 'B', 'D')
#' v2 <- c('A', 'B', 'C', 'D')
#' v3 <- c('B', 'Z', 'D', 'H')
#' mult_intersect(list(v1, v2, v3))
mult_intersect <- function(vec_list) {
  intersected <- vec_list[[1]]
  for(i in seq(2,length(vec_list))) {
    intersected <- intersect(intersected, vec_list[[i]])
  }
  return(intersected)
}
