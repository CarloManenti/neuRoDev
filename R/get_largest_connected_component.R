#' Get largest connected component
#'
#' @param Net An igraph Network
#'
#' @return The largest connected component of the netwokr
#' @export
#'
#' @examples
#' set.seed(123)
#' num_vertices <- 100
#' p <- 0.3
#' random_graph <- igraph::erdos.renyi.game(num_vertices, p, directed = FALSE)
#' get_largest_connected_component(random_graph)
get_largest_connected_component <- function(Net) {

  # Description: this function computes and returns the largest connected
  # component of a network Net

  D <- igraph::decompose(Net)
  o <- D[[which.max(unlist(lapply(D, igraph::vcount)))]]
  return(o)
}
