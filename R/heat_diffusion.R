#' heat_diffusion
#'
#' @param graph The network on which to diffuse
#' @param scores The starting scores to diffuse
#' @param alpha Alpha value to give to expm function
#' @param max_iter Maximum number of iterations
#' @inheritParams expm::expm
#'
#' @return The diffused scores
#' @export
#'
#' @examples
#' g <- igraph::make_ring(10)
#' igraph::E(g)$weight <- runif(igraph::ecount(g), 0.5, 1)
#' initial_scores <- rep(0, igraph::vcount(g))
#' initial_scores[1] <- 1
#' heat_diffusion(g, initial_scores)
heat_diffusion <- function(graph,
                           scores,
                           alpha = 0.7,
                           tol = 1e-6,
                           max_iter = 100) {

  A <- igraph::as_adjacency_matrix(graph, attr = "weight", sparse = TRUE)
  D <- Matrix::Diagonal(x = 1 / sqrt(rowSums(as.matrix(A))))
  L <- as.matrix(Matrix::Diagonal(nrow(A)) - D %*% A %*% D)
  f <- scores
  for (i in seq(1,max_iter)) {
    f_new <- expm::expm(-alpha * L) %*% scores
    if (sum(abs(f_new - f)) < tol) break
    f <- f_new
  }
  return(f[,1])

}
