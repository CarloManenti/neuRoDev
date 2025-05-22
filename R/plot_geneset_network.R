#' Plot the annotation network
#'
#' @param net The annotation network
#' @param quantile_to_remove Quantile of the edges weight distribution under
#' which edges are removed
#' @param increase The increase coefficient for the size of the edges
#' @param remove_unconnected A boolean, if TRUE unconnected nodes (after
#' edges removal) are removed.
#'
#' @return The igraph plot of the network
#' @export
#'
#' @examples
#' net <- igraph::sample_k_regular(10, 3, directed = FALSE, multiple = FALSE)
#' igraph::E(net)$weight <- sample(seq(1,10), length(igraph::E(net)),
#' replace = TRUE)
#' igraph::V(net)$size <- sample(seq(1,10), length(igraph::V(net)),
#' replace = TRUE)
#' plot_geneset_network(net)
plot_geneset_network <- function(net,
                           quantile_to_remove= NULL,
                           increase = 1,
                           remove_unconnected = TRUE) {

  if(is.null(quantile_to_remove)) {
    quantile_to_remove <- 1-(pathviewr::find_curve_elbow(data.frame('X' = seq(1, length(igraph::E(net)$weight)),
                                                                    'Y' = sort(igraph::E(net)$weight, TRUE)))/length(igraph::E(net)$weight))
  }

  threshold <- stats::quantile(igraph::E(net)$weight, quantile_to_remove)

  edges_to_remove <- igraph::E(net)[igraph::E(net)$weight < threshold]

  net <- igraph::delete_edges(net, edges_to_remove)

  igraph::E(net)$weight <- ((igraph::E(net)$weight)**2)*2

  igraph::E(net)$weight <- igraph::E(net)$weight*increase

  if(remove_unconnected) {
    net <- igraph::delete_vertices(net, igraph::V(net)[igraph::degree(net) == 0])
  }

  igraph::plot.igraph(net,
                      vertex.label.cex = ((igraph::V(net)$size)**(1/2)/2)/max((igraph::V(net)$size)**(1/2)/2),
                      vertex.label.dist = (igraph::V(net)$size)**(1/2)/4,
                      vertex.label.degree = -pi/2,
                      vertex.label.color = 'black',
                      edge.width = igraph::E(net)$weight)

}
