#' Builds the edges dataframe from an igraph network
#'
#' @param network The igraph network
#' @param layout A layout to define the X and Y coordinates of the edges
#' @param weight_quantile The quantile below which the edges will be put to NA
#' @param weights_normalization_coef A normalization coefficient to reduce the
#' width of the edges in the plots. It defaults to the number of nodes of the
#' network (with a maximum of 1000)
#'
#' @return The edges dataframe
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' M <- reference_signatures_correlation(S, refS)
#' network <- compute_correlation_network(M[,seq(3,dim(M)[2])])
#' umap_obj <- umap_graph_clustering(M)
#' build_edges_df(network, layout = umap_obj$umap_out$layout)
build_edges_df <- function(network,
                           layout,
                           weight_quantile=0.75,
                           weights_normalization_coef=NULL) {

  if(is.null(weights_normalization_coef)) {

    weights_normalization_coef <- min(1000, length(igraph::V(network)$name))

  }

  g <- igraph::as_data_frame(network)

  g$from.x <- layout[,1][match(g$from, rownames(layout))]
  g$from.y <- layout[,2][match(g$from, rownames(layout))]
  g$to.x <- layout[,1][match(g$to, rownames(layout))]
  g$to.y <- layout[,2][match(g$to, rownames(layout))]

  widths <- g$weight
  widths[which(widths < stats::quantile(widths, weight_quantile))] <- 0
  widths <- widths / max(widths)
  widths <- widths / weights_normalization_coef

  g$weight <- as.vector(widths)

  g$weight[which(g$weight == 0)] <- NA

  return(g)

}

