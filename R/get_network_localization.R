#' Given an expression matrix and a snn graph in the form of an adjacency matrix
#' it returns the network localization value for each gene
#'
#' @param expression_matrix The expression matrix
#' @param snn_graph The adjacency matrix of the SNN graph (or other graphs)
#'
#' @return
#' The network localization values for each gene
#' @export
#'
#' @examples
#' snn_graph <- matrix(sample(seq(1, 100), 10000, replace = TRUE), ncol = 100)
#' colnames(snn_graph) <- paste0('Column', seq(1, 100))
#' rownames(snn_graph) <- paste0('Column', seq(1, 100))
#' expression_matrix <- matrix(sample(seq(1, 100), 100000, replace = TRUE),
#' ncol = 100)
#' colnames(expression_matrix) <- paste0('Column', seq(1, 100))
#' rownames(expression_matrix) <- paste0('Gene', seq(1, 1000))
#' get_network_localization(expression_matrix, snn_graph)
get_network_localization <- function(expression_matrix, snn_graph) {

  if(!methods::is(snn_graph, 'Matrix') & !methods::is(snn_graph, 'matrix')) {
    snn_graph <- as.matrix(igraph::as_adjacency_matrix(snn_graph, attr = 'weight'))
  }

  if(any(Matrix::rowSums(snn_graph) == 0)) {
    na_idx <- which(Matrix::rowSums(snn_graph) == 0)
    expression_matrix <- expression_matrix[,-na_idx]
    snn_graph <- snn_graph[-na_idx,-na_idx]
  }

  expression_matrix <- as.matrix(expression_matrix)

  weight_matrix <- sweep(snn_graph, 1, Matrix::rowSums(snn_graph), "/")

  standardized_expression <- t(scale(t(expression_matrix), center = TRUE, scale = TRUE))

  Hx <- apply(standardized_expression, 1, function(gene_expression) {
    sum(gene_expression * (weight_matrix %*% gene_expression))
  })

  return(Hx)

}
