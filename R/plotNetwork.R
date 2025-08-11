#' Plot the network
#'
#' @inheritParams plotSameLayout
#' @param show_edges A boolean, if FALSE it doesn't show the edges. Defaults to
#' TRUE
#' @param ... Additional parameters for `plot_net`
#'
#' @return A ggplot of the network
#' @export
#'
#' @examples
plotNetwork <- function(net,
                        color_attr = 'SubClass',
                        label_attr = 'SubClass',
                        col_vector = NULL,
                        no_label = FALSE,
                        show_edges = TRUE,
                        max_size = 3,
                        ...) {

  if(is.null(col_vector)) {
    col_vector <- paste0(color_attr, '_colors')
  }

  if(is.character(label_attr) && length(label_attr) == 1 &&  label_attr %in% colnames(SingleCellExperiment::colData(net))) {
    label_attr <- SingleCellExperiment::colData(net)[,label_attr]
  }

  if(is.character(color_attr) && length(color_attr) == 1 &&  color_attr %in% colnames(SingleCellExperiment::colData(net))) {
    color_attr <- SingleCellExperiment::colData(net)[,color_attr]
  }

  color_attr <- as.vector(color_attr)
  label_attr <- as.vector(label_attr)

  if(is.character(col_vector) && length(col_vector) == 1 &&  col_vector %in% colnames(SingleCellExperiment::colData(net))) {
    col_vector <- SingleCellExperiment::colData(net)[,col_vector]
  }

  if(no_label) {
    label_attr <- NULL
  }

  if(length(col_vector) == 1) {
    if(is.character(col_vector)) {
      col_vector <- Polychrome::createPalette(N = length(unique(color_attr))+2, seedcolors = c('#FF0000', '#FFFFFF'))
      col_vector <- col_vector[seq(2,length(col_vector))]
      names(col_vector) <- unique(color_attr)
      col_vector <- col_vector[color_attr]
    } else {
      col_vector <- rep(col_vector, length(color_attr))
      names(col_vector) <- color_attr
    }
  }

  layout <- cbind(net$X_coord, net$Y_coord)
  rownames(layout) <- colnames(net)
  layout <- cbind(layout, Size = max(max_size, 100/dim(layout)[1]))

  edges <- net@metadata$network$edges
  if(!show_edges) {
    edges <- NULL
  }

  layout <- cbind(layout, Colors = col_vector[as.character(color_attr)])
  layout <- cbind(layout, Group = color_attr)
  if(!is.null(label_attr)) {
    layout <- cbind(layout, SubGroup = label_attr)
  }
  ord <- order(match(layout[,"Group"], names(col_vector)))
  layout <- layout[ord,]

  p <- plot_net(layout,
                color_attr = layout[,"Group"],
                label_attr = label_attr[ord],
                edges = edges,
                show_edges = show_edges,
                ...)

  return(p)

}
