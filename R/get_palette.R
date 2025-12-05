#' Get the palette
#'
#' @param net The reference network
#' @param color_attr Either SubClass or Stages (to get the palette of the SubClass
#' or Stages, respectively)
#'
#' @return The palette
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' subclass_palette <- c('A' = 'red', 'B' = 'blue', 'C' = 'green', 'D' = 'yellow')
#' net$SubClass_color <- subclass_palette[net$SubClass]
#' get_palette(net)
get_palette <- function(net, color_attr = 'SubClass') {
  if(color_attr == 'SubClass') {
    palette <- unique(net$SubClass_color)
    names(palette) <- unique(net$SubClass)
  }
  if(color_attr == 'Stages') {
    palette <- unique(net$Stages_color)
    names(palette) <- unique(net$Stages)
  }
  if(!color_attr %in% c('SubClass', 'Stages')) {
    message('color_attr must be either SubClass or Stages')
  }
  return(palette)
}
