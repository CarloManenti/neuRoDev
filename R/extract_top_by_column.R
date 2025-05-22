#' Extract the top rownames (genes) given a matrix
#'
#' @param mat A numerical matrix
#' @param n Number of top rows
#' @param decreasing A boolean to tell if the top genes should be taken in
#' increasing or decreasing order. Defaults to TRUE.
#' @param vector If the top rows should be given as a single vector. Defaults
#' to TRUE.
#' @param unique Only if vector = TRUE. Boolean to tell if only the unique
#' top rows should be returned. Defaults to TRUE.
#'
#' @return The top rownames for each column in `mat`
#' @export
#'
#' @examples
#' mat <- matrix(sample(seq(1,10000), 1000), ncol = 10)
#' colnames(mat) <- paste0('Col', seq(1, ncol(mat)))
#' rownames(mat) <- paste0('Row', seq(1, nrow(mat)))
#' extract_top_by_column(mat, n = 2)
extract_top_by_column <- function(mat,
                                  n = 5,
                                  decreasing = TRUE,
                                  vector = TRUE,
                                  unique = TRUE) {

  selected_genes <- apply(mat, 2, function(i) {rownames(mat)[order(i,
                                                                   decreasing = decreasing)[seq(1,n)]]})

  if(vector) {
    selected_genes <- as.vector(selected_genes)
    if(unique) {
      selected_genes <- unique(selected_genes)
    }
  }

  return(selected_genes)

}
