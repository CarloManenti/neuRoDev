#' Computes the pairwise Wilcoxon test and filters the results if wanted
#'
#' @inheritParams get_pairwise_wilcoxauc
#' @param filterDEset A boolean variable that indicates if the function
#' keep_exclusive_DEGs has to be run, keeping only exclusive genes for each group
#'
#' @return Returns the genes for each group that are exclusive
#' @export
#'
#' @examples
#' set.seed(123)
#' get_pairwise_wilcoxauc_DEs(expMat = matrix(runif(500,0,100), ncol = 10),
#' c(rep(c(1,2,3), each = 3), 4))
get_pairwise_wilcoxauc_DEs <- function(expMat,
                                       group,
                                       filterDEset=TRUE) {

  # Description: it computes the pairwise Wilcox differential expression for
  # each unique group in group, given the expression matrix expMat.
  # If filteredDEset=TRUE, it returns the genes for each group that are
  # exclusive

  out <- lapply(unique(group), function(i) {
    o <- get_pairwise_wilcoxauc(expMat = expMat, group = group, qG = i)
    pDEGs <- lapply(o, filter_dge_out)
    return(pDEGs)
  })

  names(out) <- unique(group)
  if(filterDEset) {
    out <- keep_exclusive_DEGs(out)
  }

  return(out)
}
