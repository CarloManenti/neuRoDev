#' Keeps only genes exclusively differentially expressed in each group
#'
#' @param pwOutList The output of get_pairwise_wilcoxauc filtered with filter_dge_out
#'
#' @return The genes for each group that are exclusive for that group
#' @export
#'
#' @examples
#' group <- c(rep(c(1,2,3), each = 3), 4)
#' out <- lapply(unique(group), function(i) {
#' o <- get_pairwise_wilcoxauc(expMat = matrix(runif(500,0,100),
#' ncol = 10), group = group, qG = i)
#' pDEGs <- lapply(o, filter_dge_out, cutFrac=0.01, cutFC=0.0, cutPadj=1)
#' return(pDEGs)
#' })
#'
#' keep_exclusive_DEGs(out)
keep_exclusive_DEGs <- function(pwOutList) {

  # Description: given a get_pairwise_wilcoxauc result and the BOG obtained
  # filtering with filter_dge_out, it returns the genes for each group that are
  # exclusive

    out <- lapply(pwOutList, function(j) {
        l <- list_to_assignment_matrix(lapply(j, function(i) {
            i[,'feature']
        }))

        if(dim(l)[1] == 0) {
            return(NULL)
        } else {
            names(which(Matrix::rowMeans(l)==1))
        }
    })
    return(out)
}
