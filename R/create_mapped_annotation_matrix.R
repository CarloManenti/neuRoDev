#' Get mapping annotation matrix
#'
#' @param mapped_obj The result of mapData
#' @param sub_idxs The indexes from which to select the annotations scores in
#' mapped_obj$NearestNeighborsAnnotation$Annotations. If NULL (default), they are
#' all kept.
#'
#' @return A matrix with the scores per annotation and for each mapped point
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' annotated_M <- reference_signatures_correlation(S, refS)
#' new_clusterS <- FC_signatures(matrix(runif(80,0,10), ncol = 4))
#' rownames(new_clusterS) <- paste0('Gene-', seq(1, dim(new_clusterS)[1]))
#' colnames(new_clusterS) <- paste0('New-', seq(1, dim(new_clusterS)[2]))
#' mapped_obj <- mapData(reference_df = annotated_M,
#' new_signatures = new_clusterS,
#' reference_signatures = refS,
#' color_attr = 'Best.Assignment')
#' create_mapped_annotation_matrix(mapped_obj)
create_mapped_annotation_matrix <- function(mapped_obj, sub_idxs = NULL) {

  annotations <- mapped_obj$NearestNeighborsAnnotation$Annotations

  if(is.null(sub_idxs)) {
    sub_idxs <- seq(1, length(annotations))
  }

  annotations <- annotations[sub_idxs]

  all_names <- unique(unlist(lapply(annotations, names)))

  mat <- sapply(annotations, function(vec) {
    out <- setNames(numeric(length(all_names)), all_names)
    out[names(vec)] <- vec
    out
  })

  return(mat)

}
