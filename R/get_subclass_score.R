#' Get subclass score for mapping
#'
#' @param mapped_obj The result of mapData
#' @param reference_df The dataframe with the annotations per reference cluster
#' @param sub_idxs The indexes from which to select the annotations scores in
#' mapped_obj$NearestNeighborsAnnotation$Annotations. If NULL (default), they are
#' all kept.
#'
#' @return A score per subclass considering all the mapped data points
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
#' get_subclass_score(mapped_obj, reference_df = annotated_M)
get_subclass_score <- function(mapped_obj, reference_df, sub_idxs = NULL) {

  annotations <- mapped_obj$NearestNeighborsAnnotation$Annotations

  if(is.null(sub_idxs)) {
    sub_idxs <- seq(1, length(annotations))
  }

  annotations <- annotations[sub_idxs]

  subclass_score <- rowsum(as.matrix(unlist(annotations)),
                           group = unlist(lapply(strsplit(names(unlist(annotations)), '.', fixed = T),
                                                 function(i) {paste0(i[seq(2, length(i))], collapse = '.')})))[names(table(unlist(lapply(strsplit(names(unlist(annotations)), '.', fixed = T),
                                                                                                 function(i) {paste0(i[seq(2, length(i))], collapse = '.')})))),1]/table(unlist(lapply(strsplit(names(unlist(annotations)), '.', fixed = T),
                                                                                                                                               function(i) {i[2]})))

  for(subclass in unique(reference_df$SubClass)) {

    if(!subclass %in% names(subclass_score)) {
      previous_names <- names(subclass_score)
      subclass_score <- c(subclass_score, 0)
      names(subclass_score) <- c(previous_names, subclass)
    }

  }

  subclass_score <- subclass_score[unique(reference_df$SubClass)]

  return(subclass_score)

}
