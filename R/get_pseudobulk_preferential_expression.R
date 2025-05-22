#' Obtain the pseudobulk preferential expression
#'
#' @inheritParams run_one_group_for_score
#'
#' @return
#' The pseudobulk preferential expression scores
#' @export
#'
#' @examples
#' expression_matrix <- matrix(sample(seq(1,100), 10000, replace = TRUE), ncol = 10)
#' rownames(expression_matrix) <- paste0('Gene.', seq(1, nrow(expression_matrix)))
#' colnames(expression_matrix) <- paste0('Column.', seq(1, ncol(expression_matrix)))
#' groups <- paste0('Group.', c(1,1,2,3,3,4,5,4,2,5))
#' get_pseudobulk_preferential_expression(expression_matrix, groups)
get_pseudobulk_preferential_expression <- function(expression_matrix,
                                                   groups,
                                                   Datatype="RNAseq",
                                                   threshold = 0.1){

  groups <- make.names(groups)

  results <- matrix(NA,
                    nrow=dim(expression_matrix)[1],
                    ncol=length(unique(groups)))

  rownames(results) <- rownames(expression_matrix)
  colnames(results) <- unique(groups)

  for (groupX in unique(groups)){
    colX <- run_one_group_for_score(expression_matrix,
                                    groups,
                                    groupX,
                                    Datatype="RNAseq",
                                    threshold = threshold)
    results[,groupX] <- colX
  }

  return(results)

}
