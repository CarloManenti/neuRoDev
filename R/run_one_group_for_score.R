#' Perform the enrichment for one group at a time
#'
#' @param expression_matrix The expression matrix to use
#' @param groups The groups in which to devide the expression matrix
#' @param groupX The wanted group
#' @param Datatype The expression_matrix type. Defaults to RNAseq.
#' @param threshold The threshold for the pvalue above which to put the coef
#' returned by limma::eBayes to 0.
#'
#' @return
#' The sum of the coefficients for the comparison between groupX and the other
#' groups.
#' @export
#'
#' @examples
#' expression_matrix <- matrix(sample(seq(1,100), 10000, replace = TRUE), ncol = 10)
#' rownames(expression_matrix) <- paste0('Gene.', seq(1, nrow(expression_matrix)))
#' colnames(expression_matrix) <- paste0('Column.', seq(1, ncol(expression_matrix)))
#' groups <- paste0('Group.', c(1,1,2,3,3,4,5,4,2,5))
#' groupX <- 'Group.1'
#' run_one_group_for_score(expression_matrix, groups, groupX)
run_one_group_for_score <- function(expression_matrix,
                                    groups,
                                    groupX,
                                    Datatype="RNAseq",
                                    threshold = 0.1) {

  f <- factor(groups, levels=unique(groups))
  design <- stats::model.matrix(~0+f)
  colnames(design) <- unique(groups)
  mycontrasts <- limma::makeContrasts(contrasts=make_group_contrasts_for_one(groupX,groups),
                                      levels=design)
  if(Datatype=="RNAseq") {
    v <- limma::voom(expression_matrix, design)
    fit <- limma::lmFit(v, design)
  }
  if(Datatype=="Array") fit <- limma::lmFit(expression_matrix, design)

  fit2 <- limma::eBayes(limma::contrasts.fit(fit, mycontrasts))
  fit.test <- stats::p.adjust(fit2$p.value, method="bonferroni")
  select <- fit.test > threshold
  fit.coef <- fit2$coef
  fit.coef[select] <- 0
  return(apply(fit.coef, 1, sum))

}
