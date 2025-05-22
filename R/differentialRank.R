#' Computes differential rank of genes in a set of signatures
#'
#' @param signatures The signatures to evaluate
#' @param membership A membership vector of the groups to analyse. One value for
#' each column of the signatures
#' @param selected_groups A subgroup to evaluate
#' @param n_top_genes The number of most differentially ranked genes to return
#'
#' @return A list of the n_top_genes most differentially expressed in one group
#' versus all the others and in all pairwise comparisons
#' @export
#'
#' @examples
#' S <- FC_signatures(matrix(runif(600000,0,30), ncol = 60))
#' S[,seq(1,20)] <- S[,1]
#' S[,seq(21,40)] <- S[,21]
#' S[,seq(41,60)] <- S[,41]
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' differentialRank(signatures = S, membership = c(rep('Group_1', 20),
#' rep('Group_2', 20),
#' rep('Group_3', 20)),
#' selected_groups = c('Group_1', 'Group_3'))
differentialRank <- function(signatures,
                             membership,
                             selected_groups=NULL,
                             n_top_genes=100) {

  if(!is.null(selected_groups)) {
    signatures <- signatures[,which(membership %in% selected_groups)]
    membership <- membership[which(membership %in% selected_groups)]
  }

  ranks <- apply(signatures, 2, rankSignatures)

  pBOG <- get_pairwise_wilcoxauc_DEs(expMat = ranks,
                                     group = membership,
                                     filterDEset = FALSE)

  pBOG <- lapply(pBOG, function(i) {
    lapply(i, function(x) {
      x[seq(1, min(n_top_genes, dim(x)[1])),'feature']
    })
  })

  memberships <- sort(unique(membership))

  wilcox_dif_out <- lapply(memberships, function(i) {

    run_wilcox_differential(ranks,
                            binGroup = ifelse(membership==i, 1, 0),
                            orderOut = TRUE)

  })

  names(wilcox_dif_out) <- memberships

  BOG <- lapply(wilcox_dif_out, function(i) {

    rownames(filter_dge_out(i))

  })

  BOG <- lapply(BOG, function(i) {i[seq(1, min(n_top_genes, length(i)[1]))]})

  return(S4Vectors::SimpleList('Pairwise' = pBOG, 'One_vs_all' = BOG))

}
