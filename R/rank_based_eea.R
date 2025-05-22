#' Ranked based differential expression analysis
#'
#' It ranks rows in matrix M based on the number of non-zero
#' entries they have in members of a group in group (done for each group).
#' Then computes the differences in (transformed) ranks between each group
#' and the mean of the others (RelRanks) and the same done pairwise
#' (RelRanks_pairs). Then, given a Fold-Change threshold of FCtresh, it
#' returns the genes that have an higher Fold-Change than FCtresh compared
#' to all the other groups (done for each group) (Gs). Eventually, it also
#' computes the Fold-Change rowSum for each gene in each group
#' (RelRanks_pairs_sum).
#'
#' @param M An expression matrix
#' @param group A membership vector for the columns of M
#' @param FCtresh A threshold of the Fold-Change
#'
#' @return A list containing the fraction of non-zero values for each group
#' (get_column_group_detection_rate), the ranks of rows based on such fraction,
#' the relative ranks compared to other groups (RelRanks), and the same but when
#' done pairwise (RelRanks_pairs), the Fold-Change row sum for each gene in
#' each group (RelRanks_pairs_sum) and the genes with an higher Fold-Change than
#' FCtresh compared to all the other groups (Gs).
#' @export
#'
#' @examples
#' rank_based_eea(matrix(sample(seq(0,1),100,
#' replace=TRUE),
#' ncol=10),
#' group=c(rep(c(1,2,3),
#' each = 3),
#' 4))
rank_based_eea <- function(M,
                           group,
                           FCtresh=1) {

  # Description: it ranks rows in matrix M based on the number of non-zero
  # entries they have in members of a group in group (done for each group).
  # Then computes the differences in (transformed) ranks between each group
  # and the mean of the others (RelRanks) and the same done pairwise
  # (RelRanks_pairs). Then, given a Fold-Change threshold of FCtresh, it
  # returns the genes that have an higher Fold-Change than FCtresh compared
  # to all the other groups (done for each group) (Gs). Eventually, it also
  # computes the Fold-Change rowSum for each gene in each group
  # (RelRanks_pairs_sum).

  detection_rate <- get_column_group_detection_rate(M,
                                                    group)

  Ranks <- apply(detection_rate, 2, rank)
  Ranks <- Ranks/apply(Ranks, 2, max)

  RelRanks <- do.call(cbind, lapply(colnames(Ranks), function(i) log2(Ranks[,i])-log2(rowMeans(Ranks[,!colnames(Ranks)%in%i]))))
  colnames(RelRanks) <- colnames(Ranks)

  RelRanks_pairs <- lapply(colnames(Ranks), function(j) {
    temp <- do.call(cbind, lapply(colnames(Ranks)[!colnames(Ranks)%in%j],
                                  function(i) {
                                    log2(Ranks[,j]) - log2(Ranks[,i])
                                  }))
    colnames(temp) <- colnames(Ranks)[!colnames(Ranks)%in%j]
    return(temp)
  })

  names(RelRanks_pairs) <- colnames(Ranks)

  Gs <- lapply(RelRanks_pairs, function(i) {
    names(which(rowMeans(i>FCtresh)==1))
  })

  RelRanks_pairs_sum <- do.call(cbind, lapply(RelRanks_pairs, rowSums))

  o <- S4Vectors::List(
    celltype.detection.rate=detection_rate,
    ranks=Ranks,
    logFC_other=RelRanks,
    RelRanks_pairs=RelRanks_pairs,
    RelRanks_pairs_sum=RelRanks_pairs_sum,
    twoFold_increased_allpairs=Gs
  )
  return(o)
}
