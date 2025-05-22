#' Ranks a vector (useful for signatures)
#'
#' @param signature A signature vector
#'
#' @return The ranks of the genes in the signature
#' @export
#'
#' @examples
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' rankedS <- apply(S, 2, rankSignatures)
#'
rankSignatures <- function(signature) {

  sortSig <- sort(signature)
  names_sig <- names(sortSig)
  rankedSig <- seq(1, length(sortSig))
  names(rankedSig) <- names_sig
  rankedSig <- rankedSig[names(signature)]
  return(rankedSig)

}
