#' Get the expression enrichment Trace
#'
#' @param net The reference network
#' @param genes Genes of which to compute the expression enrichment
#' @param nRand The number of random sampling to use (defaults to 100)
#'
#' @return An S4Vectors List with the averaged scaled logcounts of the genes (`obs`),
#' the z score of the enrichment (`z`), and the pvalue (`p`)
#' @export
#'
#' @examples
get_eTrace <- function(net,
                       genes,
                       nRand=100) {

  genes <- genes[genes%in%rownames(net)]
  multiple <- length(genes) > 1
  if(multiple) {
    y <- Matrix::colMeans(t(scale(t(as.matrix(logcounts(net)[genes,])))))
  } else {
    v <- SingleCellExperiment::logcounts(net)[genes,]
    y <- (v-mean(v))/(stats::sd(v))
  }
  RR <- do.call(rbind, lapply(1:nRand, function(i) {
    if(multiple) {
      Matrix::colMeans(t(scale(t(as.matrix(SingleCellExperiment::logcounts(net)[sample(rownames(net), length(genes)),])))))
    } else {
      v <- SingleCellExperiment::logcounts(net)[sample(rownames(net), length(genes)),]
      (v-mean(v))/stats::sd(v)
    }
  }))
  RR <- RR[which(apply(RR, 1, function(i) {!any(is.na(i))})),]
  z <- (y-Matrix::colMeans(RR))/MatrixGenerics::colSds(RR)
  p <- 2 * stats::pnorm(-abs(z))

  return(S4Vectors::List(obs=y,
                         z=z,
                         p=p))

}
