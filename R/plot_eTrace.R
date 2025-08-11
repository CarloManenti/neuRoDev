#' Plot the eTraces
#'
#' @inheritParams get_eTrace
#' @param score A score to plot on the y axis
#' @param expression_enrichment A boolean, if TRUE the expression enrichment is
#' computed and plotted (`get_eTrace` function). Defaults to FALSE
#' @param main The title of the plot
#' @param upper_colors The colors for the upper plot
#' @param lower_colors The colors for the lower plot
#' @param ylab The ylab to add. Defaults to `score`
#'
#' @return The plot of the eTraces divided into two subplots
#' @export
#'
#' @examples
plot_eTrace <- function(net,
                        genes = NULL,
                        score = NULL,
                        expression_enrichment = FALSE,
                        main=NULL,
                        upper_colors = NULL,
                        lower_colors = NULL,
                        nRand = 100, ylab = "score") {

  if(is.null(upper_colors)) {
    upper_colors <- net$Stages_colors
  }

  if(is.null(lower_colors)) {
    lower_colors <- net$SubClass_colors
  }

  if(is.null(genes)) {
    genes <- rownames(net)
  }

  if(expression_enrichment) {
    eTrace <- get_eTrace(net = net, genes = genes, nRand = nRand)
    ylab <- "expression enrichment (z)"
  } else {
    if(is.null(score)) {
      if(length(genes) == 1) {
        eTrace <- S4Vectors::List(z = SingleCellExperiment::logcounts(net)[genes,])
      } else {
        eTrace <- S4Vectors::List(z = Matrix::colMeans(SingleCellExperiment::logcounts(net)[genes,]))
      }
      ylab <- 'logcounts'
    } else {
      if(is.matrix(score)) {
        if(length(genes) == 1) {
          eTrace <- S4Vectors::List(z = score[genes,])
        } else {
          eTrace <- S4Vectors::List(z = Matrix::colMeans(score[genes,]))
        }
      } else {
        eTrace <- S4Vectors::List(z = score)
      }
    }
  }

  col_zero_line <- ifelse(all(eTrace$z >= 0), NA, 'black')

  x <- 1:ncol(net)
  nat_idx <- which(net$Stages == '6-early_infancy')[1]
  y_idx <- min(eTrace$z)+abs(min(eTrace$z))*0.05

  graphics::par(mfrow=c(2,1))
  graphics::par(mar=c(0,5,2,2))
  plot(eTrace$z, pch=21, bg=upper_colors, main=main, ylab=ylab, xaxt = 'n')
  graphics::abline(h=0, lty=2, lwd=2, col = col_zero_line)
  graphics::abline(v=nat_idx, col = 'darkgrey', lwd = 2, lty = 2)
  graphics::text(x = nat_idx+5, y = y_idx, labels = 'postnatal', pos = 4)
  graphics::text(x = nat_idx-5, y = y_idx, labels = 'prenatal', pos = 2)
  graphics::lines(stats::smooth.spline(x,eTrace$z, spar = 1), col='red', lwd=2.5)
  graphics::par(mar=c(2.5,5,0.5,2))
  plot(eTrace$z, pch=21, bg=lower_colors, ylab=ylab)
  graphics::abline(h=0, lty=2, lwd=2, col = col_zero_line)
  graphics::abline(v=nat_idx, col = 'darkgrey', lwd = 2, lty = 2)
  graphics::text(x = nat_idx+5, y = y_idx, labels = 'postnatal', pos = 4)
  graphics::text(x = nat_idx-5, y = y_idx, labels = 'prenatal', pos = 2)
  graphics::lines(stats::smooth.spline(x,eTrace$z, spar = 1), col='red', lwd=2.5)
}
