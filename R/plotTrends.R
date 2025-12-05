#' plotTrends
#'
#' @param net The reference network
#' @param pref_exp_genes The preferentially expressed genes list
#' @param sce The external object in the form of a SingleCellExperiment object
#' @param profiles The actual expression matrix used for the mapping
#' @param coldata Where to find the information of the replicates used for
#' averaging in the `colData` of `sce`
#' @param subclass If indicated as a vector it selects only a subset of
#' reference subclasses. Defaults to all subclasses
#' @param ylim The lower and upper limits of the y-axis
#' (defaults to NULL, which goes back to the automatic R plot ylim)
#' @param together Defines if the different preferentially expressed genes
#' should be plotted together (adjust `ylim` accordingly) (defaults to FALSE)
#' @param relative If TRUE, the expression is scaled from 0 to 1
#' (defaults to FALSE)
#'
#'
#' @return A plot of the expression of preferentially expressed genes in the
#' samples of `sce`
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' subclass_palette <- c('A' = 'red', 'B' = 'blue', 'C' = 'green', 'D' = 'yellow')
#' net$SubClass_color <- subclass_palette[net$SubClass]
#' net$X_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
#' net$Y_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
#' edges_from <- sample(colnames(net), size = 200, replace = TRUE)
#' edges_to <- sample(colnames(net), size = 200, replace = TRUE)
#' edges_from_x <- net$X_coord[match(edges_from, colnames(net))]
#' edges_from_y <- net$Y_coord[match(edges_from, colnames(net))]
#' edges_to_x <- net$X_coord[match(edges_to, colnames(net))]
#' edges_to_y <- net$Y_coord[match(edges_to, colnames(net))]
#' edges_weight <- sample(seq(0,1, length.out=1000), length(edges_from), replace = TRUE)
#' edges_df <- data.frame('from' = edges_from, 'to' = edges_to, 'weight' = edges_weight,
#' 'from.x' = edges_from_x,
#' 'from.y' = edges_from_y,
#' 'to.x' = edges_to_x,
#' 'to.y' = edges_to_y)
#' net@metadata$network$edges <- edges_df
#' pref_exp_genes <- list('A' = c('Gene-1', 'Gene-2'), 'B' = c('Gene-3', 'Gene-4'), 'C' = c('Gene-5', 'Gene-6'), 'D' = c('Gene-7', 'Gene-8'))
#' random_m <- matrix(sample(seq(1,10, length.out=10000), 15000*20, replace = TRUE), ncol = 20)
#' rownames(random_m) <- paste0('Gene-', seq(1,15000))
#' colnames(random_m) <- paste0('Sample-', seq(1,20))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = random_m))
#' sce$replicates <- rep(c('Rep1', 'Rep2', 'Rep3', 'Rep4'), each = 5)
#' profiles <- neuRoDev:::get_column_group_average(SingleCellExperiment::logcounts(sce), sce$replicates)
#' coldata <- 'replicates'
#' plotTrends(net, pref_exp_genes, sce, profiles, coldata)
plotTrends <- function(net,
                       pref_exp_genes,
                       sce,
                       profiles,
                       coldata,
                       subclass = NULL,
                       ylim = NULL,
                       together = FALSE,
                       relative = FALSE) {

  if(is.null(subclass)) {
    subclass <- unique(net$SubClass)
  }

  if(!together) {
    lapply(subclass, function(i) {
      g <- unlist(pref_exp_genes[names(pref_exp_genes) == i])
      g <- g[g %in% rownames(profiles)]
      if(length(g) < 2) {
        return(NULL)
      }

      y_means <- Matrix::colMeans(profiles[g, , drop = FALSE])
      single_values <- Matrix::colMeans(SingleCellExperiment::logcounts(sce)[g,])

      if(relative) {
        norm_factor <- max(y_means)
      } else {
        norm_factor <- 1
      }

      y_means <- y_means/norm_factor
      single_values <- single_values/norm_factor

      days <- colnames(profiles)
      y_sds <- sapply(days, function(d) {
        vals <- Matrix::colMeans(SingleCellExperiment::logcounts(sce)[g, SingleCellExperiment::colData(sce)[,coldata] == d,
                                                                      drop = FALSE])/norm_factor
        stats::sd(vals)
      })
      y_sds[which(is.na(y_sds))] <- 0

      ylim_range <- range(c(y_means - y_sds, y_means + y_sds))
      ylim_range_two <- range(single_values)
      ylim_range <- c(min(c(ylim_range[1], ylim_range_two[1])), max(c(ylim_range[2], ylim_range_two[2])))
      if(!is.null(ylim)) {
        ylim_range <- ylim
      }

      plot(seq_along(y_means), y_means,
           ylab = ifelse(relative, 'relative expression', 'expression'), xlab = '', xaxt = 'n',
           type = 'l', lwd = 2,
           col = net$SubClass_color[which(net$SubClass == i)][1],
           main = i, ylim = ylim_range)
      graphics::axis(side = 1, at = seq_along(y_means), labels = days, las = 2)

      x_vals <- seq_along(y_means)
      upper <- y_means + y_sds
      lower <- y_means - y_sds
      graphics::polygon(c(x_vals, rev(x_vals)),
              c(upper, rev(lower)),
              col = grDevices::adjustcolor(net$SubClass_color[which(net$SubClass == i)][1],
                                           alpha.f = 0.2),
              border = NA)

      graphics::lines(x_vals, y_means, lwd = 2,
            col = net$SubClass_color[which(net$SubClass == i)][1])

      for (j in seq_along(y_means)) {
        d <- names(y_means)[j]
        y_vals <- single_values[which(SingleCellExperiment::colData(sce)[,coldata] == d)]
        x_jitter <- jitter(rep(j, length(y_vals)), amount = 0.1)
        x_jitter[1] <- j
        show_top <- length(x_jitter)-2
        if(length(x_jitter) < 2) {
          idxs <- seq(1,length(x_jitter))
        } else {
          idxs <- c(which.min(y_vals), which.max(y_vals), sample(seq(1,length(x_jitter)), show_top))
        }
        x_jitter <- x_jitter[idxs]
        y_vals <- y_vals[idxs]
        graphics::points(x_jitter, y_vals, pch = 19, col = grDevices::adjustcolor(net$SubClass_color[which(net$SubClass == i)][1], alpha.f = 0.5), cex = 0.5)
      }

    })
  } else {

    days <- colnames(profiles)

    all_y_means <- lapply(subclass, function(i) {
      g <- unlist(pref_exp_genes[names(pref_exp_genes) == i])
      g <- g[g %in% rownames(profiles)]
      if(length(g) < 2) {
        return(NULL)
      }
      y_means <- Matrix::colMeans(profiles[g, , drop = FALSE])
      single_values <- Matrix::colMeans(SingleCellExperiment::logcounts(sce)[g,])

      if(relative) {
        norm_factor <- max(y_means)
      } else {
        norm_factor <- 1
      }

      y_means <- y_means/norm_factor
      single_values <- single_values/norm_factor

      y_sds <- sapply(days, function(d) {
        vals <- Matrix::colMeans(SingleCellExperiment::logcounts(sce)[g, SingleCellExperiment::colData(sce)[,coldata] == d, drop = FALSE])/norm_factor
        stats::sd(vals)
      })
      y_sds[which(is.na(y_sds))] <- 0
      return(c(y_means - y_sds, y_means + y_sds))

    })

    g <- unlist(pref_exp_genes[names(pref_exp_genes) == subclass[1]])
    g <- g[g %in% rownames(profiles)]
    if(length(g) < 2) {
      return(NULL)
    }

    y_means <- Matrix::colMeans(profiles[g, , drop = FALSE])
    single_values <- Matrix::colMeans(SingleCellExperiment::logcounts(sce)[g,])

    if(relative) {
      norm_factor <- max(y_means)
    } else {
      norm_factor <- 1
    }

    y_means <- y_means/norm_factor
    single_values <- single_values/norm_factor

    days <- colnames(profiles)
    y_sds <- sapply(days, function(d) {
      vals <- Matrix::colMeans(SingleCellExperiment::logcounts(sce)[g, SingleCellExperiment::colData(sce)[,coldata] == d, drop = FALSE])/norm_factor
      stats::sd(vals)
    })
    y_sds[which(is.na(y_sds))] <- 0

    ylim_range <- range(c(y_means - y_sds, y_means + y_sds))
    ylim_range_two <- range(single_values)
    ylim_range <- c(min(c(ylim_range[1], ylim_range_two[1])), max(c(ylim_range[2], ylim_range_two[2])))
    if(!is.null(ylim)) {
      ylim_range <- ylim
    }

    plot(seq_along(y_means), y_means,
         ylab = ifelse(relative, 'relative expression', 'expression'), xlab = '', xaxt = 'n',
         type = 'l', lwd = 2,
         col = net$SubClass_color[which(net$SubClass == subclass[1])][1],
         main = paste0(subclass, collapse = '-'), ylim = ylim_range)
    graphics::axis(side = 1, at = seq_along(y_means), labels = days, las = 2)

    x_vals <- seq_along(y_means)
    upper <- y_means + y_sds
    lower <- y_means - y_sds
    graphics::polygon(c(x_vals, rev(x_vals)),
            c(upper, rev(lower)),
            col = grDevices::adjustcolor(net$SubClass_color[which(net$SubClass == subclass[1])][1], alpha.f = 0.2),
            border = NA)

    graphics::lines(x_vals, y_means, lwd = 2,
          col = net$SubClass_color[which(net$SubClass == subclass[1])][1])

    for (j in seq_along(y_means)) {
      d <- names(y_means)[j]
      y_vals <- single_values[which(SingleCellExperiment::colData(sce)[,coldata] == d)]
      x_jitter <- jitter(rep(j, length(y_vals)), amount = 0.1)
      x_jitter[1] <- j
      show_top <- length(x_jitter)-2
      if(length(x_jitter) < 2) {
        idxs <- seq(1,length(x_jitter))
      } else {
        idxs <- c(which.min(y_vals), which.max(y_vals), sample(seq(1,length(x_jitter)), show_top))
      }
      x_jitter <- x_jitter[idxs]
      y_vals <- y_vals[idxs]
      graphics::points(x_jitter, y_vals, pch = 19, col = grDevices::adjustcolor(net$SubClass_color[which(net$SubClass == subclass[1])][1], alpha.f = 0.5), cex = 0.5)
    }

    for(i in subclass[seq(2,length(subclass))]) {
      g <- unlist(pref_exp_genes[names(pref_exp_genes) == i])
      g <- g[g %in% rownames(profiles)]
      if(length(g) < 2) {
        return(NULL)
      }

      y_means <- Matrix::colMeans(profiles[g, , drop = FALSE])
      single_values <- Matrix::colMeans(logcounts(sce)[g,])

      if(relative) {
        norm_factor <- max(y_means)
      } else {
        norm_factor <- 1
      }

      y_means <- y_means/norm_factor
      single_values <- single_values/norm_factor

      days <- colnames(profiles)
      y_sds <- sapply(days, function(d) {
        vals <- Matrix::colMeans(SingleCellExperiment::logcounts(sce)[g, SingleCellExperiment::colData(sce)[,coldata] == d, drop = FALSE])/norm_factor
        sd(vals)
      })
      y_sds[which(is.na(y_sds))] <- 0

      ylim_range <- range(c(y_means - y_sds, y_means + y_sds))
      ylim_range_two <- range(single_values)
      ylim_range <- c(min(c(ylim_range[1], ylim_range_two[1])), max(c(ylim_range[2], ylim_range_two[2])))
      if(!is.null(ylim)) {
        ylim_range <- ylim
      }

      graphics::lines(y_means,
            lwd = 2,
            col = net$SubClass_color[which(net$SubClass == i)][1])

      x_vals <- seq_along(y_means)
      upper <- y_means + y_sds
      lower <- y_means - y_sds
      graphics::polygon(c(x_vals, rev(x_vals)),
              c(upper, rev(lower)),
              col = grDevices::adjustcolor(net$SubClass_color[which(net$SubClass == i)][1], alpha.f = 0.2),
              border = NA)

      graphics::lines(x_vals, y_means, lwd = 2,
            col = net$SubClass_color[which(net$SubClass == i)][1])

      for (j in seq_along(y_means)) {
        d <- names(y_means)[j]
        y_vals <- single_values[which(SingleCellExperiment::colData(sce)[,coldata] == d)]
        x_jitter <- jitter(rep(j, length(y_vals)), amount = 0.1)
        x_jitter[1] <- j
        show_top <- length(x_jitter)-2
        if(length(x_jitter) < 2) {
          idxs <- seq(1,length(x_jitter))
        } else {
          idxs <- c(which.min(y_vals), which.max(y_vals), sample(seq(1,length(x_jitter)), show_top))
        }
        x_jitter <- x_jitter[idxs]
        y_vals <- y_vals[idxs]
        graphics::points(x_jitter, y_vals, pch = 19, col = grDevices::adjustcolor(net$SubClass_color[which(net$SubClass == i)][1], alpha.f = 0.5), cex = 0.5)
      }
    }

  }
}
