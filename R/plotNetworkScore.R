#' Plot a score onto the network
#'
#' @inheritParams plotNetwork
#' @param genes The genes to consider for the score. If multiple genes are given,
#' the score will be their average
#' @param score The score to plot. If NULL (default), the logcounts of the
#' pseudobulks are used
#' @param expression_enrichment A boolean, if TRUE and score is NULL, the
#' expression enrichment of the given genes is computed. Defaults to FALSE.
#' @param smooth A boolean, if TRUE (default), the scores are smoothed to consider
#' the network structure
#' @param n_nearest The number of nearest neighbors to consider for smoothing.
#' Defaults to 15
#' @param na.vec The vector to define if certain clusters are not to be shown
#' @param palette The color palette to use. Defaults to `blues9`, unless
#' expression_enrichment = TRUE, in which it defaults to `BlGrRd` (if not
#' otherwise specified by the user).
#' @param main The title of the plot
#' @param show_edges A boolean, if FALSE no edges are shown. Defaults to TRUE
#' @param stroke The stroke of the points, defaults to 0.5
#' @param fix_alpha A boolean, if TRUE (default) all points will have alpha 1,
#' otherwise it will be proportional to the score
#'
#' @return A ggplot of the network colored based on the score
#' @export
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
#' plotNetworkScore(net, 'Gene-1', smooth = FALSE)
plotNetworkScore <- function(net,
                             genes = NULL,
                             score = NULL,
                             expression_enrichment = FALSE,
                             label_attr = NULL,
                             smooth = TRUE,
                             n_nearest = 15,
                             na.vec = NULL,
                             palette = NULL,
                             main=NULL,
                             show_edges = TRUE,
                             stroke = 0.5,
                             fix_alpha = TRUE,
                             max_size = 3) {

  if(is.null(palette)) {
    if(expression_enrichment) {
      palette <- 'BlGrRd'
    } else {
      palette <- 'blues9'
    }
  }

  if (palette %in% c("greys", "inferno", "magma", "viridis", "BlGrRd", "RdYlBu", "Spectral", "blues9")) {
    palette <- switch(palette,
                      greys = grDevices::gray.colors(100),
                      inferno = viridis::inferno(500, alpha = 0.8),
                      magma = viridis::magma(500,alpha = 0.8),
                      viridis = viridis::viridis(500, alpha = 0.8),
                      BlGrRd = (grDevices::colorRampPalette(c("blue", "grey", "red")))(500),
                      Spectral = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral"))))(100),
                      RdYlBu = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100),
                      blues9 = grDevices::blues9)
  } else {
    palette <- (grDevices::colorRampPalette(palette))(500)
  }

  if(is.character(label_attr) && length(label_attr) == 1 && label_attr %in% colnames(SingleCellExperiment::colData(net))) {
    label_attr <- SingleCellExperiment::colData(net)[,label_attr]
  }

  label_attr <- as.vector(label_attr)

  if(is.null(genes)) {
    genes <- rownames(net)
  } else {
    if(!all(genes %in% rownames(net))) {
      message(paste0(genes[which(!genes %in% rownames(net))], collapse = '-'), ' not in `net`')
      genes <- genes[which(genes %in% rownames(net))]
      if(length(genes) == 0) {
        return('No gene in `net`')
      }
    }
    genes <- genes[which(genes %in% rownames(net))]
  }

  if(is.null(score)) {
    if(expression_enrichment) {
      score <- get_eTrace(net = net, genes = genes)$z
    } else {
      exp_mat <- SingleCellExperiment::logcounts(net)
      score <- exp_mat[genes,]
    }
  } else {
    if(is.matrix(score)) {
      if(!all(genes %in% rownames(score))) {
        message(paste0(genes[which(!genes %in% rownames(score))], collapse = '-'), ' not in `score`')
        genes <- genes[which(genes %in% rownames(score))]
        if(length(genes) == 0) {
          return('No gene in `score`')
        }
      }
      score <- score[genes,]
    }
  }

  if(!is.vector(score)) {
    score <- Matrix::colMeans(score)
  }

  if(is.null(names(score))) {
    names(score) <- colnames(net)
  } else {
    score <- score[colnames(net)]
  }

  edges <- net@metadata$network$edges

  if(smooth) {
    network <- igraph::graph_from_adjacency_matrix(as.matrix(net@metadata$network$adj_matrix),
                                                   weighted = TRUE)
    new_score <- network_smoothing(network = network,
                                   scores = score[igraph::V(network)$name],
                                   n_nearest_smooth = n_nearest)

    new_score <- new_score[names(score)]
    new_score[which(is.na(new_score))] <- score[which(is.na(new_score))]
    score <- new_score

  }

  if(is.null(na.vec)) {
    na.vec <- rep(1, length(score))
  }

  alpha <- min_max_normalize(score)

  min_maxed_score <- alpha

  if(fix_alpha) {
    alpha <- NULL
  }

  if(is.null(main)) {
    if(length(genes) != nrow(net)) {
      main <- paste(genes[seq(1,min(length(genes), 5))], collapse = "-")
    }
  }

  weights <- apply(edges, 1, function(i) {
    i[["weight"]] <- as.numeric(i[["weight"]]) * min(c(as.numeric(min_maxed_score[i[["from"]]]),
                                                       as.numeric(min_maxed_score[i[["to"]]])))
    return(i[["weight"]])
  })
  edges$weight <- as.numeric(weights)

  layout <- data.frame('x' = net$X_coord, 'y' = net$Y_coord)
  rownames(layout) <- colnames(net)
  if(is.null(label_attr)) {
    label_attr <- rep(NA, nrow(layout))
  }
  layout$labels <- label_attr
  layout$fill <- score
  layout$color <- score
  if(is.null(alpha)) {
    alpha <- rep(1, nrow(layout))
  }
  layout$trans <- alpha
  size <- rep(max(max_size, 100/dim(layout)[1]), dim(layout)[1])
  layout$size <- size
  layout <- layout[order(score, decreasing = FALSE),]

  custom_theme <- ggplot2::theme(
    axis.title       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    legend.key       = ggplot2::element_blank(),
    legend.title     = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.grid       = ggplot2::element_blank(),
    plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
    plot.margin      = ggplot2::margin(1, 1, 1, 1, unit = "lines"),
    plot.title       = ggplot2::element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    )
  )

  p <- ggplot2::ggplot()

  if(!is.null(edges) & show_edges) {
    p <- p + ggplot2::geom_segment(data = edges, ggplot2::aes(x = from.x,
                                                              xend = to.x,
                                                              y = from.y,
                                                              yend = to.y),
                                   colour = 'grey',
                                   na.rm = TRUE,
                                   linewidth = edges$weight*(min(1, stroke+0.5)),
                                   alpha = edges$weight/max(edges$weight))
  }

  p <- p + ggplot2::geom_point(data = layout, mapping = ggplot2::aes(x = x,
                                                                     y = y,
                                                                     color = color,
                                                                     fill = fill,
                                                                     alpha = trans),
                               size = layout$size,
                               shape = 21,
                               stroke = stroke,
                               show.legend = FALSE)

  p <- p + ggplot2::scale_fill_gradientn(colors = palette,
                                         na.value = "#CCCCCC", guide = "colourbar", aesthetics = "fill") +
    ggplot2::scale_colour_gradientn(colors = colorspace::darken(palette, 0.1), na.value = "#CCCCCC",
                                    guide = NULL,
                                    aesthetics = "colour")

  p <- p + ggplot2::theme(legend.text = ggplot2::element_text(size = 10))

  p <- p + ggplot2::scale_alpha_identity() + custom_theme +
    ggplot2::labs(fill = "cat") +
    ggplot2::ggtitle(main) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = max(15, min(30, dim(layout)[1]/20)),
                                                      face = "bold"))

  if(all(is.na(label_attr))) {
    return(p)
  }

  all_cx <- c()
  all_cy <- c()
  all_colors <- c()
  for (g in unique(label_attr)) {
    c <- layout[which(layout$labels == g), ]
    if(dim(c)[1] == 1) {
      all_cx <- c(all_cx, as.numeric(c[1]))
      all_cy <- c(all_cy, as.numeric(c[2]))
      all_colors <- c(all_colors, as.vector(c[[4]]))
    } else {
      v_x <- as.numeric(c[,1])
      v_y <- as.numeric(c[,2])
      mx <- stats::median(v_x)
      my <- stats::median(as.numeric(c[,2]))
      external_point <- c(mx, my)
      points <- apply(c[,c(1, 2)], 2, as.numeric)
      dif_x <- (points[,1] - external_point[1])^2
      dif_y <- (points[,2] - external_point[2])^2
      distances <- sqrt(rowSums(cbind(dif_x, dif_y)))
      closest_point_index <- which.min(distances)
      closest_point <- c[closest_point_index, ]
      all_cx <- c(all_cx, as.numeric(closest_point[1]))
      all_cy <- c(all_cy, as.numeric(closest_point[2]))
      all_colors <- c(all_colors, as.vector(closest_point[[4]]))
    }
  }

  if(any(is.na(all_colors))) {
    all_colors[which(is.na(all_colors))] <- "grey90"
  }

  all_colors_hex <- scales::col_numeric(
    palette = palette,
    domain = range(layout$color, na.rm = TRUE)
  )(all_colors)

  p <- p + ggrepel::geom_label_repel(ggplot2::aes(x = all_cx,
                                                  y = all_cy,
                                                  label = unique(label_attr)),
                                     color = all_colors_hex,
                                     box.padding = 0.5,
                                     max.overlaps = Inf,
                                     show.legend = FALSE)
  return(p)

}
