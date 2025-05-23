#' Plot the rank of a gene (set) across the UMAP
#'
#' @inheritParams build_edges_df
#' @param umap_obj A UMAP object as given by the umap_signature_plot function or
#' by the add_to_reference function.
#' It has to be an object that contains the following chain:
#' `umap_obj$umap_obj$umap_out$layout`
#' @param genes A vector of one or multiple genes whose ranks have to be plotted
#' @param refExpression The expression matrix of the clusters in the UMAP
#' @param newExpression The expression matrix of the new clusters
#' @param label_attr The labels that are wanted in the new UMAP
#' @param col_vector A color vector for the UMAP. If NULL, different shades of
#' blues9 will be given (the darker the color, the higher the rank of the genes,
#' so the higher the signature)
#' @param new_clusters If the umap_obj contains also new clusters (if originated
#' from add_to_reference, 'New'), and this new clusters want to be
#' plotted bigger in the UMAP, their names have to be specified here
#' @param smooth A boolean variable to define if a smoothing step will be
#' performed to define the scores of each cluster
#' @param n_nearest_smooth The number of nearest neighbors to use to smooth the
#' rank value of expression. Defaults to all nodes.
#' @param title A title to give. If NULL, the names of the genes will be plotted
#' @param relative A boolean variable to define if the colors should be relative
#' or absolute. If smooth, the result will be relative regardless of the input
#' @param n_components The number of dimensions to use in the plot
#' @param show_edges A boolean variable to define if edges have to be shown or not
#' @param rank A boolean variable to define if the values are to be ranked
#' @param na.vec A vector with length nrow(refExpression) and with NA where you
#' want the plot to be transparent and without NAs where you want to see the
#' expression
#' @return The UMAP with colors and transparency proportional to the ranks of
#' the genes
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
#' new_M <- reference_signatures_correlation(new_clusterS, refS)
#' to_reference <- add_to_reference(annotated_M, new_M, rownames(annotated_M))
#' plotGenesRank(umap_obj = to_reference$New,
#' genes = 'Gene-8',
#' refExpression = S,
#' newExpression = new_clusterS,
#' new_clusters = rownames(new_M),
#' label_attr = colnames(S))
plotGenesRank <- function(umap_obj,
                          genes,
                          refExpression,
                          newExpression=NULL,
                          label_attr=NULL,
                          col_vector=NULL,
                          new_clusters=NULL,
                          title=NULL,
                          smooth=TRUE,
                          n_nearest_smooth=NULL,
                          relative=TRUE,
                          n_components=2,
                          weight_quantile=0.75,
                          weights_normalization_coef=NULL,
                          show_edges=TRUE,
                          rank=FALSE,
                          na.vec=NULL) {

  if(smooth & is.null(n_nearest_smooth)) {
    message('The number of nearest neighbors was not set and it will
            authomatically be the whole network. If the network is divided
            into groups, it is better to consider a number of nearest neighbors
            similar to the size of the groups')
  }

  if('umap_obj' %in% names(umap_obj)) {
    umap_obj <- umap_obj$umap_obj
  }

  if(smooth) {
    relative <- TRUE
  }

  if (!is.null(newExpression)) {
    common_genes <- intersect(rownames(refExpression), rownames(newExpression))
    signatures <- cbind(refExpression[common_genes, ], newExpression[common_genes, ])
  } else {
    signatures <- refExpression
  }

  signatures <- signatures[,rownames(umap_obj$umap_out$layout)]

  genes <- genes[which(genes %in% rownames(signatures))]

  if(length(genes) == 0) {
    return('Error: the given genes are not present in the signatures')
  }

  rowname <- rownames(signatures)

  if(rank) {
    rank_signatures <- apply(signatures, 2, function(i) {

      ranks <- seq(1, length(rowname))
      idxs <- order(i)
      new_rowname <- rowname[idxs]
      names(ranks) <- new_rowname

      ranks <- ranks[rowname]
      return(ranks)

    })
  } else {
    rank_signatures <- signatures
  }

  rank_value <- apply(rank_signatures, 2, function(i) {
    if(length(genes) == 1) {
      return(i[which(rowname == genes)])
    } else {
      return(mean(i[which(rowname %in% genes)]))
    }
  })

  if(is.null(umap_obj$refined_network)) {
    network_to_use <- umap_obj$umap_knn_igraph
    edges <- NULL
    distances <- umap_obj$umap_out$refined_network$distances
  } else {
    network_to_use <- umap_obj$refined_network
    edges <- umap_obj$umap_out$refined_network$edges
  }

  if (smooth) {
    rank_value <- network_smoothing(network = network_to_use,
                                    scores = rank_value[igraph::V(network_to_use)$name],
                                    n_nearest_smooth = n_nearest_smooth)
  }

  if(dim(umap_obj$umap_out$layout)[2] != n_components) {
    warning('The given layout has a different number of dimensions to those wanted by the user.
            The fewest dimensions will be considered, but taking the first dimensions of a layout
            with an higher number of dimensions not always is the right choice. Consider using a different
            layout, with the wanted dimensions.')
  }

  layout <- umap_obj$umap_out$layout[,seq(1,(min(n_components, dim(umap_obj$umap_out$layout)[2], 3)))]

  if(is.null(na.vec)) {
    na.vec <- rep(1, length(rank_value))
  }

  if(is.null(col_vector)) {
    if(relative) {
      breakpoints <- stats::quantile(seq((min(rank_value[which(!is.na(na.vec))])*rounder(1/min(rank_value[which(!is.na(na.vec))]))),(max(rank_value[which(!is.na(na.vec))])*rounder(1/min(rank_value[which(!is.na(na.vec))])))), probs = seq(0, 1, length.out = 6))
      breakpoints[which.min(breakpoints)] <- breakpoints[which.min(breakpoints)]-1
      breakpoints[which.max(breakpoints)] <- breakpoints[which.max(breakpoints)]+1

      col_vector <- as.vector(cut(rank_value*rounder(1/min(rank_value)), breaks = breakpoints, labels = grDevices::blues9[c(3,5,6,7,9)]))
    } else {
      breakpoints <- stats::quantile(seq(1,max(rank_signatures)), probs = seq(0, 1, length.out = 6))
      breakpoints[which.min(breakpoints)] <- breakpoints[which.min(breakpoints)]-1
      breakpoints[which.max(breakpoints)] <- breakpoints[which.max(breakpoints)]+1

      col_vector <- as.vector(cut(rank_value, breaks = breakpoints, labels = grDevices::blues9[c(3,5,6,7,9)]))
    }
  }

  if(any(is.na(na.vec))) {
    col_vector[which(is.na(na.vec))] <- NA
  }

  names(col_vector) <- rank_value
  col_vector <- col_vector[order(rank_value)]

  if(is.null(new_clusters)) {
    if(n_components == 2) {
      size <- rep(max(2, 100/dim(layout)[1]), dim(layout)[1])
    } else {
      size <- rep(max(200, 10000/dim(layout)[1]), dim(layout)[1])
    }
    color_attr <- label_attr
  } else {
    if(n_components == 2) {
      size <- rep(max(2, 100/dim(layout)[1]), dim(layout)[1])
      size[which(rownames(layout) %in% new_clusters)] <- size[which(rownames(layout) %in% new_clusters)] * 2
    } else {
      size <- rep(max(200, 10000/dim(layout)[1]), dim(layout)[1])
      size[which(rownames(layout) %in% new_clusters)] <- size[which(rownames(layout) %in% new_clusters)] * 6
    }
    if((length(label_attr) != dim(layout)[1]) && (!is.null(label_attr))) {
      label_attr <- c(label_attr, new_clusters)
    }
    color_attr <- label_attr
  }

  layout <- cbind(layout, 'Size' = size)
  layout <- cbind(layout, 'Colors' = col_vector[as.character(rank_value)])

  layout <- cbind(layout, 'Group' = col_vector[as.character(rank_value)])

  if(!is.null(label_attr)) {
    layout <- cbind(layout, 'SubGroup' = color_attr)
  }

  alpha <- (rank_value - min(rank_value))/(max(rank_value) - min(rank_value))

  if(is.null(title)) {
    if(length(genes) > 1) {
      title <- paste(genes, collapse = "-")
    } else {
      title <- genes
    }
  }

  min_maxed_rank_value <- (rank_value - min(rank_value))/(max(rank_value)-min(rank_value))

  if(any(is.na(col_vector))) {
    alpha[order(rank_value)][which(is.na(col_vector))] <- 0
    min_maxed_rank_value[order(rank_value)][which(is.na(col_vector))] <- 0
  }

  weights <- apply(edges, 1, function(i) {
    i[['weight']] <- as.numeric(i[['weight']]) * min(c(as.numeric(min_maxed_rank_value[i[['from']]]), as.numeric(min_maxed_rank_value[i[['to']]])))
    return(i[['weight']])
  })

  edges$weight <- as.numeric(weights)

  if(n_components == 2) {
    if(!is.null(label_attr)) {
      p <- umap_pointsize(layout = layout,
                          color_attr = layout[,'Colors'],
                          label_attr = layout[,'SubGroup'],
                          title = title,
                          legend = FALSE,
                          only_new_points= FALSE,
                          alpha = alpha,
                          edges = edges,
                          show_edges = show_edges)
    } else {
      p <- umap_pointsize(layout = layout,
                          color_attr = layout[,'Colors'],
                          label_attr = label_attr,
                          title = title,
                          legend =FALSE,
                          only_new_points =FALSE,
                          alpha = alpha,
                          edges = edges,
                          show_edges = show_edges)
    }
  } else {
    p <- plotly::plot_ly(x = layout[,1],
                         y = layout[,2],
                         z = layout[,3],
                         type="scatter3d",
                         mode="markers",
                         colors = col_vector,
                         color = rank_value,
                         size = I(as.numeric(layout[,'Size'])))

    p$x$layoutAttrs <- list('Params' = list('title' = title,
                                            'scene' = list('xaxis' = list('title' = 'UMAP-1'),
                                                           'yaxis' = list('title' = 'UMAP-2'),
                                                           'zaxis' = list('title' = 'UMAP-3'))))
  }

  return(p)
}
