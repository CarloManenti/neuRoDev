#' Prediction of previous clusters given new ones
#'
#'  Given a UMAP object umap_obj, obtained from the function
#'  add_to_reference (the New object), and names for new_clusters,
#'  it computes the diffusion of the new_clusters labels in the UMAP network
#'  defined in umap_obj. If a col_vector is given, the plot will be coloured
#'  as such. If new_points_col is given, the points corresponding to
#'  new_clusters will be plotted twice as big and coloured as such.
#'  If legend=TRUE, a legend will be displayed instead of the labels.
#' @inheritParams umap_plot_same_layout
#' @param update_reference A boolean variable. If TRUE, the reference will be
#' updated based on the label with the highest score
#' @param reference The reference to update if update_reference = TRUE
#' @param layout A possible starting layout as returned by
#' umap_graph_clustering (defaults to use the one given in `umap_obj`)
#' @param only_best_association A boolean variable to decide if only the best
#' association is returned
#'
#' @return A list with the distance scores of the new, labelled clusters and the
#' old clusters, plus the best association and its confidence and the UMAP plot.
#' If update_reference=TRUE, also the updated reference
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
#' res <- add_to_reference(annotated_M,
#' new_M,
#' annotated_M$`Best.Assignment`)
#' umap_obj <- res
#' new_clusters <- rownames(new_M)
#' labelPrediction(umap_obj, new_clusters)
labelPrediction <- function(umap_obj,
                            new_clusters,
                            col_vector=NULL,
                            new_points_col=NULL,
                            legend=FALSE,
                            layout=NULL,
                            update_reference=FALSE,
                            reference=NULL,
                            weight_quantile=0.75,
                            weights_normalization_coef=NULL,
                            show_edges=TRUE,
                            only_best_association=FALSE) {

  old_umap_obj <- umap_obj$Original$umap_obj

  umap_obj <- umap_obj$New$umap_obj

  net <- umap_obj$refined_network
  edges <- old_umap_obj$umap_out$refined_network$edges

  d <- umap_obj$umap_out$refined_network$distances
  ind <- umap_obj$umap_out$refined_network$indexes

  if(is.null(net)) {
    net <- umap_obj$umap_knn_igraph
    edges <- NULL
    d <- umap_obj$umap_out$knn$distances
    ind <- umap_obj$umap_out$knn$indexes
  } else {
    igraph::E(net)$weights <- 1-igraph::E(net)$weight
  }

  full_c_score <- matrix(0, nrow = nrow(d), ncol = ncol(d))
  rownames(full_c_score) <- rownames(d)
  colnames(full_c_score) <- rownames(d)

  full_c_score[cbind(rep(seq(1,nrow(d)), ncol(d)),
                as.vector(ind))] <- as.vector(d)

  diag(full_c_score) <- 0

  c_score <- full_c_score[which(!rownames(full_c_score) %in% new_clusters),new_clusters]

  best_association <- unlist(lapply(seq_len(dim(c_score)[1]), function(i) {
    colnames(c_score)[which.min(c_score[i,])]
  }))

  names(best_association) <- rownames(c_score)
  best_association[which(rowMeans(c_score) == 1)] <- 'NA'

  v <- new_clusters
  names(v) <- new_clusters

  best_association <- c(best_association, v)

  if(only_best_association) {
    return(best_association)
  }

  if(is.null(col_vector)) {
    col_vector <- Polychrome::createPalette(length(new_clusters)+2+length(new_points_col),
                                c("#ffffff", "#ff0000", "#00ff00", "#0000ff", "#000000"))

    col_vector <- col_vector[seq(3,length(col_vector))]
    col_vector <- col_vector[which(!col_vector %in% c("#ffffff", "#ff0000", new_points_col))]
    col_vector <- col_vector[seq(1,length(new_clusters))]
    names(col_vector) <- new_clusters
  }

  if(!is.null(new_points_col)) {
    best_association[which(!names(best_association) %in% new_clusters)] <- paste0('Predicted-', best_association[which(!names(best_association) %in% new_clusters)])

    names(col_vector) <- paste0('Predicted-', names(col_vector))

    if(is.vector(new_points_col)) {
      new_points_col <- rep(new_points_col, length(new_clusters))
      names(new_points_col) <- new_clusters
    }
    if(is.null(names(new_points_col))) {
      names(new_points_col) <- new_clusters
    }
    col_vector <- c(col_vector, new_points_col)
  }

  if(is.null(layout)) {
    layout <- umap_obj$umap_out$layout
  }

  best_association <- best_association[rownames(layout)]

  while(length(which(best_association == 'NA')) != 0) {

    sub_c_score <- full_c_score[names(best_association)[which(best_association == 'NA')],
                                names(best_association)[which(best_association != 'NA')]]

    if(length(which(best_association == 'NA')) == 1) {
      ordered_na_association <- names(best_association)[which(best_association == 'NA')]
      best_association[ordered_na_association[1]] <- best_association[names(sub_c_score)[which.min(sub_c_score)]]
    } else {
      ordered_na_association <- rownames(sub_c_score)[order(apply(sub_c_score, 1, min))]
      best_association[ordered_na_association[1]] <- best_association[colnames(sub_c_score)[which.min(sub_c_score[ordered_na_association[1],])]]
    }

  }

  layout <- cbind(layout, 'Size' = max(2,100/dim(layout)[1]))

  if(!is.null(new_points_col)) {
    layout[which(rownames(layout) %in% new_clusters), 'Size'] <- max(2,100/dim(layout)[1]) * 2
  }

  layout <- cbind(layout, 'Colors' = col_vector[(best_association)])

  layout <- cbind(layout, 'Group' = (best_association))

  layout <- cbind(layout, 'SubGroup' = (best_association))

  layout <- layout[order(match(layout[,'Group'], (col_vector))),]

  intra_confidence <- unlist(lapply(seq(1,dim(c_score)[1]), function(i) {
    1/((c_score)[i, which.min(c_score[i,])])
  }))
  names(intra_confidence) <- rownames(c_score)
  intra_confidence <- intra_confidence/max(intra_confidence)

  inter_confidence <- unlist(lapply(seq(1,dim(c_score)[1]), function(i) {
    ref_col <- which.min(c_score[i,])
    sum_v <- 0
    for(n in seq_len(ncol(c_score))) {
      sum_v <- sum_v + abs(c_score[i, ref_col]-c_score[i, n])
    }
    return(sum_v)
  }))
  names(inter_confidence) <- rownames(c_score)
  inter_confidence <- inter_confidence/max(inter_confidence)

  confidence <- (intra_confidence + inter_confidence)/2

  layout <- layout[which(!rownames(layout) %in% new_clusters),]

  if(legend) {
    p <- umap_pointsize(layout,
                        color_attr = layout[,'Group'],
                        legend = legend,
                        label_attr = NULL,
                        new_points_col = new_points_col,
                        title = 'Prediction',
                        only_new_points = FALSE,
                        edges = edges,
                        show_edges = show_edges)
  } else {
    p <- umap_pointsize(layout,
                        color_attr = layout[,'Group'],
                        legend = legend,
                        label_attr = layout[,'SubGroup'],
                        new_points_col = new_points_col,
                        title = 'Prediction',
                        only_new_points = FALSE,
                        edges = edges,
                        show_edges = show_edges)
  }

  if(update_reference) {
    if(is.null(reference)) {
      stop('No reference provided.')
    } else {
      labeled_reference <- reference
      labeled_reference$`Predicted.Best.Assignment` <- best_association[rownames(labeled_reference)]
      labeled_reference$Confidence <- confidence
      return(S4Vectors::SimpleList('Distance_Scores' = c_score,
                                   'Best_Association' = best_association,
                                   'Confidence' = confidence,
                                   'Plot' = p,
                                   'Reference' = labeled_reference))
    }
  }

  return(S4Vectors::SimpleList('Distance_Scores' = c_score,
                               'Best_Association' = best_association,
                               'Confidence' = confidence,
                               'Plot' = p))

}
