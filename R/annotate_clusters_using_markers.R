#' annotate_clusters_using_markers
#'
#' @param umap_obj A UMAP object as given by the umap_signature_plot function or
#' by the add_to_reference function
#' (but specifying if 'New' or 'Original').
#' It has to be an object that contains the following chain:
#' `umap_obj$umap_out$layout`
#' @param marker_list A named list that contains marker genes for given cell
#' types
#' @param signatures A matrix of the signatures of the clusters in the UMAP
#'
#' @return An annotation per cluster
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
#' umap_object <- umap_signature_plot(annotated_M)
#' annotate_clusters_using_markers(umap_obj = umap_object,
#' marker_list = list('A' = c('Gene-1', 'Gene-2'),
#' 'B' = c('Gene-1', 'Gene-3', 'Gene-8'),
#' 'C' = c('Gene-6', 'Gene-3')),
#' signatures = S)
annotate_clusters_using_markers <- function(umap_obj,
                                            marker_list,
                                            signatures) {


  if('umap_obj' %in% names(umap_obj)) {
    umap_obj <- umap_obj$umap_obj
  }

  signatures <- signatures[,rownames(umap_obj$umap_out$layout)]

  marker_list <- lapply(marker_list, function(i) {
    i[which(i %in% rownames(signatures))]
  })

  rowname <- rownames(signatures)

  rank_signatures <- apply(signatures, 2, function(i) {

    ranks <- seq(1, length(rowname))
    idxs <- order(i)
    new_rowname <- rowname[idxs]
    names(ranks) <- new_rowname

    ranks <- ranks[rowname]
    return(ranks)

  })

  if(is.null(umap_obj$refined_network)) {
    network_to_use <- umap_obj$umap_knn_igraph
  } else {
    network_to_use <- umap_obj$refined_network
  }

  rank_values <- lapply(marker_list, function(genes) {

    rank_value <- apply(rank_signatures, 2, function(i) {
      if(length(genes) == 1) {
        return(i[which(rowname == genes)])
      } else {
        return(mean(i[which(rowname %in% genes)]))
      }
    })

    rank_value_dif <- igraph::page_rank(network_to_use,
                                        personalized = rank_value[igraph::V(network_to_use)$name])
    rank_value_dif <- rank_value_dif$vector
    rank_value <- rank_value_dif[names(rank_value)]

    return(rank_value)

  })

  annotations <- unlist(lapply(seq(1, length(rank_values[[1]])), function(i) {

    all_values <- unlist(lapply(rank_values, function(x) {x[i]}))

    annotation <- names(rank_values)[which.max(all_values)]

    return(annotation)

  }))

  names(annotations) <- colnames(signatures)

  return(annotations)

}
