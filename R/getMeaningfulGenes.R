#' Identifies the genes responsible for the mapping of new clusters onto the UMAP
#'
#' @param refSignatures The signatures of the clusters in the UMAP (not new ones)
#' @param newSignatures The signatures of the new clusters
#' @param umap_obj A UMAP object as given by the add_to_annotated_reference
#' function (but specifying 'New').
#' It has to be an object that contains the following chain:
#' `umap_obj$umap_obj$umap_out$layout` or `umap_obj$umap_out$layout`
#' @param new_clusters The names of the new clusters
#' @param n_genes The number of genes to be returned
#' @param threshold A threshold for the quality of the genes (hence, genes might
#' be less than n_genes in the end)
#' @param n_nearest Number of nearest neighbors to consider. It defaults to the
#' number of nearest neighbors used to construct the UMAP
#'
#' @return A list of lists, in which each element of the first list considers one
#' new cluster and it gives a list with the the impactful genes and the
#' confidence
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(100,0.1,7), ncol = 5))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' annotated_M <- reference_signatures_correlation(S, refS)
#' new_clusterS <- FC_signatures(matrix(runif(80,0,10), ncol = 4))
#' rownames(new_clusterS) <- paste0('Gene-', seq(1, dim(new_clusterS)[1]))
#' colnames(new_clusterS) <- paste0('New-', seq(1, dim(new_clusterS)[2]))
#' new_M <- reference_signatures_correlation(new_clusterS, refS)
#' to_reference <- add_to_annotated_reference(annotated_M,
#' new_M,
#' annotated_M$`Best.Assignment`)
#' getMeaningfulGenes(umap_obj = to_reference$New,
#' new_clusters = rownames(new_M),
#' refSignatures = S,
#' newSignatures = new_clusterS,
#' threshold = 0.5)
getMeaningfulGenes <- function(refSignatures,
                               newSignatures,
                               umap_obj,
                               new_clusters,
                               n_genes=20,
                               threshold=0.99,
                               n_nearest=NULL) {

  if('umap_obj' %in% names(umap_obj)) {
    umap_obj <- umap_obj$umap_obj
  }

  if(is.null(umap_obj$umap_out$refined_network)) {
    distances <- umap_obj$umap_out$knn$distances
    indexes <- umap_obj$umap_out$knn$indexes
  } else {
    distances <- umap_obj$umap_out$refined_network$distances
    indexes <- umap_obj$umap_out$refined_network$indexes
  }

  if(is.null(n_nearest)) {

      n_nearest <- dim(distances)[2]

  } else {

      n_nearest <- min(n_nearest, dim(distances)[2])

  }

  distances <- distances[,seq(1,n_nearest)]
  indexes <- indexes[,seq(1,n_nearest)]

  common_genes <- intersect(rownames(refSignatures), rownames(newSignatures))
  refSignatures <- refSignatures[common_genes,]
  newSignatures <- newSignatures[common_genes,]

  n_genes <- min(n_genes, nrow(refSignatures))

  impactful_genes <- list()

  for(c in new_clusters) {

    nn <- rownames(indexes)[indexes[c,]]
    dist <- distances[c,]
    idx <- which(!nn %in% new_clusters)
    nn <- nn[idx]
    dist <- dist[idx]

    f_signatures <- refSignatures[,nn]

    colnames(f_signatures) <- nn

    weights <- 1/dist

    mean_sig <- getSignatureMeans(signatures = f_signatures,
                                  n_genes = n_genes,
                                  weights = weights)

    scaled_ranks <- scale(mean_sig$All_ranks)
    threshold_ranks <- stats::quantile(scaled_ranks, threshold)

    sel_genes <- mean_sig$Selected_genes

    sel_genes <- intersect(sel_genes, rownames(scaled_ranks)[which(scaled_ranks > threshold_ranks)])

    confidence <- mean_sig$All_ranks/max(mean_sig$All_ranks)
    confidence <- confidence[sel_genes]

    impactful_genes[[c]] <- S4Vectors::SimpleList('Genes' = sel_genes,
                                                  'Confidence' = confidence)

  }

  return(impactful_genes)

}
