#' Get genes driving the newly mapped points
#'
#' @param new_umap_obj The umap object coming from add_to_annotated_reference,
#' under `New$umap_obj`
#' @param newSignatures The signatures of the new points
#' @param n_genes The number of top genes for each new point to be returned
#' @param gene_set A starting gene set from which to extract the top genes
#' @param threshold The signatures level threshold under which genes are not
#' considered. It's specific for signatures of each new point.
#'
#' @return A list of genes for each new point mapped
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
#' mapped <- add_to_annotated_reference(annotated_M,
#' new_M,
#' annotated_M$`Best.Assignment`)
#' getDrivingGenes(mapped$New$umap_obj, new_clusterS)
getDrivingGenes <- function(new_umap_obj,
                            newSignatures,
                            n_genes = 100,
                            threshold=0.90,
                            gene_set = NULL) {

  if('umap_obj' %in% names(new_umap_obj)) {
    new_umap_obj <- new_umap_obj$umap_obj
  }

  newSignatures <- newSignatures[which(rowSums(abs(newSignatures)) != 0),]

  if(is.null(gene_set)) {
    gene_set <- rownames(newSignatures)
  } else {
    gene_set <- intersect(gene_set, rownames(newSignatures))
  }

  new_layout <- new_umap_obj$umap_out$layout[colnames(newSignatures),c(1,2)]

  n_genes <- min(n_genes, length(gene_set))

  newSignatures <- newSignatures[,rownames(new_layout)]

  x_coord <- c()
  y_coord <- c()
  color <- c()
  group <- c()
  gene_list <- list()

  x_values <- c()
  y_values <- c()
  for(gene in gene_set) {
    weights <- newSignatures[gene,] + abs(min(newSignatures[gene,]))
    weights <- weights - mean(weights)
    weights[which(weights < 0)] <- 0

    v_x <- stats::weighted.mean(x = as.numeric(new_layout[,1]),
                                w = weights)
    if(v_x > max(as.numeric(new_layout[,1]))) {
      v_x <- max(as.numeric(new_layout[,1]))
    }
    if(v_x < min(as.numeric(new_layout[,1]))) {
      v_x <- min(as.numeric(new_layout[,1]))
    }
    x_values <- c(x_values, v_x)
    v_y <- stats::weighted.mean(as.numeric(new_layout[,2]),
                                w = weights)
    if(v_y > max(as.numeric(new_layout[,2]))) {
      v_y <- max(as.numeric(new_layout[,2]))
    }
    if(v_y < min(as.numeric(new_layout[,2]))) {
      v_y <- min(as.numeric(new_layout[,2]))
    }
    y_values <- c(y_values, v_y)
  }

  gene_points <- cbind(x_values, y_values)
  rownames(gene_points) <- gene_set
  n_ref <- cbind(as.numeric(new_layout[,1]), as.numeric(new_layout[,2]))

  x_values_original <- c()
  y_values_original <- c()
  min_distances <- c()
  idxs_dist <- c()
  for(gene in seq(1, length(gene_set))) {
    gp <- gene_points[gene,]
    dif_x <- (n_ref[,1] - gp[1]) ^ 2
    dif_y <- (n_ref[,2] - gp[2]) ^ 2
    distances <- sqrt(rowSums(cbind(dif_x, dif_y)))
    idx_dist <- which.min(distances)
    min_distances <- c(min_distances, distances[idx_dist])
    idxs_dist <- c(idxs_dist, idx_dist)
    x_values_original <- c(x_values_original, n_ref[idx_dist, 1])
    y_values_original <- c(y_values_original, n_ref[idx_dist, 2])
  }

  x_coord <- c(x_coord, x_values_original)
  y_coord <- c(y_coord, y_values_original)

  names(x_coord) <- gene_set
  names(y_coord) <- gene_set

  closest_groups <- colnames(newSignatures)[idxs_dist]
  names(closest_groups) <- gene_set
  names(min_distances) <- gene_set

  for(g in colnames(newSignatures)) {

    closest_genes <- names(closest_groups)[which(closest_groups == g)]
    distances <- min_distances[closest_genes]
    closest_genes <- closest_genes[order(distances)]
    n_new_layout <- new_layout[which(rownames(new_layout) == g),]
    filt_signatures <- newSignatures[,rownames(new_layout)[which(rownames(new_layout) == g)]]

    if (!is.vector(filt_signatures)) {
      if (dim(filt_signatures)[2] == 1) {
        filt_signatures <- filt_signatures[,1]
      }
      to_remove <- which(apply(filt_signatures, 1, max) == 0)
      if(length(to_remove) > 0) {
        filt_signatures <- filt_signatures[-which(apply(filt_signatures,
                                                        1, max) == 0), ]
      }
    } else {
      to_remove <- which(filt_signatures == 0)
      if(length(to_remove) > 0) {
        filt_signatures <- filt_signatures[-to_remove]
      }
    }

    if(is.vector(filt_signatures)) {
      sorted <- sort(filt_signatures)
      names_sorted <- names(sorted)
      ranks <- seq(1, length(sorted))
      names(ranks) <- names(sorted)
      ranks <- ranks[names(filt_signatures)]
      whole <- S4Vectors::SimpleList('Selected_genes' = names(ranks), 'All_ranks' = ranks)
    } else {
      whole <- getSignatureMeans(filt_signatures, n_genes = nrow(filt_signatures))
    }
    to_keep_genes <- whole$Selected_genes
    scaled_ranks <- scale(whole$All_ranks)
    threshold_ranks <- stats::quantile(scaled_ranks, threshold)
    to_keep_genes <- intersect(to_keep_genes, rownames(scaled_ranks)[which(scaled_ranks > threshold_ranks)])
    to_keep_genes <- to_keep_genes[which(to_keep_genes %in% gene_set)]

    if(is.null(to_keep_genes) | length(to_keep_genes) == 0) {
      next()
    }

    to_keep_genes <- intersect(closest_genes, to_keep_genes)

    to_keep_genes <- to_keep_genes[seq(1, min(length(to_keep_genes), n_genes))]

    if(length(to_keep_genes)==1) {
      if(is.na(to_keep_genes)) {
        next()}
    }

    gene_list[[g]] <- to_keep_genes
  }

  return(gene_list)

}
