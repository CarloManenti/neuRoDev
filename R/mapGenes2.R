#' Creates a gene-view over a cluster-view UMAP, in which the most impactful genes
#' for the layout of the cluster-view UMAP are displayed
#'
#' @param umap_obj A UMAP object as given by the add_to_reference
#' function (but specifying 'New').
#' It has to be an object that contains the following chain:
#' `umap_obj$umap_out$layout`
#' @param new_clusters The names of the new clusters
#' @param refSignatures The signatures of the clusters of the UMAP (excluding
#' new clusters)
#' @param col_vector The color palette to be used in the UMAP
#' @param color_attr The annotation labels to be used in the UMAP
#' @param n_genes The number of genes that have to be displayed
#' @param new_points_col The color of the new clusters in the UMAP
#' @param gene_set A vector from which to choose the displayed genes (if NULL,
#' all genes are considered)
#' @param threshold A threshold to keep only the most significant genes, even
#' if it might lead to less than n_genes to be displayed.
#' @param proportional A boolean variable that defines if the genes to be
#' displayed are proportional in number to the number of clusters in each color_attr
#' defined in color_attr. If FALSE, each color_attr in color_attr will
#' have approximately the same number of genes displayed
#' @param title The title of the plot
#'
#' @return A SimpleList with the original cluster-view and the new gene-view,
#' along with a list of the displayed genes for each color_attr
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(stats::runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(stats::runif(100,0.1,7), ncol = 5))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' annotated_M <- reference_signatures_correlation(S, refS)
#' new_clusterS <- FC_signatures(matrix(stats::runif(80,0,10), ncol = 4))
#' rownames(new_clusterS) <- paste0('Gene-', seq(1, dim(new_clusterS)[1]))
#' colnames(new_clusterS) <- paste0('New-', seq(1, dim(new_clusterS)[2]))
#' new_M <- reference_signatures_correlation(new_clusterS, refS)
#' to_reference <- add_to_reference(annotated_M,
#' new_M,
#' annotated_M$`Best.Assignment`)
#' mapGenes2(umap_obj = to_reference$New,
#' new_clusters = rownames(new_M),
#' refSignatures = S,
#' color_attr = annotated_M$`Best.Assignment`,
#' threshold = 0.05)
mapGenes2 <- function(umap_obj,
                     refSignatures,
                     color_attr,
                     new_clusters=NULL,
                     col_vector=NULL,
                     n_genes=100,
                     new_points_col="#FF0000",
                     gene_set=NULL,
                     threshold=0.95,
                     proportional=TRUE,
                     title=NULL) {

  if('umap_obj' %in% names(umap_obj)) {
    umap_obj <- umap_obj$umap_obj
  }

  if(is.null(title)) {
    title <- 'Gene View'
  }

  if(is.null(gene_set)) {
    gene_set <- rownames(refSignatures)
  } else {
    gene_set <- intersect(gene_set, rownames(refSignatures))
    n_genes <- min(n_genes, length(gene_set))
  }

  n_genes <- min(n_genes, length(gene_set))

  layout <- umap_obj$umap_out$layout[,c(1,2)]

  if(is.null(col_vector)) {
    col_vector <- Polychrome::createPalette(length(unique(color_attr))+2,
                                            c("#ffffff", "#ff0000", "#00ff00", "#0000ff", "#000000"))
    col_vector <- col_vector[seq(3,length(col_vector))]
    col_vector <- col_vector[seq(1,length(unique(color_attr)))]
    names(col_vector) <- unique(color_attr)
  }

  if(!is.null(new_clusters)) {
    ref_layout <- layout[which(!rownames(layout) %in% new_clusters),]
    new_layout <- layout[which(rownames(layout) %in% new_clusters),]
  } else {
    ref_layout <- layout
  }

  refSignatures <- refSignatures[,rownames(ref_layout)]

  ref_layout <- cbind(ref_layout, 'Group' = color_attr)
  ref_layout <- cbind(ref_layout, 'Color' = col_vector[color_attr])

  all_genes <- c()
  color <- c()
  group <- c()
  gene_list <- list()

  x_values <- c()
  y_values <- c()
  for(gene in gene_set) {
    weights <- refSignatures[gene,] + abs(min(refSignatures[gene,]))
    weights <- weights - mean(weights)
    weights[which(weights < 0)] <- 0

    v_x <- stats::weighted.mean(x = as.numeric(ref_layout[,1]),
                         w = weights)
    if(v_x > max(as.numeric(ref_layout[,1]))) {
      v_x <- max(as.numeric(ref_layout[,1]))
    }
    if(v_x < min(as.numeric(ref_layout[,1]))) {
      v_x <- min(as.numeric(ref_layout[,1]))
    }
    x_values <- c(x_values, v_x)
    v_y <- stats::weighted.mean(as.numeric(ref_layout[,2]),
                         w = weights)
    if(v_y > max(as.numeric(ref_layout[,2]))) {
      v_y <- max(as.numeric(ref_layout[,2]))
    }
    if(v_y < min(as.numeric(ref_layout[,2]))) {
      v_y <- min(as.numeric(ref_layout[,2]))
    }
    y_values <- c(y_values, v_y)
  }

  gene_points <- cbind(x_values, y_values)
  n_ref <- cbind(as.numeric(ref_layout[,1]), as.numeric(ref_layout[,2]))

  x_coord <- c()
  y_coord <- c()
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
    x_coord <- c(x_coord, n_ref[idx_dist, 1])
    y_coord <- c(y_coord, n_ref[idx_dist, 2])
  }

  names(idxs_dist) <- gene_set
  names(x_coord) <- gene_set
  names(y_coord) <- gene_set

  closest_groups <- color_attr[idxs_dist]
  names(closest_groups) <- gene_set
  names(min_distances) <- gene_set

  for(g in unique(color_attr)) {

    if(proportional) {
      n_genes_g <- max(2,round(n_genes*(table(color_attr) / sum(table(color_attr)))[g]))
    } else {
      n_genes_g <- max(2, round(n_genes/length(unique(color_attr))))
    }

    closest_genes <- names(closest_groups)[which(closest_groups == g)]

    if(length(closest_genes) == 0) {
      priority <- FALSE
      sub_n_ref <- n_ref[which(color_attr == g),]
      sub_x_coord <- c()
      sub_y_coord <- c()
      sub_min_distances <- c()
      sub_idxs_dist <- c()

      for(gene in seq(1, length(gene_set))) {
        gp <- gene_points[gene,]
        dif_x <- (sub_n_ref[,1] - gp[1]) ^ 2
        dif_y <- (sub_n_ref[,2] - gp[2]) ^ 2
        distances <- sqrt(rowSums(cbind(dif_x, dif_y)))
        idx_dist <- which.min(distances)
        sub_min_distances <- c(sub_min_distances, distances[idx_dist])
        sub_idxs_dist <- c(sub_idxs_dist, idx_dist)
      }

      names(sub_min_distances) <- gene_set
      sub_min_distances <- sort(sub_min_distances, decreasing = FALSE)

      closest_genes <- names(sub_min_distances)[seq(1, n_genes_g * 10)]

      min_distances[closest_genes] <- sub_min_distances[seq(1, n_genes_g * 10)]

    } else {
      priority <- TRUE
    }

    distances <- min_distances[closest_genes]
    closest_genes <- closest_genes[order(distances)]
    n_ref_layout <- ref_layout[which(ref_layout[,'Group'] == g),]
    filt_signatures <- refSignatures[,rownames(ref_layout)[which(ref_layout[,'Group'] == g)]]

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
    to_keep_genes <- to_keep_genes[which(to_keep_genes %in% gene_set[which(!gene_set %in% all_genes)])]

    if(is.null(to_keep_genes) | length(to_keep_genes) == 0) {
      next()
    }

    to_keep_genes <- intersect(closest_genes, to_keep_genes)

    to_keep_genes <- to_keep_genes[seq(1, min(length(to_keep_genes), n_genes_g))]

    if(length(to_keep_genes)==1) {
      if(is.na(to_keep_genes)) {
        next()}
    }

    if(priority) {
      x_coord[to_keep_genes] <- n_ref[idxs_dist[to_keep_genes], 1]
      y_coord[to_keep_genes] <- n_ref[idxs_dist[to_keep_genes], 2]
    }

    gene_list[[g]] <- to_keep_genes
    all_genes <- c(all_genes, to_keep_genes)

    if(is.vector(n_ref_layout)) {
      color <- c(color, rep(unique(n_ref_layout['Color']), length(to_keep_genes)))
      group <- c(group, rep(unique(n_ref_layout['Group']), length(to_keep_genes)))
    } else {
      color <- c(color, rep(unique(n_ref_layout[,'Color']), length(to_keep_genes)))
      group <- c(group, rep(unique(n_ref_layout[,'Group']), length(to_keep_genes)))
    }
  }

  genes_layout <- cbind(x_coord[all_genes], y_coord[all_genes])

  duplicated_indices <- duplicated(genes_layout) | duplicated(genes_layout, fromLast = TRUE)

  while (any(duplicated_indices)) {
    genes_layout[duplicated_indices, ] <- genes_layout[duplicated_indices, ] + stats::runif(sum(duplicated_indices), -0.0001, 0.0001)

    duplicated_indices <- duplicated(genes_layout) | duplicated(genes_layout, fromLast = TRUE)
  }

  rownames(genes_layout) <- all_genes

  genes_layout <- cbind(genes_layout, 'Size' = rep(max(2,100/dim(genes_layout)[1]), dim(genes_layout)[1]))
  genes_layout <- cbind(genes_layout, 'Color' = color)
  genes_layout <- cbind(genes_layout, 'Group' = group)

  if(!is.null(new_clusters)) {
    new_points_layout <- cbind(new_layout, 'Size' = rep(max(2,100/dim(genes_layout)[1])*2, dim(new_layout)[1]))
    if(length(new_points_col)==1) {
      new_points_layout <- cbind(new_points_layout, 'Color' = rep(new_points_col, dim(new_layout)[1]))
      new_points_layout <- cbind(new_points_layout, 'Group' = rep('New', dim(new_layout)[1]))
    } else {
      new_points_layout <- cbind(new_points_layout, 'Color' = new_points_col)
      new_points_layout <- cbind(new_points_layout, 'Group' = new_clusters)
    }

    final_layout <- rbind(genes_layout, new_points_layout)

    final_layout <- cbind(final_layout, 'SubGroup' = c(all_genes, new_clusters))

    p <- umap_pointsize(layout = final_layout, color_attr = final_layout[,'Group'], label_attr = final_layout[,'SubGroup'], new_points_col = NULL, legend=TRUE, title = title)
  } else {

    final_layout <- cbind(genes_layout, 'SubGroup' = all_genes)

    p <- umap_pointsize(layout = final_layout, color_attr = final_layout[,'Group'], label_attr = final_layout[,'SubGroup'], new_points_col = NULL, legend=TRUE, title = title)

  }

  return(S4Vectors::SimpleList('Cluster_View' = umap_obj$umap_plot, 'Gene_View' = p, 'Genes' = gene_list))

}
