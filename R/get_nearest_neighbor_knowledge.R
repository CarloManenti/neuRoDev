#' Get nearest neighbor knowledge
#'
#' Given an annotated reference dataframe reference_df,
#' a set of new clusters to analyze in new_clusters, a umap object as run
#' with add_to_reference, which contains the new clusters to
#' analyze, and a given annotation label, it returns information about the
#' nearest neighbors for each cluster and as a whole for that given
#' annotation label. If remove_low_confidence=TRUE, it doesn't consider in the
#' analysis low confidence neighbors, as defined with the distance_threshold
#' parameter. You can provide a col_vector palette for the final barplot or
#' boxplot.
#' @inheritParams nearest_neighbor_annotation
#' @param ylim The Y limits of the boxplot
#' @param color_attr The labels to use for the annotation
#'
#' @return A list of the whole distribution of the variable in the nearest
#' neighbors, the whole frequency and whole percentages, then also divided by
#' cluster, plus a Heatmap and a barplot if the label is a character, a boxplot
#' for the whole distribution and one per new cluster if it is numeric
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
#' umap_obj <- res$New
#' new_clusters = rownames(new_M)
#' get_nearest_neighbor_knowledge(annotated_M,
#' new_clusters,
#' umap_obj,
#' color_attr = 'Best.Assignment')
get_nearest_neighbor_knowledge <- function(reference_df,
                                           new_clusters,
                                           umap_obj,
                                           color_attr,
                                           col_vector=NULL,
                                           n_nearest=NULL,
                                           ylim=NULL,
                                           title=NULL,
                                           to_exclude=NULL,
                                           compute_means=FALSE) {

  # Description: given an annotated reference dataframe reference_df,
  # a set of new clusters to analyze in new_clusters, a umap object as run
  # with add_to_annotated_reference, which contains the new clusters to
  # analyze, and a given annotation label, it returns information about the
  # nearest neighbors for each cluster and as a whole for that given
  # annotation label. If remove_low_confidence=TRUE, it doesn't consider in the
  # analysis low confidence neighbors, as defined with the distance_threshold
  # parameter. You can provide a col_vector palette for the final barplot or
  # boxplot.

  if('umap_obj' %in% names(umap_obj)) {
    umap_obj <- umap_obj$umap_obj
  }

  if(!is.null(umap_obj$refined_network)) {
    index_matrix <- umap_obj$umap_out$refined_network$indexes
    distance_matrix <- umap_obj$umap_out$refined_network$distances
  } else {
    index_matrix <- umap_obj$umap_out$knn$indexes
    distance_matrix <- umap_obj$umap_out$knn$distances
  }

  if(!is.null(to_exclude)) {
    idxs_to_remove <- which(reference_df[,color_attr] %in% to_exclude)
  } else {
    idxs_to_remove <- NULL
  }

  n_nearest <- min(n_nearest, ((dim(distance_matrix)[2])-1))

  if(n_nearest > 0.9*dim(reference_df)[1] & compute_means == FALSE) {
    message('The number of nearest neighbours
            (computed as min(n_nearest, ((dim(distance_matrix)[2])-1)))
            is higher than 90% of the total clusters. It should be better to
            compute the means instead of the sum (compute_means = TRUE), to
            avoid biases of the base structure of the clusters')
  }

  annotation_vector <- list()
  weights_list <- list()
  all_distribution <- c()
  all_weights <- c()

  if(is.character(reference_df[,color_attr])) {
    summary_info <- list()
    summary_info_percentage <- list()
    all_info <- c()
  }

  min_n_nearest <- c()
  for(c in new_clusters) {
    if(!c %in% rownames(index_matrix)) {
      c <- paste0('New-', c)
      new_clusters[which(new_clusters == c)] <- paste0('New-', c)
    }
    all_n <- rownames(index_matrix)[index_matrix[c,seq(2,dim(index_matrix)[2])]]
    idxs <- which(!(all_n %in% new_clusters))
    min_n_nearest <- c(min_n_nearest, length(idxs))
  }

  n_nearest <- min(min(min_n_nearest), n_nearest)

  clusters_to_exclude <- c(new_clusters, rownames(reference_df)[idxs_to_remove])

  for(c in new_clusters) {

    dist <- distance_matrix[c,seq(2,dim(distance_matrix)[2])]
    all_n <- rownames(index_matrix)[index_matrix[c,seq(2,dim(distance_matrix)[2])]]
    idxs <- which(!(all_n %in% clusters_to_exclude))
    all_n <- all_n[idxs]
    dist <- dist[idxs]
    all_n <- all_n[seq(1,n_nearest)]
    dist <- dist[seq(1,n_nearest)]

    annotated_reference_f <- reference_df[all_n[which(all_n %in% rownames(reference_df))],]

    names(dist) <- annotated_reference_f[,color_attr]

    if(is.character(reference_df[,color_attr])) {
      annotations <- value_table(dist,
                                 perc = FALSE,
                                 reciprocal = TRUE,
                                 compute_means = compute_means)
      summary_info[[c]] <- annotations
      annotations_perc <- value_table(dist,
                                      perc = TRUE,
                                      reciprocal = TRUE,
                                      compute_means = compute_means)
      summary_info_percentage[[c]] <- annotations_perc
    }

    annotation_vector[[c]] <- annotated_reference_f[,color_attr]
    weights_list[[c]] <- 1/dist
    all_weights <- c(all_weights, 1/dist)
    all_distribution <- c(all_distribution, annotated_reference_f[,color_attr])
  }

  if(is.numeric(reference_df[,color_attr])) {
    summary_info <- lapply(seq(1,length(annotation_vector)), function(i) {
      weighted_summary(annotation_vector[[i]], weights_list[[i]])
    })

    names(summary_info) <- names(annotation_vector)

    out <- S4Vectors::SimpleList('Whole_Distribution' = all_distribution,
                                 'Whole_Summary' = weighted_summary(v = all_distribution,
                                                                    weights = all_weights),
                                 'Cluster_Distribution' = annotation_vector,
                                 'Cluster_Summary' = summary_info,
                                 'Number_of_NN' = n_nearest)

    graphics::boxplot(all_distribution,
            outline = FALSE,
            col = Polychrome::createPalette(1, c(grDevices::blues9[5])),
            ylab = color_attr,
            xlab = 'All',
            ylim=ylim,
            main=title)

    out[['Whole_Boxplot']] <- grDevices::recordPlot()

    if(is.null(col_vector)) {
      col_vector <- Polychrome::createPalette(length(annotation_vector),
                                  c("#ff0000", "#00ff00", "#0000ff"))
      names(col_vector) <- names(annotation_vector)
    }

    graphics::boxplot(annotation_vector,
            outline = FALSE,
            col = col_vector,
            las = 2,
            ylim=ylim,
            main=title,
            ylab = color_attr)

    out[['Cluster_Boxplot']] <- grDevices::recordPlot()

    out <- out[c('Whole_Distribution',
                 'Whole_Summary',
                 'Whole_Boxplot',
                 'Cluster_Distribution',
                 'Cluster_Summary',
                 'Cluster_Boxplot',
                 'Number_of_NN')]

    return(out)
  }

  if(is.character(reference_df[,color_attr])) {

    m <- get_matrix_from_list(summary_info_percentage)

    h1 <- ComplexHeatmap::Heatmap(get_matrix_from_list(summary_info_percentage),
                  col = grDevices::blues9,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  name = 'Neighborhood score\nAnnotation')

    if(is.null(col_vector)) {
      col_vector <- Polychrome::createPalette(ncol(m),
                                  c("#ff0000", "#00ff00", "#0000ff"))
      names(col_vector) <- colnames(m)
    }

    df <- reshape2::melt(m)

    df$Var2 <- factor(df$Var2, levels = gtools::mixedsort(unique(as.vector(df$Var2)), decreasing = TRUE))

    plot <- ggplot2::ggplot(df, ggplot2::aes(x = Var1,
                           y = value,
                           fill = Var2)) +
      ggplot2::geom_bar(stat = "identity")    +
      ggplot2::labs(x = "",
           y = "Values",
           fill = color_attr,
           title = title) +
      ggplot2::theme_minimal() +
      ggplot2::scale_fill_manual(values = col_vector[match(df$Var2, names(col_vector))]) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))

    return(S4Vectors::SimpleList('Whole_Labels' = all_distribution,
                                 'Whole_Frequency' = value_table(all_weights),
                                 'Whole_Percentage' = value_table(all_weights,
                                                                  perc = TRUE,
                                                                  reciprocal = TRUE),
                                 'Cluster_Labels' = annotation_vector,
                                 'Cluster_Frequency' = summary_info,
                                 'Cluster_Percentage' = summary_info_percentage,
                                 'Heatmap-Annotation' = h1,
                                 'Barplot-Annotation' = plot,
                                 'Number_of_NN' = n_nearest))
  }

  return(annotation_vector)

}
