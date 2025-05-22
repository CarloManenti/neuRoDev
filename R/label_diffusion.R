#' Label diffusion K-medoids like
#'
#' @param umap_obj The starting umap object
#' @param labels The starting labels vector (named with the vertex/cluster names
#' if you suspect that the network stored in the umap object doesn't have the
#' same order of vertexes as those that are in the labels vector)
#' @param max_n_iter Maximum number of iterations
#' @param max_n_equal The number of times the previous association has to be
#' equal (to a percentage dictated by `perc_equal`) to the current association
#' to stop the diffusion
#' @param perc_equal The percentage of same association between the different
#' iterations to be considered equal
#' @param n_nearest The number of nearest neighbours to consider for
#' classification. It defaults to the minimum between the number of labels / 100
#' and half of the label group with the lowest number of components.
#'
#' @return The diffused labels
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' M <- reference_signatures_correlation(S, refS)
#' umap_obj <- umap_graph_clustering(M)
#' graph <- umap_obj$refined_network
#' igraph::V(graph)$label <- paste0('Group', rep(seq(1,2), each = 5))
#' label_diffusion(umap_obj = umap_obj,
#' labels = paste0('Group', rep(seq(1,2), each = 5)))
label_diffusion <- function(umap_obj,
                            labels,
                            max_n_iter = 1000,
                            max_n_equal = 10,
                            perc_equal = 0.99,
                            n_nearest = NULL) {

  if(is.null(n_nearest)) {
    n_nearest <- min(round(length(labels)/100), ceiling(min(table(labels))/2))
  }

  if('refined_network' %in% names(umap_obj)) {
    graph <- umap_obj$refined_network
  } else {
    graph <- umap_obj$umap_knn_igraph
  }

  if(is.null(names(labels))) {
    igraph::V(graph)$label <- labels
  } else {
    igraph::V(graph)$label <- labels[igraph::V(graph)$name]
    labels <- igraph::V(graph)$label
  }

  names(labels) <- igraph::V(graph)$name
  groups <- unique(labels)

  d <- umap_obj$umap_out$refined_network$distances
  ind <- umap_obj$umap_out$refined_network$indexes

  old_best_association <- NULL

  n_equal <- 0

  n <- 1
  while(n <= max_n_iter) {

    best_association <- character()

    for(point in names(labels)) {

      sub_d <- d[point, seq(2,(n_nearest+1))]
      sub_i <- ind[point, seq(2,(n_nearest+1))]
      closest_points <- rownames(d)[sub_i]
      closest_groups <- labels[closest_points]

      names(sub_d) <- closest_groups

      closest_groups_table <- table(closest_groups)

      if(length(closest_groups_table[which(closest_groups_table == max(closest_groups_table))]) > 1) {

        names(sub_d) <- closest_groups

        sub_d <- sub_d[which(names(sub_d) %in% names(closest_groups_table)[which(closest_groups_table == max(closest_groups_table))])]

        if(max(closest_groups_table) == 1) {
          best_association <- c(best_association, names(sub_d)[which.min(sub_d)])
        } else {
          mean_by_group <- get_column_group_average(t(as.matrix(sub_d)),
                                                    names(sub_d))

          if(is.vector(mean_by_group)) {
            mean_by_group <- t(as.matrix(mean_by_group))
          }

          best_association <- c(best_association,
                                apply(mean_by_group,
                                      1,
                                      function(i) {
                                        colnames(mean_by_group)[which.min(i)]
                                      }))
        }
      } else {
        best_association <- c(best_association, names(closest_groups_table)[which.max(closest_groups_table)])
      }

    }

    names(best_association) <- names(labels)

    labels <- best_association

    if(sum(best_association == old_best_association)/length(best_association) >= perc_equal) {

      n_equal <- n_equal + 1

    }

    if(n_equal == max_n_equal) {
      break()
    }

    old_best_association <- best_association

    n <- n + 1

  }

  return(best_association)

}
