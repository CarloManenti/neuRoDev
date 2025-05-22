#' Signatures Initialization
#'
#' @param expression_matrix A single-cell RNAseq matrix or a list of matrices
#' @param algorithm The algorithm to use between 'louvain' and 'leiden' to
#' cluster cells
#' @param split A vector that defines how to split a expression_matrix into multiple ones
#' @param subsample A boolean variable that defines if each of the provided
#' expression_matrixs has to be proportionally subsampled to the same number of cells
#' @param PCs_to_use The number of Principal Components to consider for
#' FindNeighbors of Seurat
#' @param resolution The resolution to consider for FindClusters of Seurat
#' @param clusters_mv A vector that substitutes the clustering
#' (ex. celltypes if already annotated)
#' @param Ntotal The total number of cells to consider in the subsampling
#' @param BOGs A boolean variable to define if map_clusters_to_representation has
#' to compute also the BOGs
#' @param pBOGs A boolean variable to define if map_clusters_to_representation has
#' to compute also the pBOGs
#' @param rEEA A boolean variable to define if map_clusters_to_representation has
#' to compute also the rEEA
#' @param pBOGsig A boolean variable to define if map_clusters_to_representation
#' has to compute also the pBOGsig
#' @param verbose A boolean variable to define if the steps of Seurat will be
#' printed or not
#' @param condaenv The name of the conda environment in which 'leidenalg'
#' module is present, in case the user wants to use the Leiden algorithm
#'
#' @return The result of map_clusters_to_representation run on each expression_matrix after
#' splitting, subsampling, clustering if/when required
#' @export
#'
#' @examples
#' set.seed(123)
#' M <- matrix(runif(2000000,0,100), ncol=20000)
#' rownames(M) <- paste0('Gene-', seq(1, dim(M)[1]))
#' colnames(M) <- paste0('Cell-', seq(1, dim(M)[2]))
#' getClusterSignatures(M, resolution=2)
getClusterSignatures <- function(expression_matrix,
                           algorithm='louvain',
                           split=NULL,
                           subsample=FALSE,
                           PCs_to_use=NULL,
                           resolution=1,
                           clusters_mv=NULL,
                           Ntotal=2000,
                           BOGs=FALSE,
                           pBOGs=FALSE,
                           rEEA=FALSE,
                           pBOGsig=FALSE,
                           verbose=FALSE,
                           condaenv='r-reticulate') {

  if(!is.list(expression_matrix)) {
    expression_matrix <- Matrix::Matrix(expression_matrix, sparse = TRUE)
    expression_matrix <- list(expression_matrix)
    if(!is.null(split)) {
      split <- list(split)
    }
    if(!is.null(clusters_mv)) {
      clusters_mv <- list(clusters_mv)
    }
  }

  if(!is.null(split)) {
    expression_matrix <- lapply(seq_len(length(expression_matrix)), function(i) {
      matrix_split(expression_matrix[[i]], split[[i]])
    })
    expression_matrix <- unlist(expression_matrix, recursive = FALSE)
  }

  if(is.null(clusters_mv)) {
    mvs <- lapply(expression_matrix, function(i) {
      s <- Seurat::CreateSeuratObject(counts = i)
      s <- Seurat::NormalizeData(object = s,
                                 verbose = verbose)

      span <- 0.3
      span_x <- 0.3
      warns <- TRUE

      while(warns) {
        tryCatch(Seurat::FindVariableFeatures(object = s,
                                              verbose = verbose,
                                              span = span), warning = function(x) {
            span_x <<- span + 0.1
            })
        if(span_x != span) {
          span <- span_x
        } else {
          warns <- FALSE
        }
      }

      s <- Seurat::FindVariableFeatures(object = s,
                                        verbose = verbose,
                                        span = span)
      s <- Seurat::ScaleData(object = s,
                             features = Seurat::VariableFeatures(s),
                             verbose = verbose)
      numPCs <- min(50, dim(i)[2]-1)
      s <- Seurat::RunPCA(object = s,
                          features = Seurat::VariableFeatures(s),
                          npcs = numPCs,
                          approx = FALSE,
                          verbose = verbose)

      if(is.null(PCs_to_use)) {
        PCs_to_use <- (s$pca@stdev)^2
        PCs_to_use <- PCs_to_use/sum(PCs_to_use)
        PCs_to_use <- cumsum(PCs_to_use)
        PCs_to_use <- min(which(PCs_to_use >= 0.75))
      }

      PCs_to_use <- min(numPCs, PCs_to_use)

      s <- Seurat::FindNeighbors(object = s,
                                 dims = seq_len(PCs_to_use),
                                 verbose = verbose)

      if(algorithm == 'leiden') {
        reticulate::use_condaenv(condaenv)
      }

      s <- Seurat::FindClusters(object = s,
                                resolution = resolution,
                                verbose = verbose,
                                algorithm = algorithm)

      mv <- s$seurat_clusters
      mv <- as.numeric(as.vector(mv))
      names(mv) <- colnames(i)
      return(mv)
    })

    if(subsample) {

        Ntotal <- min(unlist(lapply(expression_matrix, function(i) {dim(i)[2]})), Ntotal)

        new_expression_matrix <- lapply(seq_len(length(expression_matrix)), function(i) {
            Proportional_sampling(expression_matrix[[i]], mvs[[i]], Ntotal=Ntotal)
        })

        new_mvs <- lapply(new_expression_matrix, function(i) {
          s <- Seurat::CreateSeuratObject(counts = i)
          s <- Seurat::NormalizeData(object = s,
                                     verbose = verbose)

          span <- 0.3
          span_x <- 0.3
          warns <- TRUE

          while(warns) {
            tryCatch(Seurat::FindVariableFeatures(object = s,
                                                  verbose = verbose,
                                                  span = span), warning = function(x) {
                                                    span_x <<- span + 0.1
                                                  })
            if(span_x != span) {
              span <- span_x
            } else {
              warns <- FALSE
            }
          }

          s <- Seurat::FindVariableFeatures(object = s,
                                            verbose = verbose,
                                            span = span)
          s <- Seurat::ScaleData(object = s,
                                 features = Seurat::VariableFeatures(s),
                                 verbose = verbose)
          numPCs <- min(50, dim(i)[2]-1)
          s <- Seurat::RunPCA(object = s,
                              features = Seurat::VariableFeatures(s),
                              npcs = numPCs,
                              approx = FALSE,
                              verbose = verbose)

          if(is.null(PCs_to_use)) {
            PCs_to_use <- (s$pca@stdev)^2
            PCs_to_use <- PCs_to_use/sum(PCs_to_use)
            PCs_to_use <- cumsum(PCs_to_use)
            PCs_to_use <- min(which(PCs_to_use >= 0.75))
          }

          PCs_to_use <- min(numPCs, PCs_to_use)

          s <- Seurat::FindNeighbors(object = s,
                                     dims = seq_len(PCs_to_use),
                                     verbose = verbose)

          s <- Seurat::FindClusters(object = s,
                                    resolution = resolution,
                                    verbose = verbose,
                                    algorithm = algorithm)
          mv <- s$seurat_clusters
          mv <- as.numeric(as.vector(mv))
          names(mv) <- colnames(i)
          return(mv)
        })

      summary <- lapply(seq_len(length(new_expression_matrix)), function(i) {
          map_clusters_to_representation(M = new_expression_matrix[[i]],
                                         group = new_mvs[[i]],
                                         BOGs = BOGs,
                                         pBOGs = pBOGs,
                                         rEEA = rEEA,
                                         pBOGsig = pBOGsig)
      })

      summary <- lapply(seq(1,length(summary)), function(i) {
          s <- summary[[i]]
          s$Membership <- new_mvs[[i]]
          return(s)
      })

    } else {

      summary <- lapply(seq_len(length(expression_matrix)), function(i) {
        if(length(unique(mvs[[i]])) == 1) {
          return(S4Vectors::List(group_counts = mvs[[i]],
                                 P = as.matrix(Matrix::rowMeans(expression_matrix[[i]])),
                                 S = as.matrix(Matrix::rowMeans(expression_matrix[[i]])),
                                 BOG = NULL, DE_out = NULL,
                                 pBOG = NULL,
                                 pBOGs = NULL,
                                 clust_rEEA = NULL))
        }
        map_clusters_to_representation(expression_matrix[[i]],
                                       mvs[[i]],
                                       BOGs = BOGs,
                                       pBOGs = pBOGs,
                                       rEEA = rEEA,
                                       pBOGsig = pBOGsig)
      })

      summary <- lapply(seq(1,length(summary)), function(i) {
        s <- summary[[i]]
        s$Membership <- mvs[[i]]
        return(s)
      })

    }

  } else {

    if(subsample) {

      Ntotal <- min(unlist(lapply(expression_matrix, function(i) {dim(i)[2]})), Ntotal)

      new_expression_matrix <- lapply(seq_len(length(expression_matrix)), function(i) {
        Proportional_sampling(expression_matrix[[i]], clusters_mv[[i]], Ntotal=Ntotal)
      })

      new_clusters_mv <- lapply(seq_len(length(expression_matrix)), function(i) {
        idxs <- Proportional_sampling(expression_matrix[[i]],
                                      clusters_mv[[i]],
                                      Ntotal=Ntotal,
                                      ReturnSCE = FALSE)$idx
        clusters_mv[[i]][idxs]
      })

      summary <- lapply(seq_len(length(new_expression_matrix)), function(i) {
        if(length(unique(mvs[[i]])) == 1) {
          return(S4Vectors::List(group_counts = mvs[[i]],
                                 P = as.matrix(Matrix::rowMeans(expression_matrix[[i]])),
                                 S = as.matrix(Matrix::rowMeans(expression_matrix[[i]])),
                                 BOG = NULL,
                                 DE_out = NULL,
                                 pBOG = NULL,
                                 pBOGs = NULL,
                                 clust_rEEA = NULL))
        }
        map_clusters_to_representation(new_expression_matrix[[i]],
                                       new_clusters_mv[[i]],
                                       BOGs = BOGs,
                                       pBOGs = pBOGs,
                                       rEEA = rEEA,
                                       pBOGsig = pBOGsig)
      })

      summary <- lapply(seq(1,length(summary)), function(i) {
        s <- summary[[i]]
        s$Membership <- new_clusters_mv[[i]]
        return(s)
      })

    } else {

      summary <- lapply(seq_len(length(expression_matrix)), function(i) {
        if(length(unique(clusters_mv[[i]])) == 1) {
          return(S4Vectors::List(group_counts = clusters_mv[[i]],
                                 P = as.matrix(Matrix::rowMeans(expression_matrix[[i]])),
                                 S = as.matrix(Matrix::rowMeans(expression_matrix[[i]])),
                                 BOG = NULL,
                                 DE_out = NULL,
                                 pBOG = NULL,
                                 pBOGs = NULL,
                                 clust_rEEA = NULL))
        }
        map_clusters_to_representation(expression_matrix[[i]],
                                       clusters_mv[[i]],
                                       BOGs = BOGs,
                                       pBOGs = pBOGs,
                                       rEEA = rEEA,
                                       pBOGsig = pBOGsig)
      })

      summary <- lapply(seq(1,length(summary)), function(i) {
        s <- summary[[i]]
        s$Membership <- clusters_mv[[i]]
        return(s)
      })

    }
  }

  if(length(summary) == 1) {
    summary <- summary[[1]]
  }

  return(summary)
}
