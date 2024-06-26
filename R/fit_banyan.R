#' Fit Banyan community detection model for spatially resolved single cell data
#'
#' This function allows you to fit the Bayesian multi-layer SBM for integration of spatial and gene expression data for identifying cell sub-populations
#' @param seurat_obj A Seurat object with PCA reduction and spatial coordinates. If provided, the exp and coords arguments are ignored
#' @param labels User-defined tissue architecture labels as a length-n vector. If provided, only posterior of ACC parameters will be provided and tissue architecture identification will be skipped.
#' @param exp A binary adjacency matrix encoding the gene expression network. Not used if seurat_obj is provided. 
#' @param coords_df A matrix or data frame with rows as cells and 2 columns for coordinates. Rows should be ordered the same as in exp. Not used if seurat_obj is provided. 
#' @param z_init An optional initialization for cluster indicators. Ignored if seurat_obj is provided.
#' @param K The number of sub-populations to infer
#' @param n_pcs The number of principal components to use from the Seurat object
#' @param R A length 2 vector of integers for the number of neighbors to use. 1st element corresponds to the number of neighbors in gene expression network and 2nd element for spatial.
#' @param a0 Dirichlet prior parameter (shared across all communities)
#' @param b10 Beta prior number of connections
#' @param b20 Beta prior number on non-edges
#' @param n_iter The number of total MCMC iterations to run.
#' @param burn The number of MCMC iterations to discard as burn-in. A total of n_iter - burn iterations will be saved.
#' @param verbose Whether or not to print cluster allocations at each iteration
#' @param s Louvain resolution parameter
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @importFrom Seurat GetTissueCoordinates
#' @export
#' @return A list of MCMC samples, including the MAP estimate of cluster indicators (z)
#' 
fit_banyan <- function(seurat_obj = NULL,
                       labels = NULL,
                       exp = NULL,
                       coords_df = NULL,
                       z_init = NULL,
                       K,
                       n_pcs = 16,
                       R = NULL,
                       a0 = 2,
                       b10 = 1,
                       b20 = 1,
                       n_iter = 1000,
                       burn = 100, 
                       verbose = TRUE,
                       s = 1.2)
{
  
  # user provides a Seurat object
  if(!is.null(seurat_obj))
  {
    # check PCA reductions exist in the Seurat object
    if(!is.null(seurat_obj@reductions$pca))
    {
      exp <- seurat_obj@reductions$pca@cell.embeddings[,1:n_pcs]
      
    }
    else
    {
      return("Error: No PCA reductions found in the supplied Seurat object. Try passing features through using the exp argument, or add PCA reductions to seurat_obj$reductions$pca.")
    }
    # check coordinates exist in the Seurat object
    if(length(seurat_obj@images) != 1)
    {
      return("Error: Please provide Seurat object with exactly 1 image slot to access spatial coordinates from.")
    }
    else
    {
      coords<-as.data.frame(GetTissueCoordinates(seurat_obj))
      coords <- coords[,c(2,1)]
      n <- nrow(coords)
    }
    # set number of neighbors
    if(is.null(R))
    {
      r1 <- r2 <- round(sqrt(n))
    }
    else
    {
      # user provided some invalid but non-null # of neighbors
      if(length(R) != 2)
      {
        return("Error: Please set R to either NULL or a length 2 vector of integer number of neighbors for each layer.")
      }
      else
      {
        # number of neighbors for gene expression
        r1 <- round(R[1])
        # number of neighbors for spatial
        r2 <- round(R[2])
      }
    }
    # build networks
    # gene expression
    exp <- as.matrix(exp)
    A1 <- build_knn_graph(exp,r1)
    # spatial
    coords <- as.matrix(coords)
    A2 <- build_knn_graph(coords,r2)
    
    # compile into multi-layer network
    AL <- list(A1,A2)
    zinit <- z_init   
    

  }
  else
  {
    n <- nrow(coords_df)
    
    # set number of neighbors
    if(is.null(R))
    {
      r1 <- r2 <- round(sqrt(n))
    }
    else
    {
      # user provided some invalid but non-null # of neighbors
      if(length(R) != 2)
      {
        return("Error: Please set R to either NULL or a length 2 vector of integer number of neighbors for each layer.")
      }
      else
      {
        # number of neighbors for gene expression
        r1 <- round(R[1])
        # number of neighbors for spatial
        r2 <- round(R[2])
      }
    }
    A1 <- exp
    A2 <- build_knn_graph(coords_df,r2)
    # compile into multi-layer network
    AL <- list(A1,A2)
    coords <- coords_df
    zinit <- z_init
    
  }
  
  if(!is.null(labels))
  {
    message("Computing ACC parameters based on provided tissue architecture labels")
    fit <- get_ACC(z = labels,
                   AL = AL, 
                   b10 = b10, 
                   b20 = b20)
  }
  else
  {
    fit <- fit_mlsbm(A = AL,
                            K = K,
                            z_init = zinit,
                            a0 = a0,
                            b10 = b10,
                            b20 = b20, 
                            n_iter = n_iter,
                            burn = burn,
                            verbose = verbose,
                            r = s)
    # return coordinates for plotting
    fit$coords = as.data.frame(coords)
    
    # get BIC
    lls = fit$logf
    K = fit$K
    n = length(lls)
    K_param = choose(K,2) + K + n
    bic =  (max(lls) + 2*K_param*log(n))
    fit$bic = bic 
  }
  
  # return as banyan object
  class(fit) <- "banyan"
  return(fit)
}
