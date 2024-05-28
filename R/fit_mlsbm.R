#' R/Rcpp function for fitting multilevel stochastic block model
#'
#' This function allows you to fit multilevel stochastic block models.
#' @param A An adjacency list of length L, the number of levels. Each level contains an n x n symmetric adjacency matrix.
#' @param K The number of clusters specified a priori.
#' @param z_init Initialized cluster indicators. If NULL, will initialize automatically with Louvain algorithm.
#' @param a0 Dirichlet prior parameter for cluster sizes for clusters 1,...,K.
#' @param b10 Beta distribution prior paramter for community connectivity.
#' @param b20 Beta distribution prior parameter for community connectivity.
#' @param n_iter The number of total MCMC iterations to run.
#' @param burn The number of burn-in MCMC iterations to discard. The number of saved iterations will be n_iter - burn.
#' @param verbose Whether to print a progress bar to track MCMC progress. Defaults to true.
#' @param r Resolution parameter for Louvain initialization. Sould be >= 0 and higher values give a larger number of smaller clusters.
#' @keywords SBM MLSBM Gibbs Bayesian networks 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain
#' @importFrom bluster mergeCommunities
#' @export
#' @return A list of MCMC samples, including the MAP estimate of cluster indicators (z)
#' @examples
#' data(AL)
#' # increase n_iter in practice
#' fit <- fit_mlsbm(AL,3,n_iter = 100)
fit_mlsbm <- function(A,
                      K,
                      z_init = NULL,
                      a0 = 2,
                      b10 = 1,
                      b20 = 1,
                      n_iter = 1000,
                      burn = 100,
                      verbose = FALSE,
                      r = 1.2)
{
    # Initialize parameters
    n = dim(A[[1]])[1] # number of nodes
    K0 = K # putative number of clusters
    if(is.null(z_init)) # initialize clusters?
    {
        print("Initializing using igraph")
        # use igraph
        G = igraph::graph_from_adjacency_matrix(A[[1]],mode = "undirected",diag = FALSE)
        
        # check resolutions for inits
        initialized = FALSE
        while(!initialized)
        {
            fit_init = igraph::cluster_louvain(G, resolution = r)
            zinit = fit_init$membership
            K_found = length(unique(zinit))
            if(K_found > K)
            {
                initialized = TRUE
            }
            else
            {
                r = r*1.05
            }
        }
        
        zinit = as.numeric(bluster::mergeCommunities(G,as.factor(zinit),number = K0))
        zs = remap_canonical2(zinit)
    }
    else
    {
        zs = remap_canonical2(as.numeric(z_init))
        zinit = z_init
    }
    ns = table(zs) # initial cluster counts
    pis = ns/n # initial cluster proportions
    Ps = array(0.4,c(K0,K0)) # initial connectivity matrix
    
    n_sim = n_iter - burn # number of stored simulations
    Z = array(0,c(n_sim,n)) # storage for cluster assignments
    PI = array(0,c(n_sim,K0)) # storage for cluster proportions
    PM = array(0,c(n_sim,K0,K0)) # storage for community connection params
    draw_logf_z_given_A = rep(0,n_sim) # P(z|Data) for MAP estimate of z
    
    start.time<-proc.time()
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    for (i in 1:n_iter){
        # TODO: use pointers
        # Step 1. pi
        as = update_counts(zs,K0) + a0
        pis = rdirichlet_cpp(1, as)
        
        # Step 2. z
        zs = update_z(zs,A,Ps,pis,1:K0)
        if(verbose) print(table(zs))
        if(any(update_counts(zs,K0) < 10))
        {
            zs = zinit
        }
        
        # Step 3. P
        Ps = update_P(A,zs,K0,b10,b20)
        
        # Calculating P(z|Data) to find MAP
        lz = logf_z_given_A(zs,A,pis,Ps)
        
        # Store values
        if(i > burn)
        {
            j = i - burn
            Z[j,] = zs
            PI[j,] = pis
            PM[j,,] = Ps
            draw_logf_z_given_A[j] = lz
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)
    run.time<-proc.time()-start.time
    cat("Finished MCMC after",run.time[1],"seconds")
    
    ret_list <- list(A = A,
                     K = K,
                     Z = Z,
                     PI = PI,
                     PM = PM,
                     logf = draw_logf_z_given_A,
                     z = Z[which.max(draw_logf_z_given_A),],
                     P = PM[which.max(draw_logf_z_given_A),,],
                     pi = PI[which.max(draw_logf_z_given_A),])
    class(ret_list) <- "MLSBM"
    return(ret_list)
}
