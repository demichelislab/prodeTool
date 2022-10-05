# Wrapper around ProDe call ----------------------------------------------------

#' Runs linear models on input data
#'
#' @param x
#' @param y
#' @param covs
#'
#' @return matrix with summary statistics of each fit
#' @export
#'
#' @examples
runProde <- function(
    scores_matrix, columns_data, condition,
    covariates, edge_table, cores=1, filterCtrl = T
    ){

    start <- Sys.time()
    ## 1. Fit linear models ....................................................
    message("\n\nRunning ProDe on ", nrow(scores_matrix), " genes!\n\n")
    message("[1] Running Linear models fit \t\t\t", Sys.time(), "\n")

    beta_tab <- fitLms( # fast fit linear models to each gene + covariates
        x        = columns_data,
        y        = scores_matrix,
        cond     = condition,
        covs     = covariates,
        filter   = filterCtrl
    )

    ## 2. Retrieve graph-data ..................................................

    message("[2] Collect graph data \t\t\t\t", Sys.time(), "\n")

    adj_mat <- .getAdjMatr( # get ordered and filtered adj matrix
        etab = edge_table,
        gns  = rownames(beta_tab)
    )

    g2deg <- .getDegLookup( # get lookup for gene2degree in graph
        adj_mat = adj_mat
    )

    ## 3. Compute new RRA probabilities ........................................

    message("[3] Compute background random distribution \t", Sys.time(), "\n")

    if (cores > 1){

        back_dis <- getRandomBetasPar(
            degs   = g2deg,
            n_iter = 10000,
            cores  = cores
        )

    } else {

        back_dis <- getRandomBetas(
            degs   = g2deg,
            n_iter = 10000
        )

    }

    message("[4] Compute RRA statistics and p-value \t\t", Sys.time(), "\n")

    real_betas <- getRealBetas(
        rr_v     = beta_tab[,"u"],
        adj_mat  = adj_mat,
        back_dis = back_dis,
        g2deg    = g2deg
    )

    ## 4. Collect final output .................................................

    end <- Sys.time()
    message("[5] Done !!! \t\t\t\t\t", end-start)

    out <- data.frame(
        "gene" = rownames(beta_tab),
        beta_tab,
        real_betas,
        "node_degree" = c(g2deg[rownames(beta_tab)])
    )

    return(out)

}
