#' @export
runProde <- function(prodeInput, cores=1, filterCtrl = T, n_iter=10000){

    start <- Sys.time()
    message("\n\nRunning ProDe on ", nrow(assay(prodeInput)), " genes!\n\n")
    message("[1] Running Linear models fit \t\t\t", Sys.time(), "\n")

    beta_tab <- fitLms( # fast fit linear models to each gene + covariates
        x        = designMatrix(prodeInput),
        y        = SummarizedExperiment::assay(prodeInput)
    )

    filtered <-S4Vectors::DataFrame()
    if (filterCtrl){
        keep <- which(beta_tab[,"wtm"] < 0)
        filtered <- beta_tab[-keep,]
        beta_tab <- beta_tab[ keep,]
    }

    message("[2] Subsetting adj matrix and NodesDegree \t", Sys.time(), "\n")

    adj_m <- filterAdjMatrix(prodeInput, rownames(beta_tab))

    stopifnot(.checkBetAdj(beta_tab, adj_m))

    message("[3] Compute background random distribution \t", Sys.time(), "\n")

    if (cores > 1){

        back_dis <- getRandomBetasPar(
            adj_m  = adj_m,
            n_iter = n_iter,
            cores  = cores
        )

    } else {

        back_dis <- getRandomBetas(
            adj_m  = adj_m,
            n_iter = n_iter
        )

    }

    message("[4] Compute RRA statistics and p-value \t\t", Sys.time(), "\n")

    real_betas <- getRealBetas(
        bet_tab  = beta_tab,
        adj_m    = adj_m,
        back_dis = back_dis
    )

    ## 4. Collect final output .................................................

    end <- Sys.time()
    message("[5] Done !!! \t\t\t\t\t", end-start)

    out <- S4Vectors::DataFrame(
        "gene" = rownames(beta_tab),
        beta_tab,
        real_betas,
        "node_degree" = c(Matrix::rowSums(adj_m)[rownames(beta_tab)])
    )

    prodeResults(
        out,
        adj_m,
        S4Vectors::DataFrame(filtered)
    )

}

