# Main Prode functions ---------------------------------------------------------

.fitLms <- function(x, y, covs=NULL){

    yy <- t(y)

    dat <- data.frame(
        yy, x
    )

    formula <- stats::as.formula(
        paste0(
            "yy~",
            paste(
                covs,
                collapse="+"
            )
        )
    )

    fits <- summary(stats::lm(
        formula,
        data=dat
    ))

    coefs  <- lapply(fits, function(fit)
        fit[[4]][1,]
    )

    mm <- do.call(rbind, coefs)
    mm <- cbind(mm, "u"=rank(mm[,3])/nrow(mm))
    rownames(mm) <- colnames(yy)
    return(mm)

}

getRandomBetas <- function(degs, n_iter){

    sapply(sort(unique(degs)), function(ll){

        .runBetaRan(ll, n_iter)

    }, simplify=T)

}


getRandomBetasPar <- function(degs, n_iter, cores=1){

    do.call(rbind,
            parallel::mclapply(sort(unique(degs)), function(ll){

        .runBetaRan(ll, n_iter)

    }, mc.cores=cores))

}

getRealBetas <- function(rr_v, adj_mat, back_dis, g2deg){

        all_ecdfs <- .getAllEcdfs(back_dis, g2deg)

        pbs <- apply(adj_mat, 1, function(x){

            score <- .runBetaReal(x, rr_v)
            p.val <- .runBetaPval(x, score, all_ecdfs)

            c(score, p.val)

        })

        cbind(
            "rra_score" = pbs[1,],
            "p.value"   = pbs[2,],
            "fdr"       = stats::p.adjust(pbs[2,], "fdr")
        )
}

# Wrapper ----------------------------------------------------------------------

runProde <- function(ds, dm, covs, gr, cores=1){

    start <- Sys.time()
    ## 1. Fit linear models ....................................................
    cat(paste0("\n\nRunning ProDe on ", nrow(ds), " genes!\n\n"))
    cat(paste0("[1] Running Linear models fit \t\t\t", Sys.time(), "\n"))

    beta_tab <- .fitLms( # fast fit linear models to each gene + covariates
        x    = dm,
        y    = ds,
        covs = covs
    )

    ## 2. Retrieve graph-data ..................................................

    cat(paste0("[2] Collect graph data \t\t\t\t", Sys.time(), "\n"))

    adj_mat <- .getAdjMatr( # get ordered and filtered adj matrix
        etab = gr,
        gns  = rownames(beta_tab)
    )

    g2deg <- .getDegLookup( # get lookup for gene2degree in graph
        adj_mat = adj_mat
    )

    ## 3. Compute new RRA probabilities ........................................

    cat(paste0("[3] Compute background random distribution \t", Sys.time(), "\n"))

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

    cat(paste0("[4] Compute RRA statistics and p-value \t\t", Sys.time(), "\n"))

    real_betas <- getRealBetas(
        rr_v     = beta_tab[,"u"],
        adj_mat  = adj_mat,
        back_dis = back_dis,
        g2deg    = g2deg
    )

    ## 4. Collect final output .................................................

    end <- Sys.time()
    cat(paste0("[5] Done !!! \t\t\t\t\t", end-start))

    out <- data.frame(
        "gene" = rownames(beta_tab),
        beta_tab,
        real_betas,
        "node_degree" = c(g2deg[rownames(beta_tab)])
    )

    return(out)

}


















