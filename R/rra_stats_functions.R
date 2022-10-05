# Helper functions -------------------------------------------------------------

.getAllEcdfs <- function(back_dis, g2deg){
    all_ecdfs <- apply(back_dis, 2, function(x){
        stats::ecdf(x)
    })
    names(all_ecdfs) <- sort(unique(g2deg))
    all_ecdfs
}

.runBetaReal <- function(x, rr_v){
    re <- sort(rr_v[which(x!=0)], method="quick")
    ps <- stats::pbeta(re, 1:length(re), length(re) - 1:length(re) + 1)
    min(ps, na.rm=T)
}

.runBetaPval <- function(x, score, all_ecdfs){
    all_ecdfs[[as.character(sum(x!=0))]](score)
}

.runBetaRan <- function(ll, n_iter){
    ran <- .colSort(.vRunif(rep(ll, n_iter)))
    pb  <- stats::pbeta(ran, 1:ll, ll - 1:ll + 1)
    mm  <- apply(pb, 2, min)
    return(mm)
}

# Main functions ---------------------------------------------------------------

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
getRandomBetas <- function(degs, n_iter){

    sapply(sort(unique(degs)), function(ll){

        .runBetaRan(ll, n_iter)

    }, simplify=T)

}


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
getRandomBetasPar <- function(degs, n_iter, cores=1){

    do.call(cbind,
            parallel::mclapply(sort(unique(degs)), function(ll){

        .runBetaRan(ll, n_iter)

    }, mc.cores=cores))

}

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




















