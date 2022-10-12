# Helper functions -------------------------------------------------------------

.runBetaRan1 <- function(ll, n_iter, nn){
    ran <- matrix(.vRunif(rep(ll, n_iter)), nrow=1)
    #ran <- matrix(.vSample(nn, rep(ll, n_iter)), nrow=1)/nn
    pb  <- stats::pbeta(ran, 1:ll, ll - 1:ll + 1)
    return(t(pb))
}

.runBetaRan <- function(ll, n_iter, nn){
    ran <- .colSort(.vRunif(rep(ll, n_iter)))
    #ran <- .colSort(.vSample(nn, rep(ll, n_iter)))/nn
    pb  <- stats::pbeta(ran, 1:ll, ll - 1:ll + 1)
    mm  <- apply(pb, 2, min)
    return(mm)
}

# .getAllEcdfs <- function(back_dis, g2deg){
#     all_ecdfs <- apply(back_dis, 2, function(x){
#         stats::ecdf(x)
#     })
#     names(all_ecdfs) <- sort(unique(g2deg))
#     all_ecdfs
# }

.runBetaReal <- function(x, rr_v){
    re <- sort(rr_v[which(x!=0)], method="quick")
    ps <- stats::pbeta(re, 1:length(re), length(re) - 1:length(re) + 1)
    min(ps, na.rm=T)
}

# .runBetaPval <- function(x, score, all_ecdfs){ # one can also do mean > x
#     all_ecdfs[[as.character(sum(x!=0))]](score)
# }

.runBetaPval <- function(x, score, back_dis){ # empricial p-valeu
    mean(back_dis[,as.character(sum(x!=0))] <= score)
}

# Main functions ---------------------------------------------------------------

#' Runs linear models on input data
#'
#' @param degs
#' @param n_iter
#' @param nn
#'
#' @return matrix with summary statistics of each fit
#' @export
#'
#' @examples
getRandomBetas <- function(degs, n_iter, nn){

    sapply(sort(unique(degs)), function(ll){
        if (ll == 1){
            .runBetaRan1(ll, n_iter, nn)
        } else {
            .runBetaRan(ll, n_iter, nn)
        }
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
getRandomBetasPar <- function(degs, n_iter, nn, cores=1){

    do.call(cbind,
            parallel::mclapply(sort(unique(degs)), function(ll){

        if (ll == 1){
            .runBetaRan1(ll, n_iter, nn)
        } else {
            .runBetaRan(ll, n_iter, nn)
        }

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

        colnames(back_dis) <- sort(unique(g2deg))
        pbs <- apply(adj_mat, 1, function(x){

            score <- .runBetaReal(x, rr_v)
            p.val <- .runBetaPval(x, score, back_dis)

            c(score, p.val)

        })

        cbind(
            "rra_score" = pbs[1,],
            "p.value"   = pbs[2,],
            "fdr"       = stats::p.adjust(pbs[2,], "fdr")
        )
}




















