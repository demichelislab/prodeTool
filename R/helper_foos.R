# helper functions

.getDegLookup <- function(adj_mat){
    rowSums(adj_mat > 0)
}

.getAdjMatr <- function(etab, gns){
    (unclass(table(etab))[gns,gns] > 0)*1
}

.vRunif <- Vectorize(function(n){
    runif(n)
}, "n")

.colSort <- function(x){
    apply(x, 2, sort, method="quick")
}

.runBetaRan <- function(ll, n_iter){
    ran <- .colSort(.vRunif(rep(ll, n_iter)))
    pb  <- stats::pbeta(ran, 1:ll, ll - 1:ll + 1)
    mm  <- apply(pb, 2, min)
    return(mm)
}

.getAllEcdfs <- function(back_dis, g2deg){
    all_ecdfs <- apply(back_dis, 2, function(x){
        stats::ecdf(x)
    })
    names(all_ecdfs) <- sort(unique(g2deg))
}

.runBetaReal <- function(x, rr_v){
    re <- sort(rr_v[which(x!=0)], method="quick")
    ps <- stats::pbeta(re, 1:length(re), length(re) - 1:length(re) + 1)
    min(ps, na.rm=T)
}

.runBetaPval <- function(x, score, all_ecdfs){
    all_ecdfs[[as.character(sum(x!=0))]](score)
}


