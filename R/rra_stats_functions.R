# Helper functions -------------------------------------------------------------

.runBetaRan1 <- function(ll, n_iter){
    ran <- matrix(.vRunif(rep(ll, n_iter)), nrow=1)
    #ran <- matrix(.vSample(nn, rep(ll, n_iter)), nrow=1)/nn
    pb  <- stats::pbeta(ran, 1:ll, ll - 1:ll + 1)
    return(t(pb))
}

.runBetaRan <- function(ll, n_iter){
    ran <- .colSort(.vRunif(rep(ll, n_iter)))
    #ran <- .colSort(.vSample(nn, rep(ll, n_iter)))/nn
    pb  <- stats::pbeta(ran, 1:ll, ll - 1:ll + 1)
    mm  <- apply(pb, 2, min)
    return(mm)
}

.runBetaReal <- function(x, rr_v){
    re <- sort(rr_v[which(x!=0)], method="quick")
    ps <- stats::pbeta(re, 1:length(re), length(re) - 1:length(re) + 1)
    min(ps, na.rm=T)
}

.runBetaPval <- function(x, score, back_dis){ # empricial p-valeu
    mean(back_dis[,as.character(sum(x!=0))] <= score)
}

# Main functions ---------------------------------------------------------------


getRandomBetas <- function(adj_m, n_iter){

    degs <- Matrix::rowSums(adj_m)
    sapply(sort(unique(degs)), function(ll){
        if (ll == 1){
            .runBetaRan1(ll, n_iter)
        } else {
            .runBetaRan(ll, n_iter)
        }
    }, simplify=T)
}


getRandomBetasPar <- function(adj_m, n_iter, cores=1){

    degs <- Matrix::rowSums(adj_m)
    do.call(cbind,
            parallel::mclapply(sort(unique(degs)), function(ll){

        if (ll == 1){
            .runBetaRan1(ll, n_iter)
        } else {
            .runBetaRan(ll, n_iter)
        }

    }, mc.cores=cores))

}

getRealBetas <- function(bet_tab, adj_m, back_dis){

        rr_v <- rank(bet_tab[,"t value"])/nrow(bet_tab)
        degs <- Matrix::rowSums(adj_m)
        colnames(back_dis) <- sort(unique(degs))
        pbs <- apply(adj_m, 1, function(x){

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




















