# Helper functions -------------------------------------------------------------

.vRunif <- Vectorize(function(n){
  runif(n)
}, "n")

# .vSample <- Vectorize(function(n, s){
#   sample(n, s)
# }, "s")

.colSort <- function(x){
  apply(x, 2, sort, method="quick")
}

# .runRhoRan1 <- function(ll, n_iter){
#     ran <- matrix(.vRunif(rep(ll, n_iter)), nrow=1)
#     pb  <- stats::pbeta(ran, 1:ll, ll - 1:ll + 1)
#     return(t(pb))
# }

.runRhoRan <- function(ll, n_iter){
    ran <- .colSort(.vRunif(rep(ll, n_iter)))
    pb  <- stats::pbeta(ran, 1:ll, ll - 1:ll + 1)
    mm  <- apply(pb, 2, min)
    return(mm)
}

.runRhoReal <- function(x, rr_v){
    re <- sort(rr_v[which(x!=0)], method="quick")
    ps <- stats::pbeta(re, seq_along(re), length(re) - seq_along(re) + 1)
    min(ps, na.rm=T)
}

.runRhoPval <- function(x, score, back_dis){ # empricial p-value
    mean(back_dis[,as.character(sum(x!=0))] <= score)
}

.runRhoPvalFitDistr <- function(x, score, back_par){ # empricial p-valeu
    params <- unlist(back_par[(sum(x!=0)-1),1:2])
    stats::pweibull(score, params[1], params[2])
}

.computeFinalScore <- function(u_gene, u_neig){

  log(u_gene*u_neig)

}

# Main functions ---------------------------------------------------------------

# Computes n_iter rhos with randomly sampled ranks from a unifrom distribution
# for each node degree computed from adjacency matrix
getRandomRhos <- function(adj_m, n_iter){

    degs <- Matrix::rowSums(adj_m)

    back_dis <- vapply(sort(unique(degs)), function(ll){
        .runRhoRan(ll, n_iter)
    }, rep(0.1, n_iter))

    colnames(back_dis) <- sort(unique(degs))
    return(back_dis)
}

# Same function as getRandomRhos() but with chance to run in parallel
getRandomRhosPar <- function(adj_m, n_iter, cores=1){

    degs <- Matrix::rowSums(adj_m)

    back_dis <- do.call(cbind,
      parallel::mclapply(sort(unique(degs)), function(ll){
        .runRhoRan(ll, n_iter)
      }, mc.cores=cores)
    )

    colnames(back_dis) <- sort(unique(degs))
    return(back_dis)

}

# Computes rhos given the input ranks, obtained from bet_tab. Works as
# getRealRhosFitDistr() except that p-values are computed according to
# background distributions obtained from getRandomRhos().
getRealRhos <- function(bet_tab, adj_m, back_dis, scaledEst){

    if (scaledEst){
      rr_v <- rank(bet_tab[,"t value"])/nrow(bet_tab)
    } else {
      rr_v <- rank(bet_tab[,"Estimate"])/nrow(bet_tab)
    }

    pbs <- apply(adj_m, 1, function(x){
        score <- .runRhoReal(x, rr_v)
        p.val <- .runRhoPval(x, score, back_dis)
        c(score, p.val)
    })

    u_n <- rank(pbs[2,], na.last = 'keep') /
      sum(!is.na(pbs[2,]))

    cbind(
        "rra_score"     = pbs[1,],
        "rra_p"         = pbs[2,],
        "rra_fdr"       = stats::p.adjust(pbs[2,], "fdr"),
        "u_gene"        = rr_v,
        "u_neig"        = u_n,
        'score'         = .computeFinalScore(rr_v, u_n)
    )
}

#' Computes \eqn{\rho} p-values according to RRA.
#'
#' @details For each gene's neighborhood, this function computes the rho value
#' as described by the RRA algorithm (Kolde et al., 2012) and a p-value
#' leveraging an input list of fitted background rhos distributions on randomly
#' shuffled ranks, depending on neighborhood size. The lower the p-value,
#' the more genes in neighborhood are distributed towards low ranks.
#'
#' @param bet_tab a table of linear models fit statistics as computed by \code{fitLms()}.
#' @param adj_m an adjacency matrix with rows and columns corresponding to rows of \code{beta_tab}.
#' @param back_par a table of parameters of background fitted distributions.
#' @param scaledEst logical, whether to used scaled coefficients or not to computee gene-level signal.
#' @returns matrix object with statistics of RRA for each gen in beta_tab.
getRealRhosFitDistr <- function(bet_tab, adj_m, back_par, scaledEst){

    if (scaledEst){
      rr_v <- rank(bet_tab[,"t value"], na.last = 'keep')/sum(!is.na(bet_tab[,"t value"]))
    } else {
      rr_v <- rank(bet_tab[,"Estimate"], na.last = 'keep')/sum(!is.na(bet_tab[,"t value"]))
    }

    pbs <- apply(adj_m, 1, function(x){

        score <- .runRhoReal(x, rr_v)
        p.val <- .runRhoPvalFitDistr(x, score, back_par)

        c(score, p.val)

    })

    #degs <- Matrix::rowSums(adj_m)
    #mr <- (adj_m %*% rr_v)[,1]/degs

    u_n <- rank(pbs[2,], na.last = 'keep') /
      sum(!is.na(pbs[2,]))

    cbind(
        #"mean_rank"     = mr,
        "rra_score"     = pbs[1,],
        "rra_p"         = pbs[2,],
        "rra_fdr"       = stats::p.adjust(pbs[2,], "fdr"),
        "u_gene"        = rr_v,
        "u_neig"        = u_n,
        'score'         = .computeFinalScore(rr_v, u_n)
    )
}
