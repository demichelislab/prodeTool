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



