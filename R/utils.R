# ------------------------------------------------------------------------------

# helper functions

.getDegLookup <- function(adj_mat){
    rowSums(adj_mat > 0)
}

.checkAdjMatr <- function(all_gns_adj, all_gns_sco){

    if (!all(all_gns_sco%in%all_gns_adj)){
        stop("Not all genes from score matrix are in edge table!")
    }

}

# .filterMatr <- function(mm, gns){
#
#     if (!nrow(mm)==ncol(mm)){
#         mis <- gns[which(!gns%in%colnames(mm))]
#         mm  <- cbind(mm, rep(0, nrow(mm)))
#         colnames(mm)[ncol(mm)] <- mis
#     }
#
#     mm[gns,gns]
#
# }
#
# .getAdjMatr <- function(etab, gns){
#
#     mm <- Matrix((unclass(table(etab)) > 0)*1)
#
#     .checkAdjMatr(
#         rownames(mm),
#         gns
#     )
#
#     mm <- .filterMatr(
#         mm,
#         gns
#     )
#
#     return(mm)
# }

.getAdjMatr <- function(etab, gns){

    gns2n <- 1:length(gns)
    names(gns2n) <- gns

    idxm <- cbind(
        c(as.numeric(gns2n[etab[,1]]), as.numeric(gns2n[etab[,2]]), 1:length(gns)),
        c(as.numeric(gns2n[etab[,2]]), as.numeric(gns2n[etab[,1]]), length(gns):1)
    )

    idxm <- idxm[complete.cases(idxm), ]

    mm <- Matrix(0, length(gns), length(gns))
    colnames(mm) <- rownames(mm) <- gns
    mm[idxm] <- 1
    return(mm)
}


.filter0DegAdj <- function(adj_mat, deg){ # this will be a method of object
    adj_mat[-which(deg == 0), -which(deg == 0)]
}

.filter0DegBet <- function(beta_tab, deg){
    beta_tab[-which(deg == 0),]
}

.filterDeg <- function(deg){
    deg[-which(deg == 0)]
}

.checkBetAdjDeg <- function(beta_tab, adj_mat, deg){
    all(rownames(beta_tab) == rownames(adj_mat)) &
    all(rownames(beta_tab) == colnames(adj_mat)) &
    all(rownames(beta_tab) == names(deg))
}

.vRunif <- Vectorize(function(n){
    runif(n)
}, "n")

.vSample <- Vectorize(function(n, s){
    sample(n, s)
}, "s")

.colSort <- function(x){
    #apply(x, 2, sort, method="quick")
    Rfast::colSort(x)
}


.inputCheck <- function(
    score_matrix, col_data, edge_table
){

    if (!is.matrix(score_matrix)){
        stop("Input score table is not a matrix object.")
    }
    if (any(is.na(score_matrix))){
        stop("Input score matrix contains NA values.")
    }
    if (any(is.infinite(score_matrix))){
        stop("Infinite values in input score matrix are not allowed.")
    }
    if (!is.numeric(score_matrix)){
        stop("Input score matrix is not numeric.")
    }
    if (is.null(colnames(score_matrix))){
        stop("Input score matrix must contain unique samples identifiers as colnames.")
    }
    if (is.null(rownames(score_matrix))){
        stop("Input matrix must contain unique genes identifiers as rownames.")
    }

    # Input col_data -------------------------------------------------------

    if (ncol(score_matrix) != nrow(col_data)){
        stop("Input scores matrix and colDat have different numbers of samples.")
    }
    if (!all(rownames(col_data) %in% colnames(score_matrix))){
        stop("Sample info file has different sample ids than scores matrix.")
    }

    # # Input condition ----------------------------------------------------------
    #
    # if (!condition%in%colnames(samples_info)){
    #     stop("Condition must be a column name of column_data.")
    # }
    # if (any(is.na(samples_info[[condition]]))){
    #     stop("Condition has NA values.")
    # }
    # if (!length(unique(samples_info[[condition]]))==2){
    #     stop("Inspected condition has more than 2 levels.")
    # }
    #
    # # Input covariates ---------------------------------------------------------
    #
    # if (!all(covariates%in%colnames(samples_info))){
    #     stop("Some covariates are not present in samples_info.")
    # }
    #
    # if (any(vapply(covariates, function(cc) length(unique(samples_info[[cc]]))<2, T))){
    #     stop("Some covariates have less than two levels.")
    # }
    #
    # Input edge table ---------------------------------------------------------

    if (dim(edge_table)[2] != 2){
        stop("Edge table has more than two columns.")
    }

    # all_gns <- unique(c(edge_table))
    #
    # if (!all(rownames(scores_matrix)%in%all_gns)){
    #     stop("Not all gene names of score matrix are in edge_table.")
    # }
    # if (!all(all_gns%in%rownames(scores_matrix))){
    #     stop("Some genes in edge_table are not in score matrix.")
    # }

}
