# Utils for adjacency matrix ---------------------------------------------------

.checkAdjMatr <- function(all_gns_adj, all_gns_sco){

    if (!all(all_gns_sco%in%all_gns_adj)){
        stop("Not all genes from score matrix are in edge table!")
    }

}

.getAdjMatr <- function(etab, gns){

    # gns2n <- setNames(1:length(gns), gns)
    #
    # idxm <- cbind(
    #     c(as.numeric(gns2n[etab[,1]]), as.numeric(gns2n[etab[,2]]),1:length(gns)),
    #     c(as.numeric(gns2n[etab[,2]]), as.numeric(gns2n[etab[,1]]),1:length(gns))
    # )

    idxm <- cbind(
        c(    match(etab[,1], gns),     match(etab[,2], gns), seq_along(gns)),
        c(    match(etab[,2], gns),     match(etab[,1], gns), seq_along(gns))
    )

    idxm <- idxm[stats::complete.cases(idxm), ] # removes those score matrix genes not present  in edge list

    mm <- Matrix::Matrix(0, length(gns), length(gns))
    colnames(mm) <- rownames(mm) <- gns
    mm[idxm] <- 1
    return(mm)

}

.filter0DegAdj <- function(adj_mat, deg){
    adj_mat[-which(deg == 1), -which(deg == 1)]
}

.filter0DegBet <- function(beta_tab, deg){
    beta_tab[-which(deg == 1),]
}

.filterAdjMatrix <- function(filterCtrl, prodeInput, fit_tab){

    # 1. Filtering based on mean of wild-type models ...........................
    filtered <- data.frame()

    if (filterCtrl){
        keep     <- which(
            .getCtrlMean(
                SummarizedExperiment::assay(prodeInput),
                designMatrix(prodeInput)
            ) < 0
        )
        filtered <- S4Vectors::bindROWS(filtered, list(fit_tab[-keep,]))
        fit_tab  <- fit_tab[ keep,]
    }

    adj_m    <- subsetAdjMatrix(prodeInput, rownames(fit_tab)) # update adj matrix

    # 2. Filtering based on degree == 1 (meaning no connections) ...............
    dds      <- Matrix::rowSums(adj_m)

    if (any(dds == 1)){

        if (length(which(dds == 1)) == 1){
            tmp <- t(fit_tab[which(dds == 1),])
            rownames(tmp) <- rownames(fit_tab)[which(dds == 1)]
            filtered <- S4Vectors::bindROWS(filtered, list(tmp))
        } else {
            filtered <- S4Vectors::bindROWS(filtered, list(fit_tab[which(dds == 1),]))
        }

        adj_m    <- .filter0DegAdj(adj_m, dds)
        fit_tab  <- .filter0DegBet(fit_tab, dds)

    }

    if (nrow(adj_m) < 1){
        stop('No genes remained after filtering. ',
             'Please check if input edge table has ',
             'same genes as row-names of score matrix')
    }

    if (nrow(adj_m)/nrow(adjMatrix(prodeInput)) <= 0.5){
        warning('More than 50% of the genes have ',
                'been filtered out as not present in edge table!')
    }

    list(
        filtered  = S4Vectors::DataFrame(filtered),
        fit_tab   = fit_tab,
        adjMatrix = adj_m
    )

}

.checkBetAdj <- function(beta_tab, adj_mat, deg){
    all(rownames(beta_tab) == rownames(adj_mat)) &
    all(rownames(beta_tab) == colnames(adj_mat))
}

.findCommonIndexes <- function(ci_splits, weights, etab, wdir='<'){
  
  thrs <- seq(min(enet$score), max(enet$score), length.out = ci_splits)[
  1:(ci_splits -1)]
  
  ids <- lapply(1:length(thrs), function(i){
    which(weights >= thrs[i])
  })
  
  all_gns <- Reduce(intersect, lapply(ids, function(kk){
    etab_loc <- etab[kk,-3]
    unique(unlist(etab_loc))
  }))
  
  etab <- etab[which(etab[,1] %in% all_gns & etab[,2] %in% all_gns),]
  
  list(
    etab, ids
  )
  
}

.getCI <- function(points){
  # Compute the mean and standard error
  mean_value <- mean(points)
  stderr <- sd(points) / sqrt(length(points))
  
  # Compute the 95% CI (using t-distribution)
  alpha <- 0.05  # 95% CI -> alpha = 0.05
  t_critical <- qt(1 - alpha / 2, df = length(points) - 1)  # t-critical value
  
  # Lower and upper bounds of the CI
  ci_lower <- mean_value - t_critical * stderr
  ci_upper <- mean_value + t_critical * stderr
  
  c(ci_lower, ci_upper)
}

.mergeResTable <- function(res_tab, gns){
  res_tab[match(gns, res_tab$gene),] -> res_tab
  return(res_tab)
}

# Checks for input data --------------------------------------------------------

.inputCheck <- function(
    score_matrix, col_data, design=NULL, edge_table, edge_weights=NULL
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

    # Input design -------------------------------------------------------------

    if (!is.null(design)){

        if (!inherits(design,"formula")){
            stop("Input design should be a formula, please pass it to as.formula()")
        }

        if (diff(dim(attr(stats::terms(design), "factors"))) != 0){
            stop(
            paste0("Please check your formula, no left-side term is needed.\n",
                 "Formula should be of the form: '~covariates'."))
        }

        if (!all(colnames(attr(stats::terms(design), "factors")) %in% colnames(col_data))){
            stop("Some variables in design formula are not in column data")
        }

        if (!all(colnames(attr(stats::terms(design), "factors")) %in% colnames(col_data))){
            stop("Some variables in design formula are not in column data")
        }

        ff <- colnames(attr(stats::terms(design), "factors"))[ncol(attr(stats::terms(design), "factors"))]
        if (!length(unique(col_data[,ff])) == 2){
            stop("Grouping variable does not have two groups")
        }

        if (any(is.na(unique(col_data[,ff])))){
            stop("Some of the values in grouping variable are NAs")
        }

        if (! all(col_data[,ff] %in% c(1,0))){
            stop("Grouping variable is not 0,1 encoded, please modify it")
        }

    }

    # Input edge table ---------------------------------------------------------

    if (dim(edge_table)[2] != 2){
        stop("Edge table has more than two columns.")
    }
  
    # Input edge weigths 
  
    if (!is.null(edge_weights)){

      if (!is.numeric(edge_weights)){
        stop('Provided edge weights are not numeric')
      }

      if (length(edge_weights)!=nrow(edge_table)){
        stop('Number of weights and rows of edge table do not match')
      }
      
      if (any(is.na(edge_weights))){
        stop('Some weights contain NAs. Please remove those edges and weights.')
      }

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
