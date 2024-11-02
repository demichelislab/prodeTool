#' Constructs input object for PRODE run.
#'
#' @details This function computes an adjacency matrix from edge list and filters
#'    out all genes that are not present in the edge-list or score matrix. The resulting
#'    \linkS4class{prodeInput} object inherits from \linkS4class{SummarizedExperiment}.
#' @param design this is the design formula required for NICE score computation. In the case of
#'    NIE score computation, this can be ignored (default set to \code{NULL}). \code{design} has to
#'    be a \code{formula} object. Last variable in formula will be considered as the group
#'    variable for the linear model fit, requiring 1 for case group and 0 for control group.
#' @param score_matrix matrix of input scores with genes on rows
#'    and samples on columns. No \code{NA} values are allowed. Gene names are required
#'    as \code{rownames}. Sample names are required as \code{colnames} and need to be
#'    identical to \code{col_data} column names.
#' @param col_data data.frame of column data, as described by \linkS4class{SummarizedExperiment}. Row names
#'    need to be identical to score-matrix column names and correspond to samples.
#'    Every sample-level included in design formula (if provided) should be included in \code{col_data}.
#'    The group variable, encoding the sample groups that need to be compared for NICE score
#'    computation, has to be a binary integer vector, encoding 1 for the case-group and 0 for the
#'    control-group.
#' @param edge_table matrix of interactions between gene pairs. This matrix displays
#'    two columns for each gene in the pair. Genes not present in \code{score_matrix} row names
#'    will be discarded.
#'
#' @returns An object of class \linkS4class{prodeInput}.
#'
#' @export
getProdeInput <- function(score_matrix, col_data, edge_table, design=NULL, 
                          edge_weights=NULL){

    ## Input check .............................................................

    .inputCheck(
        score_matrix = score_matrix,
        col_data     = col_data,
        edge_table   = edge_table,
        design       = design, 
        edge_weights = edge_weights
    )

    if (is.null(design)){

        message(paste0('NOTE: No formula has been given - ',
                       'Preparing Input for NIE scores computation',
                       ' (essentiality analysis)'))

        mod <- 'NIE_score'

        design <- stats::model.matrix.default(
            stats::as.formula("~1"),
            as.data.frame(col_data)
        )

    } else {

        message(paste0('NOTE: Formula has been given - ',
                       'Preparing Input for NICE scores computation',
                       ' (context-essentiality analysis)'))

        mod <- 'NICE_score'

        design <- stats::model.matrix.default(
            design,
            as.data.frame(col_data)
        )

    }

    ## Get adjMatr .............................................................

    adj_m <- .getAdjMatr(
        etab = edge_table,
        gns  = rownames(score_matrix)
    )

    stopifnot(.checkBetAdj(score_matrix, adj_m))

    ## Class Construction ......................................................

    if (!is.null(edge_weights)){
       message(paste0('Set up for confidence intervals across Edges Weights...'))
    }
    
    newProdeInput(
        score_matrix = score_matrix,
        col_data     = col_data,
        design       = design,
        adjMatrix    = adj_m,
        modality     = mod, 
        weights      = edge_weights, 
        etab         = edge_table
    )

}












