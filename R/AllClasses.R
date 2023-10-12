# Class prode Input ============================================================

#' @import methods
#' @export
.prodeInput <- setClass(
    Class = "prodeInput",
    slots = representation(
        design     = "ANY",
        adjMatrix  = "Matrix",
        modality   = "character"
    ),
    contains="SummarizedExperiment"
)

# Object Constructor -----------------------------------------------------------

#' @import methods
#' @export
getProdeInput <- function(score_matrix, col_data, edge_table, design=NULL){

    ## Input check .............................................................

    .inputCheck(
        score_matrix = score_matrix,
        col_data     = col_data,
        edge_table   = edge_table,
        design       = design
    )

    if (is.null(design)){

        message(paste0('NOTE: No formula has been given - ',
                       'Preparing Input for NIE scores computation',
                       ' (essentiality analysis)'))

        mod <- 'NIE_score'

        design <- stats::model.matrix.default(
            as.formula("~1"),
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

    se <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(score_matrix = score_matrix),
        colData    = S4Vectors::DataFrame(col_data)
    )

    .prodeInput(
        se,
        design     = design,
        adjMatrix  = adj_m,
        modality   = mod
    )

}

# Class prode results ==========================================================

#' @import methods
#' @export
prodeResults <- setClass(
    Class = "prodeResults",
    slots = representation(
        adjMatrix     = "Matrix",
        filteredData  = "DFrame"
    ),
    contains = "DFrame"
)

# Object Constructor -----------------------------------------------------------

#' @import methods
#' @export
prodeResults <- function(output_df, adj_m, filteredData){

    new("prodeResults", output_df, adjMatrix=adj_m, filteredData=filteredData)

}






