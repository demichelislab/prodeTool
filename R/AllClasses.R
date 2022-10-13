# Class prode Input ============================================================

#' @import methods
#' @export
.prodeInput <- setClass(
    Class = "prodeInput",
    slots = representation(
        design     = "ANY",
        adjMatrix  = "Matrix"
    ),
    contains="SummarizedExperiment"
)

# Object Constructor -----------------------------------------------------------

#' @import methods
#' @export
getProdeInput <- function(score_matrix, col_data, design, edge_table){

    # TODO: implement all checks on col_data
    # TODO: remove the adjacency matrix

    ## Input check .............................................................

    .inputCheck(
        score_matrix = score_matrix,
        col_data     = col_data,
        edge_table   = edge_table
    )

    design <- stats::model.matrix.default(design, as.data.frame(col_data))

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
        adjMatrix  = adj_m
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






