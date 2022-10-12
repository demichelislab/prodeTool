# All classes ==================================================================

.prodeInput <- setClass(
    Class = "prodeInput",
    slots = representation(
        adjMatrix  = "Matrix",
        degLookup  = "integer"
    ),
    contains="SummarizedExperiment"
)

# Object Constructor ...........................................................

getProdeInput <- function(score_matrix, col_data, edge_table){

    # TODO: implement all checks on col_data

    .inputCheck(
        score_matrix = score_matrix,
        col_data     = col_data,
        edge_table   = edge_table
    )

    adj_m <- .getAdjMatr(
        etab = edge_table,
        gns  = rownames(score_matrix)
    )

    degLookup <- .getDegLookup(adj_m)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(score_matrix = score_matrix),
        colData    = S4Vectors::DataFrame(col_data)
    )

    .prodeInput(
        se,
        adjMatrix  = adj_m,
        degLookup  = degLookup
    )

}








