# Class prodeInput =============================================================
#' \code{prodeInput} class
#'
#' @description prodeInput class extends \linkS4class{SummarizedExperiment} class
#' including adjacency matrix computed by \code{getProdeInput} function from
#' input \code{edge_list} and \code{modality}, a character string reporting
#' the type of score that needs to be computed by PRODE (\code{NIE_score} or
#' \code{NICE_score}).
#' @import methods
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
.prodeInput <- setClass(
    Class = "prodeInput",
    slots = representation(
        design     = "ANY",
        adjMatrix  = "Matrix",
        modality   = "character", 
        weights    = "ANY", 
        etab       = "ANY"
    ),
    contains="SummarizedExperiment"
)

# prodeInput Class Constructor.................................................

newProdeInput <- function(score_matrix, col_data, design, adjMatrix, modality, weights=NULL, etab=NULL){

    se <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(score_matrix = score_matrix),
        colData    = S4Vectors::DataFrame(col_data)
    )

    .prodeInput(
        se,
        design     = design,
        adjMatrix  = adjMatrix,
        modality   = modality, 
        weights    = weights, 
        etab       = etab
    )

}

# Class prodeResults ===========================================================

#' \code{prodeResults} class
#'
#' @description \code{prodeResults} class is an object which extends DataFrame
#' object (reporting results of PRODE run), by including \code{adjMatrix} as the
#' adjacency matrix resulting after PRODE run (in case some genes have been filtered out) and
#' \code{filterdData}, a DataFrame object containing scores and statistics of genes that have been filtered
#' out during the analyisis (depending on the analyses settings).
#' @slot results = DataFrame object, reporting the results of PRODE run.
#' @slot adjMatrix = Matrix object, an adjacency matrix as computed by \code{\link{getProdeInput}}.
#' @slot filteredData = a \code{\link[S4Vectors]{DataFrame}} object collecting genes that have been filtered out during the analysis.
#' @import methods
#' @export
.prodeResults <- setClass(
    Class = "prodeResults",
    slots = representation(
        results       = "DFrame",
        adjMatrix     = "Matrix",
        filteredData  = "DFrame"
    ),
    contains = "DFrame"
)

# prodeResults Class Constructor................................................

newProdeResults <- function(fit_tab, rra_tab, adjMatrix, modality, filtered){

    df <- S4Vectors::DataFrame(
        "gene" = rownames(fit_tab),
        fit_tab,
        rra_tab
    )

    colnames(df)[ncol(df)] <- modality

    .prodeResults(
        results = df,
        adjMatrix = adjMatrix,
        filteredData = S4Vectors::DataFrame(filtered)
    )

}




