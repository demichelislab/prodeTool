#' @import methods
#' @include allClasses.R
NULL

# Accessors for ProdeInput =====================================================

#' adjMatrix
#'
#' @description This function allows to access the adjacency matrix in a \linkS4class{prodeInput}
#' object.
#' @rdname adjMatrix
#' @usage NULL
#' @export
#' @returns and object of class Matrix.
setGeneric("adjMatrix", function(object) standardGeneric("adjMatrix"))

#' @rdname adjMatrix
#' @aliases adjMatrix
#' @usage NULL
setMethod("adjMatrix", signature = "prodeInput",
          definition = function(object) object@adjMatrix)

#' designMatrix
#'
#' @description This function allows to access the design matrix in a \linkS4class{prodeInput}
#' object.
#' @rdname designMatrix
#' @usage NULL
#' @export
#' @returns and object of class Matrix.
setGeneric("designMatrix", function(object) standardGeneric("designMatrix"))

#' @rdname designMatrix
#' @aliases designMatrix
#' @usage NULL
setMethod("designMatrix", signature = "prodeInput",
          definition = function(object) object@design)

#' modality
#'
#' @description This function allows to access modality in a \linkS4class{prodeInput}
#' object.
#' @returns a character string.
#' @rdname modality
#' @usage NULL
#' @export
setGeneric("modality", function(object) standardGeneric("modality"))

#' @rdname modality
#' @aliases modality
#' @usage NULL
setMethod("modality", signature = "prodeInput",
          definition = function(object) object@modality)

#' @rdname subsetAdjMatrix
#' @usage NULL
setGeneric("subsetAdjMatrix", function(object, ...) standardGeneric("subsetAdjMatrix"))

#' @rdname subsetAdjMatrix
#' @aliases subsetAdjMatrix
#' @usage NULL
setMethod("subsetAdjMatrix", signature = "prodeInput",
          definition = function(object, nms) object@adjMatrix[nms,nms])

# Accessors for ProdeResults ===================================================

#' results
#'
#' @description This function allows to access results in a \linkS4class{prodeResults}
#' object. Details on columns are reported in \code{\link{runProde}}
#' @returns a  DataFrame object.
#' @rdname results
#' @usage NULL
#' @export
setGeneric("results", function(object) standardGeneric("results"))

#' @rdname results
#' @aliases results
#' @usage NULL
setMethod("results", signature = "prodeResults",
          definition = function(object) as.data.frame(object))

#' adjMatrixResults
#'
#' @description This function allows to access the adjacency matrix in a \linkS4class{prodeResults}
#' object.
#' @rdname adjMatrixResults
#' @usage NULL
#' @export
#' @returns and object of class Matrix.
setGeneric("adjMatrixResults", function(object) standardGeneric("adjMatrixResults"))

#' @rdname adjMatrixResults
#' @aliases adjMatrixResults
#' @usage NULL
setMethod("adjMatrixResults", signature = "prodeResults",
          definition = function(object) object@adjMatrix)

#' filteredGenes
#'
#' @description This function allows to access the genes filtered out in a \linkS4class{prodeResults}
#' object.
#' @returns a DataFrame object.
#' @rdname filteredGenes
#' @usage NULL
#' @export
setGeneric("filteredGenes", function(object) standardGeneric("filteredGenes"))

#' @rdname filteredGenes
#' @aliases filteredGenes
#' @usage NULL
setMethod("filteredGenes", signature = "prodeResults",
          definition = function(object) object@filteredData)
