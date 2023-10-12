#' @import methods
#' @include AllClasses.R
NULL

# Accessors for ProdeInput =====================================================

#' @rdname adjMatrix
#' @usage NULL
#' @export
setGeneric("adjMatrix", function(object) standardGeneric("adjMatrix"))

#' @rdname adjMatrix
#' @aliases adjMatrix
#' @usage NULL
setMethod("adjMatrix", signature = "prodeInput",
          definition = function(object) object@adjMatrix)

#' @rdname designMatrix
#' @usage NULL
#' @export
setGeneric("designMatrix", function(object) standardGeneric("designMatrix"))

#' @rdname designMatrix
#' @aliases designMatrix
#' @usage NULL
setMethod("designMatrix", signature = "prodeInput",
          definition = function(object) object@design)

#' @rdname modality
#' @usage NULL
#' @export
setGeneric("modality", function(object) standardGeneric("modality"))

#' @rdname modality
#' @aliases modality
#' @usage NULL
setMethod("modality", signature = "prodeInput",
          definition = function(object) object@modality)

#' @rdname filterAdjMatrix
#' @usage NULL
#' @export
setGeneric("filterAdjMatrix", function(object, ...) standardGeneric("filterAdjMatrix"))

#' @rdname filterAdjMatrix
#' @aliases filterAdjMatrix
#' @usage NULL
setMethod("filterAdjMatrix", signature = "prodeInput",
          definition = function(object, nms) object@adjMatrix[nms,nms])

# Accessors for ProdeResults ===================================================

#' @rdname results
#' @usage NULL
#' @export
setGeneric("results", function(object) standardGeneric("results"))

#' @rdname results
#' @aliases filterNodesDegree
#' @usage NULL
setMethod("results", signature = "prodeResults",
          definition = function(object) as.data.frame(object))

#' @rdname adjMatrixResults
#' @usage NULL
#' @export
setGeneric("adjMatrixResults", function(object) standardGeneric("adjMatrixResults"))

#' @rdname adjMatrixResults
#' @aliases adjMatrixResults
#' @usage NULL
setMethod("adjMatrixResults", signature = "prodeResults",
          definition = function(object) object@adjMatrix)

#' @rdname filteredGenes
#' @usage NULL
#' @export
setGeneric("filteredGenes", function(object) standardGeneric("filteredGenes"))

#' @rdname filteredGenes
#' @aliases filterNodesDegree
#' @usage NULL
setMethod("filteredGenes", signature = "prodeResults",
          definition = function(object) object@filteredData)
