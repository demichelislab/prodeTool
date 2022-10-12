#' @import methods
#' @include class.R
NULL

#' @rdname adjMatrix
#' @usage NULL
#' @export
setGeneric("adjMatrix", function(object) standardGeneric("adjMatrix"))

#' @rdname adjMatrix
#' @aliases adjMatrix
#' @usage NULL
setMethod("adjMatrix", signature = "prodeInput",
          definition = function(object) object@adjMatrix)

#' @rdname NodeDegreeLookup
#' @usage NULL
#' @export
setGeneric("NodeDegreeLookup", function(object) standardGeneric("NodeDegreeLookup"))

#' @rdname NodeDegreeLookup
#' @aliases adjMatrix
#' @usage NULL
setMethod("NodeDegreeLookup", signature = "prodeInput",
          definition = function(object) object@degLookup)
