#' Split a string
#'
#' @param x A character vector with one element.
#' @param split  What to split on.
#'
#' @return A character vector
#' @export
#'
#' @examples
#' x <- "a,b,c,d,e"
#' strsplit1(x, split = ",")
strsplit1 <- function(x, split) {
    strsplit(x, split = split)[[1]]
}

#' Split a string
#'
#' @param string A character vector with, at most, one element.
#' @inheritParams stringr::str_split
#'
#' @return A character vector.
#' @export
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' str_split_one(x, pattern = ",")
#' str_split_one(x, pattern = ",", n = 2)
#'
#' y <- "192.168.0.1"
#' str_split_one(y, pattern = stringr::fixed("."))
str_split_one <- function(string, pattern, n = Inf) {
    stopifnot(is.character(string), length(string) <= 1)
    if (length(string) == 1) {
        stringr::str_split(string = string, pattern = pattern, n = n)[[1]]
    } else {
        character()
    }
}
