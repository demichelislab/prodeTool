# helper functions -------------------------------------------------------------

.fitModels <- function(yy, x){
    summary(stats::lm(
        yy~x
    ))
}

.extractCoefs <- function(fits){
    lapply(fits, function(fit)
        fit[[4]][nrow(fit[[4]]),]
    )
}

.addStats <- function(y,x){
    rowMeans(y[,which(x[,ncol(x)] == 0)])
}

# Main function ----------------------------------------------------------------

#' Runs linear models on given as input a prode object.
#'
#' @param x
#' @param y
#' @param covs
#'
#' @return matrix with summary statistics of each fit
#' @export
#'
#' @examples
fitLms <- function(x, y){

    # Fit all models
    yy      <- t(y)
    fits    <- .fitModels(yy,x)
    coefs   <- .extractCoefs(fits)
    wtm     <- .addStats(y,x)

    # Prepare output matrix
    mm <- do.call(rbind, coefs)

    mm <- cbind(
        mm,
        "wtm" = wtm
    )

    rownames(mm) <- colnames(yy)

    return(mm)

}



