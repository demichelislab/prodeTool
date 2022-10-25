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

.computeStats <- function(y, x, extendedStats){

    ctrl_mean <- rowMeans(y[,which(x[,ncol(x)] == 0)], na.rm=T)

    if (extendedStats) {
        case_mean <- rowMeans(y[,which(x[,ncol(x)] == 1)], na.rm=T)
        ctrl_sd   <- apply(y[,which(x[,ncol(x)] == 0)], 1,  function(x) sd(x, na.rm=T))
        case_sd   <- apply(y[,which(x[,ncol(x)] == 1)], 1, function(x) sd(x, na.rm=T))
        ctrl_n    <- sum(x[,ncol(x)] == 0)
        case_n    <- sum(x[,ncol(x)] == 1)

        cbind(
            "ctrl_mean" = ctrl_mean,
            "case_mean" = case_mean,
            "ctrl_sd"   = ctrl_sd,
            "case_sd"   = case_sd,
            "ctrl_n"    = ctrl_n,
            "case_n"    = case_n
        )
    } else {

        mm <- matrix(ctrl_mean, ncol=1)
        colnames(mm) <- "ctrl_mean"
        mm
    }

}

# Main function ----------------------------------------------------------------

#' Fit linear models on input data
#'
#' @param x a model.matrix object containing variables included in the linear model.
#'     This object is expected to have rows corresponding to the score-matrix columns.
#' @param y a score-matrix with a number of columns equal to x.
#' @param extendedStats a logical if extended stats need to be computed.
#'
#' @return a matrix object with results of linear model fitting procedure.
#' @export
fitLms <- function(x, y, extendedStats=F){

    # Fit all models
    yy      <- t(y)
    fits    <- .fitModels(yy,x)
    coefs   <- .extractCoefs(fits)
    info    <- .computeStats(y,x,extendedStats)

    # Prepare output matrix
    mm <- do.call(rbind, coefs, extendedStats)

    mm <- cbind(
        mm,
        info
    )

    rownames(mm) <- colnames(yy)

    return(mm)

}



