# helper functions -------------------------------------------------------------

.prepareDataForFit <- function(x, yy){
    data.frame(
        yy, x
    )
}

.prepareFormula <- function(cond, covs, yy){
    stats::formula(
        paste0(
            "yy~",
            paste(
                c(cond,
                  covs),
                collapse="+"
            )
        )
    )
}

.fitModels <- function(formula, dat, yy){
    summary(stats::lm(
        formula,
        data=dat,
        method="qr"
    ))
}

.extractCoefs <- function(fits){
    lapply(fits, function(fit)
        fit[[4]][2,]
    )
}

.addStats <- function(x,y,cv){
    rowMeans(y[,which(x[,cv] == 0)])
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
fitLms <- function(x, y, cond, covs=NULL, filter=T){

    # Fit all models
    yy      <- t(y)
    dat     <- .prepareDataForFit(x,yy)
    form    <- .prepareFormula(cond, covs, yy)
    fits    <- .fitModels(form, dat, yy)
    coefs   <- .extractCoefs(fits)
    wtm     <- .addStats(x,y,cond)

    # Prepare output matrix
    mm <- do.call(rbind, coefs)

    mm <- cbind(
        mm,
        "wtm" = wtm
    )

    rownames(mm) <- colnames(yy)

    if (filter){

        mm <- mm[which(mm[,"wtm"] < 0),]

    }

    mm <- cbind(mm, "u"=rank(mm[,3])/nrow(mm))

    return(mm)

}



