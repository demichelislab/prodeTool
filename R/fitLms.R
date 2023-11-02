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

.extractCoefsEss <- function(fits){
  lapply(fits, function(fit)
    fit[[4]][1,]
  )
}

.getCtrlMean <- function(y,x){

  rowMeans(y[,which(x[,ncol(x)] == 0)], na.rm=T)

}

# Computes extended statistics for each group in NICE scores computation.
.computeStats <- function(y, x, extendedNICEStats){

    ctrl_mean <- rowMeans(y[,which(x[,ncol(x)] == 0)], na.rm=T)
    case_mean <- rowMeans(y[,which(x[,ncol(x)] == 1)], na.rm=T)
    ctrl_sd   <- apply(y[,which(x[,ncol(x)] == 0)], 1,  function(x) sd(x, na.rm=T))
    case_sd   <- apply(y[,which(x[,ncol(x)] == 1)], 1, function(x) sd(x, na.rm=T))
    ctrl_n    <- sum(x[,ncol(x)] == 0, na.rm=T)
    case_n    <- sum(x[,ncol(x)] == 1, na.rm=T)

    cbind(
        "ctrl_mean" = ctrl_mean,
        "case_mean" = case_mean,
        "ctrl_sd"   = ctrl_sd,
        "case_sd"   = case_sd,
        "ctrl_n"    = ctrl_n,
        "case_n"    = case_n
    )

}

# Main function ----------------------------------------------------------------

#' Fit linear models on each row of score-matrix
#'
#' @param x a model.matrix object containing variables included in the linear model.
#'     This object is expected to have rows corresponding to the score-matrix columns.
#' @param y a score-matrix with a number of columns equal to x, no NA values are allowed.
#' @param extendedNICEStats a logical if extended stats need to be computed.
#' @returns matrix object with statistics of linear model fits.
#' @examples
#' \dontrun{
#'
#'     y <- matrix(rnorm(100), 10, 10)
#'
#'     col_data <- data.frame('gr' = rep.int(c(1,0), c(5,5)))
#'
#'     x <- stats::model.matrix.default(
#'           stats::as.formula("~gr"),
#'           col_data
#'     )
#'
#'     fitLms(x, y)
#' }
#'
fitLms <- function(x, y, extendedNICEStats=F){

    if (any(is.na(y))){
        stop('Matrix contains NAs, this would cause',
             ' altered linear model estimates.')
    }

    if (any(is.na(x))){
        warning('There are some NAs in the design matrix: ',
                'rows will be excluded from model fit.')
    }

    # Fit all models
    yy      <- t(y)
    fits    <- .fitModels(yy,x)
    coefs   <- .extractCoefs(fits)

    # Prepare output matrix
    mm <- do.call(rbind, coefs)

    if (extendedNICEStats){ # If NICE score is computed

      if (all(x[,ncol(x)] != 0)){

          warning('Computing extended stats for NIE scores will ',
                  'not produce any estimate for control group.')

      }

      info    <- .computeStats(y,x)
      mm <- cbind(
        mm,
        info
      )

    }

    rownames(mm) <- colnames(yy)

    return(mm)

}



