#' Parameters of background distributions
#'
#' Parameters of Weibull distributions fitted on \eqn{rho} values computed
#' by RRA algorithm through random sampling of input ranks. This object covers
#' neighborhood sizes from 1 to 1999.
#'
#' @format A data frame with 1999 rows and 3 variables:
#' \describe{
#'   \item{alpha}{Alpha parameter of Weibull distribution}
#'   \item{beta}{Beta parameter of Weibull distribution}
#'   \item{degree}{Corresponding nodee degree (neighborhood size)}
#' }
"ww"
