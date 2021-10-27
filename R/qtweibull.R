#' Quantile Function for Truncated Weibull
#'
#' This is an attempt to provide a more numerically stable
#' version of the truncated Weibull distribution, especially
#' for more extreme truncation points. Underneith the hood it
#' relies on \code{uniroot}, which may issue errors. This should
#' be considered the least stable of the functions associated
#' with this distribution in this package.
#'
#' @param p (Numeric) Vector of probabilities
#' @param a (Numeric) Vector of lower truncation points
#' @param b (Numeric) Vector of upper truncation points
#' @param shape (Numeric) Vector of Weibull shape parameters
#' @param scale (Numeric) Vector of Weibull scale parameters
#' @param log.p (Logical) Whether p represents p or \eqn{log(p)}
#' @param ...
#'
#' @return Numeric vector of quantile function values
#' @export
#'
#' @examples
#' qtweibull(0.5, 1, 5, 1.5, 2.5)
qtweibull = function(p, a=0, b=Inf, shape, scale, log.p=FALSE, ...){
  log_p = if(isFALSE(log.p)) log(p) else p
  res = a^shape - scale^shape * log1mpexp(log_p, a, b, shape, scale)

  res^(1/shape)
}
