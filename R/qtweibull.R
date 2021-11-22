#' Quantile Function for Truncated Weibull
#'
#' This is an attempt to provide a more numerically stable
#' version of the truncated Weibull distribution, especially
#' for more extreme truncation points. Underneith the hood it
#' relies on \code{uniroot}, which may issue errors. For this reason \code{qtweibull}
#' should be considered the least stable component of this package.
#' It's also about 50 times slower than \code{stats::qweibull} in cases
#' where they both apply.
#'
#' @param p (Numeric) Vector of probabilities
#' @param a (Numeric) Vector of lower truncation points
#' @param b (Numeric) Vector of upper truncation points
#' @param shape (Numeric) Vector of Weibull shape parameters (default: 0)
#' @param scale (Numeric) Vector of Weibull scale parameters (default: Inf)
#' @param log.p (Logical) Whether p represents p or \eqn{log(p)} (default: FALSE)
#' @param ...
#'
#' @return (Numeric) vector of quantile function values
#' @importFrom VGAM log1mexp
#' @export
#'
#' @examples
#' qtweibull(0.5, 1, 5, 1.5, 2.5)
qtweibull = function(p, shape, scale, a=0, b=Inf, log.p=FALSE, ...){
  stopifnot(
    exprs = {
      is.numeric(p); is.numeric(a); is.numeric(b);
      is.numeric(shape); is.numeric(scale); all(a >= 0);
      all(b > a); all(shape > 0); all(scale > 0);
      is.logical(log.p);
      all(is.finite(shape)); all(is.finite(scale));
      all(is.finite(a));
      length(p) >= 1L; length(shape) >= 1L; length(scale) >= 1L;
      length(a) >= 1L; length(b) >= 1L; length(log.p) == 1L
    }
  )

  if(isTRUE(log.p)){
    stopifnot(exprs = {all(p <= 0)})
  } else {
    stopifnot(exprs = {all(p >= 0); all(p <= 1)})
  }

  vec_args = check_args(
    log_p = if(isFALSE(log.p)) log(p) else p,
    a = a,
    b = b,
    shape = shape,
    scale = scale,
    expand = TRUE
  )

  log1mpexp_eval = do.call("mapply", c(list(FUN = function(log_p, a, b, shape, scale) log1mpexp(log_p, a, b, shape, scale)), vec_args))

  res = a^shape - scale^shape * log1mpexp_eval

  res^(1/shape)
}
