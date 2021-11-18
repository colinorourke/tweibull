#' Mass function of truncated Weibull distribution
#'
#' @param x Numeric, value at which to evaluate mass function
#' @param shape Numeric, Weibull shape parameter
#' @param scale Numeric, Weibull scale parameter
#' @param log Logical, whether to return log mass function value
#' @param a Numeric, lower truncation point
#' @param b Numeric, upper truncation point
#'
#' @return Numeric value
#' @export
#' @importFrom VGAM log1mexp
#' @import checkmate
#'
#' @examples
#' dtweibull(3, 0.5, 1.5, 2)
dtweibull = function(x, shape, scale=1, a=0, b=Inf, log=FALSE){
  assertNumeric(a, lower=0)
  assertNumeric(b, lower=a)
  assertNumeric(shape, lower=0)
  assertNumeric(scale, lower=0)

  args = check_args(
    x = x,
    shape = shape,
    scale = scale,
    a = a,
    b = b
  )

  lns = vapply(arg_list, length, 1L)

  ld = function(x, shape, scale, a, b){
    log(shape) - shape * log(scale) + log(x) * (shape - 1) - (a/scale)^shape - ldenom(a, b, shape, scale)
  }

  low = rep_len(x < a, lns)
  hi = rep_len(x > b, lns)
  mid = rep_len(!low & !hi, lns)

  log_d = numeric(max(lns))
  if(any(low)) log_d[low] = -Inf
  if(any(hi)) log_d[hi] = 0
  if(any(mid)){
    log_d[mid] = do.call("ld", list_select(x = args, ind = mid))
  }

  if(isTRUE(log)) log_d else exp(log_d)
}
