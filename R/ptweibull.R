#' Cumulative Distribution Function for Truncated Weibull
#'
#' @param q Numeric, value at which to evaluate CDF
#' @param shape Numeric, Weibull shape parameter
#' @param scale Numeric, Weibull scale parameter
#' @param a Numeric, lower truncation point
#' @param b Numeric, upper truncation point
#' @param lower.tail Logical, whether to compute lower or upper tail probabilities
#' @param log.p Logical, whether to return log probabilities
#'
#' @return Numeric value(s)
#' @import checkmate
#' @export
#'
#' @examples
#' ptweibull(2, shape=1.5, a=1)
ptweibull = function(q, shape, scale=1, a=0, b = Inf, lower.tail=TRUE, log.p=FALSE){
  assertNumeric(a, lower=0, finite = TRUE)
  assertNumeric(b, lower=a)
  assertNumeric(shape, lower=0)
  assertNumeric(scale, lower=0)
  assertNumeric(q)
  assertLogical(lower.tail)

  arg_list = check_args(
    q = q,
    shape = shape,
    scale = scale,
    a = a,
    b = b,
    lower.tail = lower.tail
  )

  lns = vapply(arg_list, length, 1L)

  pfun = function(q, shape, scale, a, b, lower.tail){
    lt_ind = as.numeric(lower.tail)

    lt_ind * ldenom(a, q, shape, scale) + (1 - lt_ind) * ldenom(q, b, shape, scale) - ldenom(a, b, shape, scale)
  }

  low = rep_len(q <= a, lns)
  hi = rep_len(q >= b, lns)
  mid = rep_len(!low & !hi, lns)

  logp = numeric(max(lns))
  logp[low] = ifelse(lower.tail, -Inf, 0)
  logp[hi] = ifelse(lower.tail, 0, -Inf)
  logp[mid] = do.call("pfun", list_select(x = arg_list, ind = mid))

  if(isTRUE(log.p)) logp else exp(logp)
}
