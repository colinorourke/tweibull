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
#' @importFrom VGAM log1mexp
#' @export
#'
#' @examples
#' ptweibull(2, shape=1.5, a=1)
ptweibull = function(q, shape, scale=1, a=0, b = Inf, lower.tail=TRUE, log.p=FALSE){
  stopifnot(
    exprs = {
      is.numeric(q); is.numeric(shape); is.numeric(scale); is.numeric(a); is.numeric(b); is.logical(lower.tail); is.logical(log.p);
      length(log.p) == 1L;
      q >= 0; shape > 0; scale > 0; a >= 0; all(b > a)
    }
  )

  arg_list = check_args(
    q = q,
    shape = shape,
    scale = scale,
    a = a,
    b = b,
    lower.tail = lower.tail
  )

  lns = vapply(arg_list, length, 1L)
  max_lns = max(lns)

  pfun = function(q, shape, scale, a, b, lower.tail){
    lt_ind = as.numeric(lower.tail)

    ldenom = log1mexp((b^shape - a^shape)/scale^shape)

    lt_ind * (log1mexp((q^shape - a^shape) / scale^shape) - ldenom) +
      (1 - lt_ind) * ((a^shape - q^shape) / scale^shape + log1mexp((b^shape - q^shape)/scale^shape) - ldenom)
  }

  low = rep_len(q <= a, max_lns)
  hi = rep_len(q >= b, max_lns)
  mid = rep_len(!low & !hi, max_lns)

  logp = numeric(max_lns)
  if(any(low)) logp[low] = ifelse(lower.tail, -Inf, 0)
  if(any(hi)) logp[hi] = ifelse(lower.tail, 0, -Inf)
  if(any(mid)) logp[mid] = do.call("pfun", list_select(x = arg_list, ind = mid))

  if(isTRUE(log.p)) logp else exp(logp)
}
