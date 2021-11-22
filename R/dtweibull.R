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
#'
#' @examples
#' dtweibull(3, 0.5, 1.5, 2)
dtweibull = function(x, shape, scale=1, a=0, b=Inf, log=FALSE){
  stopifnot(
    exprs = {
      is.numeric(x); length(x) >= 1L;
      is.numeric(shape); length(shape) >= 1L; all(shape > 0); all(is.finite(shape));
      is.numeric(scale); length(scale) >= 1L; all(scale > 0); all(is.finite(scale));
      is.numeric(a); length(a) >= 1L; all(a >= 0); all(is.finite(a));
      is.numeric(b); length(b) >= 1L; all(b > a);
      is.logical(log); length(log) == 1L
    }
  )

  arg_list = check_args(x = x, shape = shape, scale = scale, a = a, b = b)

  lns = vapply(arg_list, length, 1L)

  ld = function(x, shape, scale, a, b){
    a_shape = a ^ shape
    scale_shape = scale ^ shape
    log(shape) + (shape - 1) * log(x) - shape * log(scale) - (x^shape - a_shape) / scale_shape - log1mexp((b^shape - a_shape) / scale_shape)
  }

  max_lns = max(lns)

  low = rep_len(x < a, max_lns)
  hi = rep_len(x > b, max_lns)
  mid = rep_len(!low & !hi, max_lns)

  log_d = numeric(max_lns)
  if(any(low)) log_d[low] = -Inf
  if(any(hi)) log_d[hi] = 0
  if(any(mid)){
    log_d[mid] = do.call("ld", list_select(x = arg_list, ind = mid))
  }

  if(isTRUE(log)) log_d else exp(log_d)
}
