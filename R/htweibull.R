#' Truncated Weibull hazard function
#'
#' @param x Numerics. Value at which to evaluate hazard function.
#' @param shape Numerics. Shape parameter.
#' @param scale Numerics. Scale parameter (default: 1).
#' @param a Numerics. Lower truncation point (default: 0).
#' @param b Numerics. Upper truncation point (default: Inf).
#' @param log Logical. Return log-hazard (default: FALSE).
#'
#' @return Numerics. Hazard function value.
#' @export
#'
#' @examples
#' htweibull(5, 2.5)
htweibull = function(x, shape, scale = 1, a = 0, b = Inf, log = FALSE){
  stopifnot(
    exprs = {
      is.numeric(x); length(x) >= 1L; is.finite(x);
      is.numeric(shape); length(shape) >= 1L; all(shape > 0); all(is.finite(shape));
      is.numeric(scale); length(scale) >= 1L; all(scale > 0); all(is.finite(scale));
      is.numeric(a); length(a) >= 1L; all(a >= 0); all(is.finite(a));
      is.numeric(b); length(b) >= 1L; all(b > a);
      is.logical(log); length(log) == 1L
    }
  )

  args_lst = check_args(x = x, shape = shape, scale = scale, a = a, b = b)

  lns = vapply(args_lst, length, 1L)

  max_lns = max(lns)

  x_lte_b = rep_len(x <= b, max_lns)

  x_gte_a = rep_len(x >= a, max_lns)

  log_h = numeric(max_lns)

  log_h[!x_lte_b] = Inf

  log_h[!x_gte_a] = -Inf

  log_h[x_gte_a & x_lte_b] = evalq(log(shape) - log(scale) + (shape - 1) * (log(x) - log(scale)) - log1mexp((b / scale)^shape - (x / scale)^shape), envir = list_select(x = args_lst, ind = x_gte_a & x_lte_b))

  if(isTRUE(log)) log_h else exp(log_h)
}

#' Truncated Weibull cumulative hazard
#'
#' @param x (Numerics) Place to evaluate cumulative hazard
#' @param shape (Numerics) Weibull shape parameter
#' @param scale (Numerics) Weibull scale parameter (default: 1)
#' @param a (Numerics) Lower truncation point (default: 0)
#' @param b (Numerics) Upper truncation point (default: Inf)
#' @param log (Logical) Return log-cumhaz (default: FALSE)
#'
#' @return (Numerics) Value of cumulative hazard
#' @export
#'
#' @examples
#' Htweibull(3, 2.5, 1.5, 1, 5)
Htweibull = function(x, shape, scale = 1, a = 0, b = Inf, log=FALSE){
  stopifnot(
    exprs = {
      is.numeric(x); is.numeric(shape); is.numeric(scale); is.numeric(a); is.numeric(b); is.logical(log);
      length(x) >= 1L; length(shape) >= 1L; length(scale) >= 1L; length(a) >= 1L; length(b) >= 1L; length(log) == 1L;
      all(!is.infinite(x)); all(shape > 0); all(!is.infinite(a)); all(b > a);
    }
  )

  args = check_args(x = x, shape = shape, scale = scale, a = a, b = b)

  lns = sapply(args, length)

  max_lns = max(lns)

  h_res = numeric(max_lns)

  x_lte_a = rep_len(x <= a, max_lns)
  x_gte_b = rep_len(x >= b, max_lns)

  h_res[x_lte_a] = 0
  h_res[x_gte_b] = Inf

  h_res[!x_lte_a & !x_gte_b] = evalq((x^shape - a^shape) / scale^shape + log1mexp((b^shape - a^shape) / scale^shape) - log1mexp((b^shape - x^shape) / scale^shape), envir = list_select(x = args, ind = !x_lte_a & !x_gte_b))

  if(isTRUE(log)) log(h_res) else h_res
}

