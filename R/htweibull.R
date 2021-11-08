#' Truncated Weibull hazard function
#'
#' @param x Numerics. Value at which to evaluate hazard function.
#' @param shape Numerics. Shape parameter.
#' @param scale Numerics. Scale parameter (default: 1).
#' @param a Numerics. Lower truncation point (default: 0).
#' @param b Numerics. Upper truncation point (default: Inf).
#' @param log Logical. Return log-hazard (default: FALSE).
#'
#' @return Numerics. Hazard value.
#' @export
#'
#' @examples
#' htweibull(5, 2.5)
htweibull = function(x, shape, scale = 1, a = 0, b = Inf, log = FALSE){
  stopifnot(
    x >= 0, shape > 0, scale > 0, a >= 0, b > a, is.logical(log),
    length(log) == 1L
  )

  args_lst = make_args(x = x, shape = shape, scale = scale, a = a, b = b)

  x_lte_b = x <= b

  x_gte_a = x >= a

  log_h = numeric(length(args_lst$x))

  log_h[!x_lte_b] = Inf

  log_h[!x_gte_a] = -Inf

  log_h[x_gte_a & x_lte_b] = evalq(log(shape) - log(scale) + (shape - 1) * (log(x) - log(scale)) - log1mexp((b / scale)^shape - (x / scale)^shape), envir = list_select(x = args_lst, ind = x_gte_a & x_lte_b))

  if(isTRUE(log)) log_h else exp(log_h)
}
