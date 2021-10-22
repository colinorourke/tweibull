#' Compute Denominator for a Truncated Weibull Distribution
#'
#' This is an internal function.
#'
#' @param a Lower interval endpoint
#' @param b Upper interval endpoint
#' @param shape Shape parameter
#' @param scale Scale parameter
#' @importFrom VGAM log1mexp
#' @keywords internal
#'
#' @return A numeric vector
ldenom = function(a, b, shape, scale){
  log1mexp((b^shape - a^shape)/scale^shape) - (a/scale)^shape
}

#' Select Elements from List of Vectors
#'
#' This function selects specified elements from vectors contained in a list.
#' For vectors of length 1 it just returns that value. The assumption when
#' this was written is that all vectors have a common length or length 1.
#'
#' @param ... Names arguments to make into a list
#' @param x List to select elements from
#' @param ind Logical indicators of which elements to choose
#'
#' @return Modified list
#' @keywords internal
list_select = function(..., x, ind){
  dots = list(...)
  lapply(if(missing(x)) dots else x, function(x, ind) if(length(x) > 1L) x[ind] else x, ind = ind)
}
