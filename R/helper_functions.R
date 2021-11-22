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


#' Numerically compute the value of log(1 - p + p * exp(-x))
#'
#' Due to numerical accuracy issues with directly computing this
#' function, a different approach is used that is hopefully
#' more robust. The approach uses uniroot, and therefore
#' may fail. This is an internal function to the package and no
#' real error checking is performed on its inputs.
#'
#' @param p (Numeric) Numeric 0 <= p <= 1
#' @param shape (Numeric) Weibull shape parameter
#' @param scale (Numeric) Weibull scale parameter
#' @param a (Numeric) Lower truncation point
#' @param b (Numeric) Upper truncation point
#' @param ... Arguments that get passed along to `uniroot`
#'
#' @return A numeric value
#' @importFrom VGAM log1mexp
#' @importFrom utils modifyList
#'
#' @keywords internal
log1mpexp = function(log_p, a, b, shape, scale, ...){
  x = (a^shape - b^shape) / scale^shape

  if(is.infinite(log_p) && log_p < 0)
    return(0)

  if(log_p == 0){
    return(x)
  }

  if(is.infinite(x) && x < 0){
    return(log1mexp(-log_p))
  }

  root_fun = function(y){
    log1mexp(-y) - log_p - log1mexp(-x)
  }

  extra_args = modifyList(
    list(tol = .Machine$double.eps),
    list(...)
  )

  root = do.call(
    what = "uniroot",
    args = c(
      list(
        f = root_fun,
        lower = x,
        upper = 0,
        extendInt = "downX"),
        extra_args
      )
    )

  root$root
}

#' Check arguments to function
#'
#' This function does some basic checks of the arguments going
#' into one of the truncated Weibull functions. It then returns
#' a list containing those arguments.
#'
#' @param ... Numerics. Truncated Weibull parameters.
#'
#' @param expand Logical. Whether to expand arguments to same length (default: FALSE).
#'
#' @return list
#' @keywords internal
check_args = function(...,expand = FALSE){
  args_lst = list(...)

  if(length(names(args_lst)) != length(args_lst)) stop("All arguments to 'make_args' must be named", call. = FALSE)

  lens = vapply(args_lst, length, 1L)

  if(!all(lens) > 0L) stop("No arguments can have 0 length", call. = FALSE)
  if(!all(lens == 1L) & length(unique(lens)) > 2L) stop("Arguments must have same length or be of length 1", call. = FALSE)

  max_len = max(lens)

  if(isTRUE(expand)){
    lapply(
      args_lst,
      function(x) {if(length(x) < max_len) rep_len(x, max_len) else x}
    )
  } else {
    args_lst
  }
}
