#' Compute trapezoidal integraion weights for grid vector tt
#'
#' If f is a vector representing discrete functional observations at the points
#' in tt, then sum(w*f) is the trapezoidal integral approximation
#'
#' @param tt observation grid for the functional data
#'
#' @return vector of same length as tt giving weights for trapezoidal integration,
#'

getTrapzVec <- function(tt) {

  if (!is.numeric(tt) || !is.vector(tt)) {
    stop("tt must be a numeric vector")
  }

  if (length(tt) < 2) {
    stop("tt must have at least 2 points")
  }

  if (any(!is.finite(tt))) {
    stop("tt contains NA/NaN/Inf")
  }

  if (any(diff(tt) <= 0)) {
    stop("tt must be strictly increasing")
  }

  spc <- diff(tt)
  w <- 0.5*(c(0, spc) + c(spc, 0))

  return(w)

}
