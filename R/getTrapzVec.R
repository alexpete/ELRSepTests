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

  spc <- diff(tt)
  w <- 0.5*(c(0, spc) + c(spc, 0))

}
