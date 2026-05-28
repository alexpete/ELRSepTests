#' Draw a null-transformed bootstrap sample
#'
#' Modifes the marginal basis expansion to match a particular separability structure and appropriate bootstrap sample
#' @param MBExp list output of getMBExp on original data
#' @param nullType one of 'ParSep', 'WkSep', 'Sep' indicating the null hypothesis
#'                   under which to draw the bootstrap sample
#' @param mnBoot integer no larger than n giving the size of bootstrap samples
#'               to be drawn - an m-out-of-n bootstrap.
#'               for a full bootstrap
#'
#' @return 3D array of null-bootstrapped two-way data with dimensions mnBoot-by-p-by-s
#'
drawBoot <- function(MBExp, nullType, mnBoot){

  # Core Input Checks
  if (!is.list(MBExp)) {
    stop("MBExp must be a list")
  }

  required_names <- c("Psi", "Phi", "scrs", "Lambda")
  if (!all(required_names %in% names(MBExp))) {
    stop("MBExp is missing required components: Psi, Phi, scrs, Lambda")
  }

  # Finite Value Checks

  check_finite <- function(x, name) {
    if (any(!is.finite(x))) {
      stop(paste(name, "contains NA, NaN, or Inf values"))
    }
  }

  check_finite(MBExp$Psi, "Psi")
  check_finite(MBExp$Phi, "Phi")
  check_finite(MBExp$scrs, "scrs")
  check_finite(MBExp$Lambda, "Lambda")

  # Dimension consistency checks
  # Psi dimensions
  if (nrow(MBExp$Psi) <= 0 || ncol(MBExp$Psi) <= 0) {
    stop("Psi has invalid dimensions")
  }

  # Phi dimensions
  if (nrow(MBExp$Phi) <= 0 || ncol(MBExp$Phi) <= 0) {
    stop("Phi has invalid dimensions")
  }

  # Bootstrap Parameter Checks
  n <- nrow(MBExp$scrs)

  if (!is.numeric(mnBoot) ||
      length(mnBoot) != 1 ||
      is.na(mnBoot) ||
      !is.finite(mnBoot) ||
      abs(mnBoot - round(mnBoot)) > 1e-10 ||
      mnBoot <= 0) {

    stop("mnBoot must be a single finite positive integer value")
  }

  if (mnBoot <= 0) {
    stop("mnBoot must be greater than 0")
  }

  if (mnBoot > n) {
    stop("mnBoot cannot exceed number of observations n")
  }

  p <- nrow(MBExp$Psi)
  q <- ncol(MBExp$Psi)
  s <- nrow(MBExp$Phi)
  L <- ncol(MBExp$Phi)
  n <- nrow(MBExp$scrs)
  scrs <- MBExp$scrs

  if (ncol(MBExp$scrs) != q * L) {
    stop("scrs columns must equal ncol(Psi) * ncol(Phi)")
  }

  valid_types <- c("ParSep", "WkSep", "Sep")

  if (!is.character(nullType) || length(nullType) != 1) {
    stop("nullType must be a single character string")
  }

  if (!(nullType %in% valid_types)) {
    stop("nullType must be one of: ParSep, WkSep, Sep")
  }

  # Lambda Specific Checks
  # Lambda must be a matrix
  if (!is.matrix(MBExp$Lambda)) {
    stop("Lambda must be a matrix")
  }

  # Diagonal must be positive for scaling
  if (any(diag(MBExp$Lambda) <= 0)) {
    stop("Lambda diagonal must contain strictly positive values")
  }

  switch(nullType,
         ParSep = { # sample across different functional bases independently
           scrsBoot <- lapply(1:L, function(l){
             bootInd <- sample.int(n, size = mnBoot, replace = TRUE)
             lind <- (l - 1)*q + 1:q
             return(scrs[bootInd, lind]) # sample column block jointly
           })
           scrsBoot <- do.call(cbind, scrsBoot)
         },
         WkSep = { # sample all coefficients independently
           scrsBoot <- sapply(1:(q*L), function(r){
             bootInd <- sample.int(n, size = mnBoot, replace = TRUE)
             return(scrs[bootInd, r])
           })
         },
         Sep = { # sample all (rescaled) coefficients independently
           LamDiag <- diag(MBExp$Lambda) # estimated score variances
           tmpMat <- matrix(LamDiag, nrow = q, ncol = L)
           LamSepDiag <- kronecker(colSums(tmpMat), rowSums(tmpMat))/sum(tmpMat)
           scrsBoot <- sapply(1:(q*L), function(r){
             bootInd <- sample.int(n, size = mnBoot, replace = TRUE)
             return(scrs[bootInd, r]*sqrt(LamSepDiag[r]/LamDiag[r]))
           })
         }
        )

  Xboot <- aperm(array(tcrossprod(kronecker(MBExp$Phi, MBExp$Psi), scrsBoot),
                       dim = c(p, s, mnBoot)), c(3, 1, 2))

  return(Xboot)
}
