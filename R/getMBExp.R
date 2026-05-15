#' Compute components of marginal tensor basis expansion of two-way data
#'
#' Uses functionality of covsep package to find the marginal bases, then expands
#' the data covariance in its tensor basis
#'
#' @param X n-by-M1-by-M2 array of data values.  The first direction indexes the
#'          n observations, with `X[i,,]` consisting of a M1-by-M2 discretization of the
#'          two-way data for the i-th observational unit
#' @param tt1 optional observation grid if the second direction of X indexes
#'            to functional data of length M1 (default is 1:M1)
#' @param tt2 optional observation grid if the second direction of X indexes
#'            to functional data of length M2 (default is 1:M2)
#' @param J number of eigenvalues to extract in the first direction. Will attempt
#'          to find J eigenvalues, but will only retain those that are positive.
#'          Default is M1, but this will almost always be too many.
#' @param L number of eigenvalues to extract in the second direction. Will attempt
#'          to find L eigenvalues, but will only retain those that are positive.
#'          Default is M2, but this will almost always be too many.
#' @param useFVE logical indicating whether or not to use cumulative FVEs to determine
#'               J and L.  If true, the values of `JFVE` and `LFVE` that achieve at
#'               least 1 - (1 - FVEthres)/2 are computed, and the number of eigenvalues
#'                extracted are `Jext = max(J, JFVE)` and `Lext = max(L, LFVE)` (default is FALSE)
#' @param FVEthres threshold for overall FVE required (default is 0.99)
#'
#' @return list with five elements
#'            Psi - M1-by-J matrix of basis estimates for the first partial trace
#'                  operator (second direction is traced out), where J is the
#'                  number of positive eigenvalues retained
#'            Phi - M2-by-L matrix of basis estimates for the second partial trace
#'                  operator (first direction is traced out), where L is the
#'                  number of positive eigenvalues retained
#'            scrs - n-by-JL array of scores with columns indexed in mirror
#'                   dictionary order (1, 1), (2, 1),...(1, L),...(J, L).  That is,
#'                   if column m is associated with pair (j, l), then `scrs[i, m]`
#'                   is the score of (centred) `X[i,,]` in the direction
#'                   of the tensor product of `Psi[j,]` and `Phi[,l]`
#'            Xbar - M1-by-M2 cross-sectional mean of X across first index
#'            Lambda - JL-by-JL covariance matrix with row and column indices
#'                     corresponding to the column ordering of scrs
#'            cumFVE - List of length two containing cumulative FVEs for first and
#'                     second directions.  The length of each list is equal to the
#'                     number of positive eignevalues of the corresponding partial
#'                     trace operator and may be greater than J or L
#' @export

getMBExp <- function(X, tt1 = 1:dim(X)[[2]], tt2 = 1:dim(X)[[3]], J = length(tt1),
                     L = length(tt2), useFVE = FALSE, FVEthres = 0.99) {

  if(!is.array(X) || length(dim(X)) != 3) stop('X must be a 3D array')

  n <- dim(X)[[1]]; M1 <- dim(X)[[2]]; M2 <- dim(X)[[3]]

  if(length(tt1) != M1 || length(tt2) != M2){
    stop('Lengths of tt1 and tt2 must match second and third dimensions of X')
  }

  if(J > M1 || L > M2){
    warning('J/L cannot be larger than M1/M2 - restting to default')
    J <- M1; L <- M2
  }

  ## Rescale X if either of the dimensions are functional
  d1 <- diff(tt1); d2 <- diff(tt2)
  if(any(d1 != 1)){
    w1 <- getTrapzVec(tt1)
  } else {
    w1 <- rep(1, M1)
  }

  if(any(d2 != 1)){
    w2 <- getTrapzVec(tt2)
  } else {
    w2 <- rep(1, M2)
  }

  Xbar <- apply(X, c(2, 3), mean)
  Xc <- sweep(X, c(2, 3), Xbar)

  # Y is rescaled and centered X to do PCA appropriately if either direction is functional
  if(any(c(d1, d2) != 1)){
    W <- tcrossprod(sqrt(w1), sqrt(w2))
    Y <- Xc * aperm(array(W, dim = c(M1, M2, n)), c(3, 1, 2))
  } else {
    Y <- Xc
  }

  ## Get the partial eigendecompositions
  PartialTrKers <- covsep::marginal_covariances(Y)
  KerMerc <- lapply(PartialTrKers, eigen, symmetric = TRUE) # Mercer decomposition of kernels

  ## Get cumulative FVEs
  cumFVE <- lapply(KerMerc, function(K) cumsum(K$values[K$values > 0]/sum(K$values[K$values > 0])))
  names(cumFVE) <- c('Dim1', 'Dim2')

  if(useFVE){
    th <- 1 - (1 - FVEthres)/2
    JFVE <- min(which(cumFVE[[1]] > th))
    LFVE <- min(which(cumFVE[[2]] > th))
    J <- min(max(J, JFVE), length(cumFVE$Dim1)) # smaller of max(J, JFVE) and number of positive eigenvalues
    L <- min(max(L, LFVE), length(cumFVE$Dim2)) # smaller of max(L, LFVE) and number of positive eigenvalues
  } else {
    J <- min(J, length(cumFVE$Dim1)) # smaller of J and number of positive eigenvalues
    L <- min(L, length(cumFVE$Dim2)) # smaller of L and number of positive eigenvalues
  }


  ## Get eigenvectors
  Psi <- KerMerc[[1]]$vectors[, 1:J]
  Phi <- KerMerc[[2]]$vectors[, 1:L]

  ## Get coefficients of product basis expansion
  scrs <- t(matrix(aperm(Y, c(2, 3, 1)), nrow = M1 * M2)) %*% kronecker(Phi, Psi)
  Lambda <- cov(scrs)*(n - 1)/n

  if(any(d1 != 1)) Psi <- Psi * (1/sqrt(w1))
  if(any(d2 != 1)) Phi <- Phi * (1/sqrt(w2))

  return(list('Psi' = Psi, 'Phi' = Phi, 'scrs' = scrs, 'Xbar' = Xbar, 'Lambda' = Lambda, 'cumFVE' = cumFVE))
}
