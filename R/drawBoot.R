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

  p <- nrow(MBExp$Psi)
  q <- ncol(MBExp$Psi)
  s <- nrow(MBExp$Phi)
  L <- ncol(MBExp$Phi)
  n <- nrow(MBExp$scrs)
  scrs <- MBExp$scrs

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
