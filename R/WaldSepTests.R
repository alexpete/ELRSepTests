#' Conduct Wald-type separability tests for a sample of two-way data
#'
#' Wald tests of Aston et al. (full studentization) and Lynch and Chen (no studentization)
#' for assessing separability structures using bootstrap calibration
#'
#' @param X n-by-M1-by-M2 array of data values.  The first direction indexes the
#'          n observations, with `X[i,,]` consisting of a M1-by-M2 discretization of the
#'          two-way data for the i-th observational unit
#' @param tt1 optional observation grid if the second direction of X indexes
#'            to functional data of length M1 (default is 1:M1)
#' @param tt2 optional observation grid if the second direction of X indexes
#'            to functional data of length M2 (default is 1:M2)
#' @param JTest number of projections to test in the first direction, cannot be more
#'              than M1 (default is 2)
#' @param LTest number of projections to test in the second direction, cannot be more
#'              than M2 (default is 2)
#' @param nullHyp subset of c('ParSep', 'WkSep', 'Sep') indicating which null
#'                  hypotheses to test (default is all)
#' @param B number of bootstrap samples for computing P values (default is 500)
#'
#' @return list with the following elements
#'         - tStats: test statistic value(s) - squared length of (studentized)
#'                   vector of Wald statistics corresponding to specified null
#'                   (note that these will BE NA for testing separability since
#'                   covsep doesn't return them)
#'         - bootRes:
#'            * tStats: list of length B vector of bootstrap test statistics, one vector per test
#'                      (note that these will be NA for testing separability since covsep doesn't return them)
#'            * bootPval: bootstrap p-values, one per test
#' @export

WaldSepTests <- function(X, tt1 = 1:dim(X)[[2]], tt2 = 1:dim(X)[[3]],
                         JTest = 2L, LTest = 2L,
                         nullHyp = c('ParSep', 'WkSep', 'Sep'),
                         B = 500L){

  # Perform checks

  if(!is.array(X) || length(dim(X)) != 3) stop('X must be a 3D array')

  n <- dim(X)[[1]]; M1 <- dim(X)[[2]]; M2 <- dim(X)[[3]]

  if(length(tt1) != M1 || length(tt2) != M2){
    stop('Lengths of tt1 and tt2 must match second and third dimensions of X')
  }

  if(JTest > M1 || LTest > M2){
    warning('J/L cannot be larger than M1/M2 - restting to default')
    J <- M1; L <- M2
  }

  nullList <- c('ParSep', 'WkSep', 'Sep')
  tmp <- nullHyp %in% nullList
  if(!all(tmp)){
    if(!any(tmp)){
      warning('All provided values in nullHyp are invalid - resetting to default')
      nullHyp <- c('ParSep', 'WkSep', 'Sep')
    } else {
      warning('Removing invalid elements of nullHyp')
      nullHyp <- nullHyp[tmp]
    }
  }

  # Set up outputs

  numTests <- length(nullHyp)
  tStats <- rep(NA, numTests)
  names(tStats) <- nullHyp

  # Bootstrap results
  bootRes <- vector(mode = "list", length = 2)
  names(bootRes) <- c('tStats', 'bootPval')
  bootRes$tStats <- matrix(NA, nrow = B, ncol = numTests)
  colnames(bootRes$tStats) <- nullHyp
  bootRes$bootPval <- rep(NA, numTests)
  names(bootRes$bootPval) <- nullHyp

  # Get lambda estimates for performing weak Wald Tests

  MBE <- getMBExp(X = X, tt1 = tt1, tt2 = tt2, J = JTest, L = LTest)
  J <- ncol(MBE$Psi)
  L <- ncol(MBE$Phi)

  if(JTest > J || LTest > L){
    warning('At least one of JTest/LTest is too large - decreasing to largest feasible value')
    JTest <- min(JTest, J)
    LTest <- min(LTest, L)
  }

  Lambda <- MBE$Lambda

  # Run Partial Weak Tests (Lynch and Chen modification)
  doPS <- 'ParSep' %in% nullHyp
  doWS <- 'WkSep' %in% nullHyp
  doS <- 'Sep' %in% nullHyp

  testCur <- 0 # tracks which null hypothesis we are doing

  if(doPS || doWS){

    if(doPS){
      lamP <- getPSlambdas(Lambda, JTest, LTest)
      lamBootP <- matrix(NA, nrow = JTest * JTest * LTest * (LTest - 1) / 2, ncol = B)
    }

    if(doWS){
      lamW <- getWSlambdas(Lambda, JTest, LTest)
      lamBootW <- matrix(NA, nrow = LTest * JTest * (JTest - 1) / 2, ncol = B)
    }

    for(b in 1:B){
      XBoot <- X[sample.int(n, n = n, replace = TRUE),,]
      lam <- getMBExp(X = XBoot, tt1 = tt1, tt2 = tt2, J = JTest, L = LTest)$Lambda
      if(doPS) lamBootP[, b] <- getPSlambdas(lam, JTest, LTest)
      if(doWS) lamBootW[, b] <- getWSlambdas(lam, JTest, LTest)
    }

    if(doPS){

      testCur <- testCur + 1

      tStats[testCur] <- n * sum(lamP * lamP)
      bootRes$tStats[, testCur] <- sapply(1:B, \(b) n * sum((drop(lamBootP[, b]) - lamP)^2))
      bootRes$bootPval[testCur] <- (sum(bootRes$tStats[, testCur] > tStats[testCur]) + 1) / (B + 1)

    }

    if(doWS){

      testCur <- testCur + 1

      tStats[testCur] <- n * sum(lamW * lamW)
      bootRes$tStats[, testCur] <- sapply(1:B, \(b) n * sum((drop(lamBootW[, b]) - lamW)^2))
      bootRes$bootPval[testCur] <- (sum(bootRes$tStats[, testCur] > tStats[testCur]) + 1) / (B + 1)

    }
  }

  if(doS){

    testCur <- testCur + 1

    sepTest <- covsep::empirical_bootstrap_test(Data = X, L1 = JTest, L2 = LTest,
                                                  studentize = 'full', B = B, verbose = FALSE)
    bootRes$bootPval[testCur] <- sepTest
  }

  return(list('tStats' = tStats, 'bootRes' = bootRes))
}


## helper functions to extract elements of Lambda corresponding to Partial/Weak Separability

getPSlambdas <- function(Lambda, J, L){

  v <- lapply(1:(L - 1L), function(l){
    u <- lapply((l + 1L):L, function(m){
      w <- Lambda[((l - 1L) * J + 1L):(l * J),((m - 1L) * J + 1L):(m * J)]
      return(c(t(w)))
    })
    return(u)
  })

  return(unlist(v, use.names = FALSE))

}

getWSlambdas <- function(Lambda, J, L){

  v <- lapply(1:L, function(l){
    ind <- ((l - 1L) * J + 1L):(l * J)
    u <- Lambda[ind, ind, drop = FALSE]
    u[row(u) < col(u)]
  })

  return(unlist(v, use.names = FALSE))

}
