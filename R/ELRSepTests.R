#' Conduct empirical likelihood ratio tests for a sample of two-way data
#'
#' Empirical likleihood ratio tests for nested separability structures using null-transformed bootstrap calibration
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
#' @param LSmax for finding the MELE under separability, the maximum number of line
#'                search steps to take (default = 100L)
#' @param ELctrl an object of class ControlEL defining options for EL estimation;
#'               see el_control (default is el_control(maxit = 10000L, maxit_l = 100L))
#' @param FVEthres threshold to determine number of basis functions to use when drawing
#'                 bootstrap samples that are similar to the observed data (default is 0.99)
#'
#' @return list with the following elements
#'         - tStats: test statistic value(s) - minus twice the constrained empirical log-likelihood ratio
#'                   corresponding to specified null and alternative
#'         - loglr - log empirical likelihood ratio value for each type of separability
#'         - ELoptInfo: EL optimization information for each test (see optim slot from fitted ELT object).
#'                      For testing weak separability and separability, there is information for both the
#'                      constrained null fit and specific alternative
#'         - bootRes: bootstrap calibration info list (if cal = 'bootstrap', otherwise is NA)
#'            * tStats: list of length B vector of bootstrap test statistics, one vector per test
#'            * loglr: list of length B, each giving the log empirical likelihood ratio for each type of separability
#'            * ELoptInfo: list of lists EL optimization information for each bootstrap, same structure as
#'                         ELoptInfo above for each element
#'            * bootPval: bootstrap p-values, one per test
#' @export

ELRSepTests <- function(X, tt1 = 1:dim(X)[[2]], tt2 = 1:dim(X)[[3]],
                        JTest = 2L, LTest = 2L,
                        nullHyp = c('ParSep', 'WkSep', 'Sep'),
                        B = 500L, LSmax = 100L,
                        ELctrl = melt::el_control(maxit = 10000L, maxit_l = 1000L),
                        FVEthres = 0.99){

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

  ELctrlBoot <- ELctrl;
  ELctrlBoot@verbose = FALSE

  # Set up outputs

  numTests <- length(nullHyp)
  tStats <- rep(NA, numTests);
  ELoptInfo <- loglr <- vector(mode = "list", length = numTests)
  names(tStats) <- names(loglr) <- names(ELoptInfo) <- nullHyp

  # Bootstrap results
  bootRes <- vector(mode = "list", length = 4)
  names(bootRes) <- c('tStats', 'loglr', 'ELoptInfo', 'bootPval')
  bootRes$tStats <- bootRes$loglr <- bootRes$ELoptInfo <- vector(mode = "list", length = numTests)
  bootRes$bootPval <- rep(NA, length = numTests)
  names(bootRes$tStats) <- names(bootRes$loglr) <- names(bootRes$ELoptInfo) <- names(bootRes$bootPval) <- nullHyp

  # Get augmented scores for performing EL Tests

  MBE <- getMBExp(X = X, tt1 = tt1, tt2 = tt2, useFVE = TRUE, FVEthres = FVEthres)
  J <- ncol(MBE$Psi)
  L <- ncol(MBE$Phi)

  if(JTest > J || LTest > L){
    warning('At least one of JTest/LTest is too large - decreasing to largest feasible value')
    JTest <- min(JTest, J)
    LTest <- min(LTest, L)
  }

  scrsAug <- getELTestData(MBE$scrs, J, L, JTest, LTest) # get products of scores
  nVarPS <- LTest * (LTest - 1) * JTest^2 / 2 # number of variables related to ParSep
  nVarWS <- LTest * JTest * (JTest - 1) / 2 # number of variables related to WkSep
  nVarS <- LTest * JTest # number of variables related to Sep
  nVarTot <- nVarPS + nVarWS + nVarS

  indCur <- 0 # tracks which test in nullHyp is currently being conducted

  if('ParSep' %in% nullHyp){ # Run Test for Partial Separability

    indCur <- indCur + 1
    elPS <- melt::el_mean(scrsAug[, 1:nVarPS], rep(0, nVarPS), control = ELctrl)
    convPS <- elPS@optim$convergence

    tStats[indCur] <- ifelse(convPS, elPS@statistic, Inf)
    loglr[[indCur]] <- c(ifelse(convPS, elPS@loglr, -Inf), 0)
    names(loglr[[indCur]]) <- c('ParSep', 'Unconstrained')
    ELoptInfo$ParSep <- elPS@optim

    if(convPS){ # only run bootstrap if the convexity condition holds
      elPSBootTests <- lapply(1:B, function(b){
        Xboot <- drawBoot(MBE, 'ParSep')
        MBEboot <- getMBExp(Xboot, tt1 = tt1, tt2 = tt2, J = JTest, L = LTest)
        scrsAugBoot <- getELTestData(MBEboot$scrs, J = JTest, L = LTest)
        elPSBoot <- melt::el_mean(scrsAugBoot[, 1:nVarPS], rep(0, nVarPS), control = ELctrlBoot)
        return(elPSBoot)
      })
      bootRes$tStats$ParSep <- sapply(elPSBootTests, function(tst) tst@statistic)
      bootRes$loglr$ParSep <- t(sapply(elPSBootTests, function(tst) c(tst@loglr, 0)))
      colnames(bootRes$loglr$ParSep) <- c('ParSep', 'Unconstrained')
      bootRes$EloptInfo$ParSep <- lapply(elPSBootTests, function(tst) tst@optim)
      bootRes$bootPval[indCur] <- (sum(bootRes$tStats$ParSep >= tStats[1]) + 1)/(B + 1)
      } else {
      bootRes$tStats$ParSep <- NA
      bootRes$loglr$ParSep <- NA
      bootRes$ELoptInfo$ParSep <- "No Bootstrap Performed - violation of the convexity condition and likelihood of zero"
      bootRes$bootPval[indCur] <- 0
      }
    }

  if('WkSep' %in% nullHyp){ # Run Test for Weak Separability

    indCur <- indCur + 1

    if(!exists('elPS')){ # in case partially separability has not been tested
      elPS <- melt::el_mean(scrsAug[, 1:nVarPS], rep(0, nVarPS), control = ELctrl)
      convPS <- elPS@optim$convergence
    }

    elWS <- melt::el_mean(scrsAug[, 1:(nVarPS + nVarWS)], par = rep(0, nVarPS + nVarWS),
                  control = ELctrl)
    convWS <- elWS@optim$convergence

    tStats[indCur] <- ifelse(convPS, ifelse(convWS, elWS@statistic - elPS@statistic, Inf), NA)
    loglr[[indCur]] <- c(ifelse(convWS, elWS@loglr, -Inf), ifelse(convPS, elPS@loglr, -Inf))
    names(loglr[[indCur]]) <- c('WkSep', 'ParSep')
    ELoptInfo$WkSep <- elWS@optim

    if(convWS){
      elWSBootTests <- lapply(1:B, function(b){
        Xboot <- drawBoot(MBE, 'WkSep')
        MBEboot <- getMBExp(Xboot, tt1 = tt1, tt2 = tt2, J = JTest, L = LTest)
        scrsAugBoot <- getELTestData(MBEboot$scrs, J = JTest, L = LTest)
        elPSBoot <- melt::el_mean(scrsAugBoot[, 1:nVarPS], par = rep(0, nVarPS), control = ELctrlBoot)
        elWSBoot <- melt::el_mean(scrsAugBoot[, 1:(nVarPS + nVarWS)], par = rep(0, nVarPS + nVarWS),
                            control = ELctrlBoot)
        return(list(elWSBoot, elPSBoot))
      })
      bootRes$tStats$WkSep <- sapply(elWSBootTests, function(tst) tst[[1]]@statistic - tst[[2]]@statistic)
      bootRes$loglr$WkSep <- t(sapply(elWSBootTests, function(tst) c(tst[[1]]@loglr, tst[[2]]@loglr)))
      colnames(bootRes$loglr$WkSep) <- c('WkSep', 'ParSep')
      bootRes$ELoptInfo$WkSep <- lapply(elWSBootTests, function(tst) lapply(tst, function(l) l@optim))
      bootRes$bootPval[indCur] <- (sum(bootRes$tStats$WkSep >= tStats[indCur]) + 1)/(B + 1)
    } else {
      bootRes$tStats$WkSep <- NA
      bootRes$loglr$WkSep <- NA
      bootRes$EloptInfo$WkSep <- "No Bootstrap Performed - violation of the convexity condition and likelihood of zero"
      bootRes$bootPval[indCur] <- ifelse(convPS, 0, NA)
    }
  }

  if('Sep' %in% nullHyp){ # Run Test for Separability

    indCur <- indCur + 1

    # Get starting value based on weakly separable fit

    if(!exists('elWS')){
      elWS <- melt::el_mean(scrsAug[, 1:(nVarPS + nVarWS)], rep(0, nVarPS + nVarWS), control = ELctrl)
      convWS <- elWS@optim$convergence
    }
    SepInd <- (nVarPS + nVarWS + 1):nVarTot
    tmp <- matrix(apply(scrsAug[, SepInd], 2, weighted.mean, w = expm1(elWS@logp) + 1),
                  nrow = JTest, ncol = LTest)
    aStrt <- sum(tmp) # trace of estimated Lambda matrix under weak separability
    gammaStrt <- rowSums(tmp)/aStrt
    betaStrt <- colSums(tmp)/aStrt

    elS <- SepELOpt(scrsAug, JTest, LTest, gammaStrt, betaStrt, aStrt,
                    mOtr = ELctrl@maxit, mInr = ELctrl@maxit_l,
                    tolOtr = ELctrl@tol, tolInr = ELctrl@tol_l,
                    LSmax = LSmax, verb = ELctrl@verbose)
    if(length(elS@optim$convergence) == 1){ # constrained optimization not performed because initial point has likelihood 0
      convS <- FALSE
    } else {
      convS <- elS@optim$convergence$FinalInr
    }

    tStats[indCur] <- ifelse(convWS, ifelse(convS, elS@statistic - elWS@statistic, Inf), NA)
    loglr[[indCur]] <- c(ifelse(convS, elS@loglr, -Inf), ifelse(convWS, elWS@loglr, -Inf))
    names(loglr[[indCur]]) <- c('Sep', 'WkSep')
    ELoptInfo$Sep <- elS@optim


    if(convS){
      elSBootTests <- lapply(1:B, function(b){
        Xboot <- drawBoot(MBE, 'Sep')
        MBEboot <- getMBExp(Xboot, tt1 = tt1, tt2 = tt2, J = JTest, L = LTest)
        scrsAugBoot <- getELTestData(MBEboot$scrs, J = JTest, L = LTest)
        elWSBoot <- melt::el_mean(scrsAugBoot[, 1:(nVarPS + nVarWS)],
                            par = rep(0, nVarPS + nVarWS), control = ELctrlBoot)

          # Now optimize over separable
          tmpBoot <- matrix(apply(scrsAugBoot[, SepInd], 2, weighted.mean, w = expm1(elWSBoot@logp) + 1),
                            nrow = JTest, ncol = LTest)
          aStrtBoot <- sum(tmpBoot) # trace of estimated Lambda matrix under weak separability
          gammaStrtBoot <- rowSums(tmpBoot)/aStrtBoot
          betaStrtBoot <- colSums(tmpBoot)/aStrtBoot

          elSBoot <- SepELOpt(scrsAugBoot, JTest, LTest, gammaStrtBoot, betaStrtBoot, aStrtBoot,
                              mOtr = ELctrl@maxit, mInr = ELctrl@maxit_l,
                              tolOtr = ELctrl@tol, tolInr = ELctrl@tol_l,
                              LSmax = LSmax, verb = FALSE)
          list(elSBoot, elWSBoot)
        })
      bootRes$tStats$Sep <- sapply(elSBootTests, function(tst) tst[[1]]@statistic - tst[[2]]@statistic)
      bootRes$loglr$Sep <- t(sapply(elSBootTests, function(tst) c(tst[[1]]@loglr, tst[[2]]@loglr)))
      colnames(bootRes$loglr$Sep) <- c('Sep', 'WkSep')
      bootRes$ELoptInfo$Sep <- lapply(elSBootTests, function(tst) lapply(tst, function(l) l@optim))
      bootRes$bootPval[indCur] <- (sum(bootRes$tStats$Sep >= tStats[indCur]) + 1)/(B + 1)
    } else {
      bootRes$tStats$Sep <- NA
      bootRes$loglr$Sep <- NA
      bootRes$EloptInfo$Sep <- "No Bootstrap Performed - violation of the convexity condition and likelihood of zero"
      bootRes$bootPval[indCur] <- ifelse(convWS, 0, NA)
    }
  }

  return(list('tStats' = tStats, 'loglr' = loglr, 'ELoptInfo' = ELoptInfo, 'bootRes' = bootRes))

}
