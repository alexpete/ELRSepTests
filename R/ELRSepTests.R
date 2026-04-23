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
#' @param JTest integer or increasing integer vector in 1:M1 indicating the first
#'              direction eigenfunctions to test, must be the same length as LTest.
#'              Default is 2L
#' @param LTest integer or increasing integer vector in 1:M2 indicating the second
#'              direction eigenfunctions to test, must be the same length as JTest.
#'              Default is 2L
#' @param nullHyp subset of c('ParSep', 'WkSep', 'Sep') indicating which null
#'                  hypotheses to test (default is all)
#' @param B number of bootstrap samples for computing P values (default is 500)
#' @param LSmax for finding the MELE under separability, the maximum number of line
#'                search steps to take (default = 100L)
#' @param ELctrl an object of class ControlEL defining options for EL estimation;
#'               see el_control (default is el_control(maxit = 10000L, maxit_l = 100L))
#' @param FVEthres threshold to determine number of basis functions to use when drawing
#'                 bootstrap samples that are similar to the observed data (default is 0.99)
#' @param thin logical indicator for thinning output.  If thin = FALSE, then all
#'             output below is returned; if thin = TRUE, the default, then tStatsBoot is set to NA.
#'
#' @return list with the following elements
#'         - tStats: matrix of empirical likelihood test statistic values, with
#'                   rows corresponding to different null hypotheses and columns
#'                   to distinct pairs `(J, L)`
#'         - ELoptInfo: List of EL optimization information for each type of separability
#'                      that corresponds to a null or alternative according to nullHyp.
#'                      See optim slot from fitted ELT object for details.
#'                      Each element itself is a list corresponding to information
#'                      for fitting a given type of separability with secondary index
#'                      corresponding to distinct pairs `(J, L)`.
#'         - bootPval: matrix of bootstrap p-values, with different null hypotheses
#'                     indexing the rows and distinct pairs `(J, L)` indexing columns
#'         - tStatsBoot: list of bootstrap test statistics, `tStatsBoot[[a]]` is a matrix
#'                      of bootstrap test statistics for the corresponding null
#'                      hypothesis with rows indexed by bootstrap samples and
#'                      columns indexed by distinct pairs `(J, L)`. Set to NA if
#'                      thin = TRUE
#' @export

ELRSepTests <- function(X, tt1 = 1:dim(X)[[2]], tt2 = 1:dim(X)[[3]],
                        JTest = 2L, LTest = 2L,
                        nullHyp = c('ParSep', 'WkSep', 'Sep'),
                        B = 500L, LSmax = 100L,
                        ELctrl = melt::el_control(maxit = 10000L, maxit_l = 1000L),
                        FVEthres = 0.99,
                        thin = TRUE){

  # Perform checks

  if(!is.array(X) || length(dim(X)) != 3) stop('X must be a 3D array')

  n <- dim(X)[[1]]; M1 <- dim(X)[[2]]; M2 <- dim(X)[[3]]

  if(length(tt1) != M1 || length(tt2) != M2){
    stop('Lengths of tt1 and tt2 must match second and third dimensions of X')
  }

  if(length(JTest) != length(LTest)){
    stop('LTest and JTest must have the same number of elements')
  }

  if(!(all(is.integer(c(JTest, LTest))) && all(c(JTest, LTest) > 0))){
    stop('LTest and JTest must contain only positive integer values')
  }

  if(any(JTest > M1) || any(LTest > M2)){
    stop('Some values in JTest and LTest are too large')
  }

  if(any(diff(JTest) <=0 ) || any(diff(LTest) <= 0)){
    stop('Value in JTest and LTest must be increasing')
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

  ELctrlBoot <- ELctrl
  ELctrlBoot@verbose = FALSE

  # Get Marginal Basis Expansion

  JTestMax <- max(JTest)
  LTestMax <- max(LTest)
  MBE <- getMBExp(X = X, tt1 = tt1, tt2 = tt2, useFVE = TRUE)
  J <- ncol(MBE$Psi)
  L <- ncol(MBE$Phi)

  if(all(JTest > J) || all(LTest > L)){
    stop('Value in JTest and LTest are all too large as fewer components will explain over 99 percent of the variability')
  }
  if(JTestMax > J || LTestMax > L){
    warning('At least one element JTest/LTest is too large after marginal basis expansions computed - removing infeasible values')
    JTest <- JTest[JTest <= J]
    LTest <- LTest[LTest <= L]
    if(length(JTest) != length(LTest)){
      JTest <- JTest[1:min(length(JTest), length(LTest))]
      LTest <- LTest[1:min(length(JTest), length(LTest))]
    }
    JTest <- unique(c(JTest, J)) # add J if missing
    LTest <- unique(c(LTest, L)) # add L if missing
    JTestMax <- max(JTest) # should be J
    LTestMax <- max(LTest) # should be L
  }

  # Set up outputs

  numTests <- length(nullHyp)
  dimsList <- sapply(1:length(LTest), \(a) paste0('(', JTest[a], ', ', LTest[a], ')'))
  q <- length(dimsList) # = length(JTest) and = length(LTest)
  tStats <- matrix(NA, nrow = numTests, ncol = q)
  rownames(tStats) <- nullHyp
  colnames(tStats) <- dimsList
  bootPval <- tStats

  # Bootstrap results storage
  if(!thin){
    tStatsBoot <- vector(mode = "list", length = numTests)
    names(tStatsBoot) <- nullHyp
  }

  # Get necessary EL Fits

  scrsAug <- getELTestData(MBE$scrs, J, L, JTestMax, LTestMax) # get products of scores
  scrsInd <- getScrsInd(JTestMax, LTestMax) # for indexing columns of scrsAug
  ELRes <- getELFits(scrsAug, nullHyp, JTest, LTest, scrsInd, ELctrl, LSmax)
  ELstats <- ELRes$ELstats
  conv <- ELRes$conv

  if('ParSep' %in% nullHyp){ # Get Test Statistic and Run Bootstrap under ParSep Null

    tStats['ParSep', ] <- ELstats['ParSep', ]

    tStatsBootPS <- vapply(1:B, function(b){
      # Set up outputs as if all likelihoods are zero
      val <- rep(NA, q)

      if(any(conv$ParSep)){ # Run bootstraps for values (J, L) with non-zero likelihood
        convInd <- which(conv$ParSep)
        JTestMaxBoot <- max(JTest[convInd])
        LTestMaxBoot <- max(LTest[convInd])

        Xboot <- drawBoot(MBE, 'ParSep') # Null bootstrap under ParSep
        MBEboot <- getMBExp(Xboot, tt1 = tt1, tt2 = tt2, J = JTestMaxBoot, L = LTestMaxBoot)

        scrsAugBoot <- getELTestData(MBEboot$scrs, J = JTestMaxBoot, L = LTestMaxBoot)
        scrsIndBoot <- getScrsInd(JTestMaxBoot, LTestMaxBoot) # for indexing columns of scrsAugBoot

        ELResBoot <- getELFits(scrsAugBoot, 'ParSep', JTest[convInd], LTest[convInd], scrsIndBoot, ELctrlBoot, LSmax)
        val[convInd] <- ELResBoot$ELstats['ParSep', ]
      }

      return(val)
    }, FUN.VALUE = numeric(q))
    tStatsBootPS <- matrix(tStatsBootPS, nrow = B, ncol = q, byrow = TRUE, dimnames = list(1:B, dimsList))

    bootPval['ParSep', ] <- ifelse(conv$ParSep, sapply(1:q, \(jl) (sum(tStatsBootPS[, jl] > tStats['ParSep', jl]) + 1) / (B + 1)), 0)
    if(!thin) tStatsBoot$ParSep <- tStatsBootPS
  }

  if('WkSep' %in% nullHyp){ # Run Test for Weak Separability

    tStats['WkSep', ] <- ifelse(conv$ParSep, ifelse(conv$WkSep, ELstats['WkSep', ] - ELstats['ParSep', ], Inf), NA)

    tStatsBootWS <- vapply(1:B, function(b){
      # Set up outputs as if all likelihoods are zero
      val <- rep(NA, q)

      if(any(conv$WkSep)){ # Run bootstraps for values (J, L) with non-zero likelihood
        convInd <- which(conv$WkSep)
        JTestMaxBoot <- max(JTest[convInd])
        LTestMaxBoot <- max(LTest[convInd])

        Xboot <- drawBoot(MBE, 'WkSep') # Null bootstrap under WkSep
        MBEboot <- getMBExp(Xboot, tt1 = tt1, tt2 = tt2, J = JTestMaxBoot, L = LTestMaxBoot)

        scrsAugBoot <- getELTestData(MBEboot$scrs, J = JTestMaxBoot, L = LTestMaxBoot)
        scrsIndBoot <- getScrsInd(JTestMaxBoot, LTestMaxBoot) # for indexing columns of scrsAugBoot

        ELResBoot <- getELFits(scrsAugBoot, 'WkSep', JTest[convInd], LTest[convInd], scrsIndBoot, ELctrlBoot, LSmax)
        val[convInd] <- ifelse(ELResBoot$conv$ParSep, ELResBoot$ELstats['WkSep', ] - ELResBoot$ELstats['ParSep', ], Inf)
      }

      return(val)
    }, FUN.VALUE = numeric(q))
    tStatsBootWS <- matrix(tStatsBootWS, nrow = B, ncol = q, byrow = TRUE, dimnames = list(1:B, dimsList))

    bootPval['WkSep', ] <- ifelse(conv$WkSep, sapply(1:q, \(jl) (sum(tStatsBootWS[, jl] > tStats['WkSep', jl]) + 1) / (B + 1)), 0)

    if(!thin) tStatsBoot$WkSep <- tStatsBootWS
  }

  if('Sep' %in% nullHyp){ # Run Test for Separability

    tStats['Sep', ] <- ifelse(conv$WkSep, ifelse(conv$Sep, ELstats['Sep', ] - ELstats['WkSep', ], Inf), NA)

    tStatsBootS <- vapply(1:B, function(b){
      # Set up outputs as if all likelihoods are zero
      val <- rep(NA, q)

      if(any(conv$Sep)){ # Run bootstraps for values (J, L) with non-zero likelihood
        convInd <- which(conv$Sep)
        JTestMaxBoot <- max(JTest[convInd])
        LTestMaxBoot <- max(LTest[convInd])
        Xboot <- drawBoot(MBE, 'Sep') # Null bootstrap under Sep
        MBEboot <- getMBExp(Xboot, tt1 = tt1, tt2 = tt2, J = JTestMaxBoot, L = LTestMaxBoot)

        scrsAugBoot <- getELTestData(MBEboot$scrs, J = JTestMaxBoot, L = LTestMaxBoot)
        scrsIndBoot <- getScrsInd(JTestMaxBoot, LTestMaxBoot) # for indexing columns of scrsAugBoot

        ELResBoot <- getELFits(scrsAugBoot, 'Sep', JTest[convInd], LTest[convInd], scrsIndBoot, ELctrlBoot, LSmax)
        val[convInd] <- ifelse(ELResBoot$conv$WkSep, ELResBoot$ELstats['Sep', ] - ELResBoot$ELstats['WkSep', ], Inf)
      }

      return(val)
    }, FUN.VALUE = numeric(q))
    tStatsBootS <- matrix(tStatsBootS, nrow = B, ncol = q, byrow = TRUE, dimnames = list(1:B, dimsList))
    bootPval['Sep', ] <- ifelse(conv$Sep, sapply(1:q, \(jl) (sum(tStatsBootS[, jl] > tStats['Sep', jl]) + 1) / (B + 1)), 0)

    if(!thin) tStatsBoot$Sep <- tStatsBootS
  }

  res <- list(tStats = tStats, bootPval = bootPval, ELoptInfo = ELRes$ELoptInfo)
  if(!thin) res$tStatsBoot <- tStatsBoot
  return(res)

}

## helper function to get empirical likelihood fits from products of scores

getELFits <- function(scrsAug, nullHyp, JTest, LTest, scrsInd, ELctrl, LSmax){

  dimsList <- sapply(1:length(LTest), \(a) paste0('(', JTest[a], ', ', LTest[a], ')'))
  q <- length(dimsList) # = length(JTest) and = length(LTest)

  fitHyp <- getFitHyp(nullHyp)
  conv <- vector("list", length(fitHyp))
  ELstats <- matrix(NA, nrow = length(fitHyp), ncol = q,
                    dimnames = list(fitHyp, dimsList))
  ELoptInfo <- vector("list", length(fitHyp))
  names(ELoptInfo) <- names(conv) <- fitHyp

  if('ParSep' %in% fitHyp){
    ELFitsPS <- lapply(1:q, \(jl){ # fit EL for each (J, L)
      jlInds <- which((pmax(scrsInd[, 1], scrsInd[, 2]) <= JTest[jl]) & # j and k no bigger than JTest[jl]
                       (pmax(scrsInd[, 3], scrsInd[, 4]) <= LTest[jl]) & # l and m no bigger than LTest[jl]
                       (scrsInd[, 3] < scrsInd[, 4])) # l < m keeps only values related to PS
      scrsAugCur <- scrsAug[, jlInds] # extract relevant columns
      return(melt::el_mean(scrsAugCur, rep(0, ncol(scrsAugCur)), control = ELctrl))
    })
    conv$ParSep <- sapply(ELFitsPS, \(el) el@optim$convergence)
    ELstats['ParSep', ] <- sapply(1:q, \(jl) ifelse(conv$ParSep[jl], ELFitsPS[[jl]]@statistic, Inf))
    ELoptInfo$ParSep <- lapply(ELFitsPS, \(el) el@optim)
  }

  if('WkSep' %in% fitHyp){
    ELFitsWS <- lapply(1:q, \(jl){ # fit EL for each (J, L)
      jlInds <- which((pmax(scrsInd[, 1], scrsInd[, 2]) <= JTest[jl]) & # j and k no bigger than JTest[jl]
                       (pmax(scrsInd[, 3], scrsInd[, 4]) <= LTest[jl]) & # l and m no bigger than LTest[jl]
                       (((scrsInd[, 3] < scrsInd[, 4])) | # all (j, k) are valid if l < m
                          ((scrsInd[, 3] == scrsInd[, 4]) & scrsInd[, 1] < scrsInd[, 2]))) # j < k if l = m
      scrsAugCur <- scrsAug[, jlInds] # extract relevant columns
      return(melt::el_mean(scrsAugCur, rep(0, ncol(scrsAugCur)), control = ELctrl))
    })
    conv$WkSep <- sapply(ELFitsWS, \(el) el@optim$convergence)
    ELstats['WkSep', ] <- sapply(1:q, \(jl) ifelse(conv$WkSep[jl], ELFitsWS[[jl]]@statistic, Inf))
    ELoptInfo$WkSep <- lapply(ELFitsWS, \(el) el@optim)
  }

  if('Sep' %in% fitHyp){

    ELFitsS <- lapply(1:q, \(jl){
      jlInds <- which((pmax(scrsInd[, 1], scrsInd[, 2]) <= JTest[jl]) & # j and k no bigger than JTest[jl]
                       (pmax(scrsInd[, 3], scrsInd[, 4]) <= LTest[jl])) # l and m no bigger than LTest[jl]
      scrsAugCur <- scrsAug[, jlInds] # extract relevant columns
      SepInds <- (ncol(scrsAugCur) - JTest[jl] * LTest[jl] + 1):ncol(scrsAugCur) # last columns are for Sep

      if(conv$WkSep[jl]){
        w <- expm1(ELFitsWS[[jl]]@logp) + 1
      } else {
        w <- rep(1/ncol(scrsAugCur), length.out = ncol(scrsAugCur))
      }
      tmp <- matrix(colSums(scrsAugCur[, SepInds] * w), nrow = JTest[jl], ncol = LTest[jl])
      aStrt <- sum(tmp) # trace of estimated Lambda matrix under weak separability
      gammaStrt <- rowSums(tmp)/aStrt
      betaStrt <- colSums(tmp)/aStrt

      return(SepELOpt(scrsAugCur, JTest[jl], LTest[jl], gammaStrt, betaStrt, aStrt,
                      mOtr = ELctrl@maxit, mInr = ELctrl@maxit_l,
                      tolOtr = ELctrl@tol, tolInr = ELctrl@tol_l,
                      LSmax = LSmax, verb = ELctrl@verbose))
    })

    conv$Sep <- sapply(ELFitsS, \(el){
      if(length(el@optim$convergence) == 1){ # constrained optimization not performed because initial point has likelihood 0
        return(FALSE)
      } else {
      return(el@optim$convergence$FinalInr)
      }
    })
    ELstats['Sep', ] <- sapply(1:q, \(jl) ifelse(conv$Sep[jl], ELFitsS[[jl]]@statistic, Inf))
    ELoptInfo$Sep <- lapply(ELFitsS, \(el) el@optim)

  }

  return(list(ELstats = ELstats, conv = conv, ELoptInfo = ELoptInfo))

}

# Helper to determine which separability hypotheses to fit given nullHyp

getFitHyp <- function(nullHyp) {

  allHyp <- c("ParSep", "WkSep", "Sep")
  fitHyp <- nullHyp

  if ("WkSep" %in% nullHyp) {
    fitHyp <- union(fitHyp, "ParSep")
  }

  if ("Sep" %in% nullHyp) {
    fitHyp <- union(fitHyp, "WkSep")
  }

  allHyp[allHyp %in% fitHyp]
}

## helper function to index columns of scrs corresponding to separability ordering

getScrsInd <- function(J, L) {

  indList <- list()
  ctr <- 1L

  ## 1. l < m
  if (L >= 2L) {
    for (l in 1:(L - 1L)) {
      for (m in (l + 1L):L) {
        for (j in 1:J) {
          for (k in 1:J) {
            indList[[ctr]] <- c(j = j, k = k, l = l, m = m)
            ctr <- ctr + 1L
          }
        }
      }
    }
  }

  ## 2. l = m, j < k
  if (J >= 2L) {
    for (l in 1:L) {
      for (j in 1:(J - 1L)) {
        for (k in (j + 1L):J) {
          indList[[ctr]] <- c(j = j, k = k, l = l, m = l)
          ctr <- ctr + 1L
        }
      }
    }
  }

  ## 3. l = m, j = k  (diagonal terms)
  for (l in 1:L) {
    for (j in 1:J) {
      indList[[ctr]] <- c(j = j, k = j, l = l, m = l)
      ctr <- ctr + 1L
    }
  }

  if (length(indList) == 0L) {
    indMat <- matrix(integer(0), nrow = 0L, ncol = 4L)
    colnames(indMat) <- c("j", "k", "l", "m")
    return(indMat)
  }

  indMat <- do.call(rbind, indList)
  storage.mode(indMat) <- "integer"
  indMat
}





