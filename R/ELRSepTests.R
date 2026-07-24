#' Conduct empirical likelihood ratio tests for a sample of two-way data
#'
#' Empirical likleihood ratio tests for nested separability structures using null-transformed bootstrap calibration
#'
#' @param X n-by-M1-by-M2 array of data values.  The first direction indexes the
#'          n observations, with `X[i,,]` consisting of a M1-by-M2 discretization of the
#'          two-way data for the i-th observational unit
#' @param tt1 optional observation grid if the second direction of X indexes
#'            to functional data of length M1 (default is NULL)
#' @param tt2 optional observation grid if the second direction of X indexes
#'            to functional data of length M2 (default is NULL)
#' @param JTest integer or nondecreasing integer vector in 1:M1 indicating the first
#'              direction eigenfunctions to test, must be the same length as LTest.
#'              Default is 2L
#' @param LTest integer or nondecreasing integer vector in 1:M2 indicating the second
#'              direction eigenfunctions to test, must be the same length as JTest.
#'              Default is 2L
#' @param nullHyp subset of c('ParSep', 'WkSep', 'Sep') indicating which null
#'                  hypotheses to test (default is all)
#' @param B number of bootstrap samples for computing P values (default is 500)
#' @param mnBoot integer no larger than n giving the size of bootstrap samples
#'               to be drawn - an m-out-of-n bootstrap.  Default is `mnBoot = n`
#'               for a full bootstrap
#' @param calibration subset of c('elBoot', 'quadForm') indicating which bootstrap
#'               calibration method(s) to compute (default is `'elBoot'` only, matching
#'               earlier package versions). `'elBoot'` calibrates by re-running EL
#'               optimization on each null-constrained bootstrap sample. `'quadForm'`
#'               instead estimates `V0` (and, for the separable null, `U0`) once from
#'               the full sample and, for each bootstrap replicate, plugs the bootstrap
#'               score-product mean vector into the fixed quadratic forms of Theorem 3
#'               of the manuscript - avoiding EL/constrained re-optimization (including
#'               the nested outer optimization otherwise required under the separable
#'               null) inside the bootstrap loop entirely. The observed statistics in
#'               `tStats` do not depend on `calibration`. If both methods are requested,
#'               they are computed from the *same* bootstrap draw at each replicate (no
#'               extra resampling cost), and `bootPval`/`tStatsBoot` are then returned as
#'               named lists keyed by calibration method; with a single calibration
#'               method (the default), their structure is unchanged from earlier package
#'               versions (see `Value`).
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
#'                     indexing the rows and distinct pairs `(J, L)` indexing columns.
#'                     If `length(calibration) > 1`, this is instead a named list of
#'                     such matrices, one per element of `calibration` (e.g.
#'                     `bootPval$elBoot`, `bootPval$quadForm`)
#'         - tStatsBoot: list of bootstrap test statistics, `tStatsBoot[[a]]` is a matrix
#'                      of bootstrap test statistics for the corresponding null
#'                      hypothesis with rows indexed by bootstrap samples and
#'                      columns indexed by distinct pairs `(J, L)`. Set to NA if
#'                      thin = TRUE. If `length(calibration) > 1`, this is instead a
#'                      named list keyed by calibration method, each element having
#'                      this same list-by-null-hypothesis structure (e.g.
#'                      `tStatsBoot$quadForm$ParSep`)
#'         - V0: only present when `'quadForm' %in% calibration`. List (indexed by
#'               distinct pairs `(J, L)`) of full-sample V0 estimates (see Theorem 3
#'               of the manuscript)
#'         - U0: only present when `'quadForm' %in% calibration` and `'Sep' %in% nullHyp`.
#'               List (indexed by distinct pairs `(J, L)`) of full-sample U0 estimates
#' @export

ELRSepTests <- function(X, tt1 = NULL, tt2 = NULL,
                        JTest = 2L, LTest = 2L,
                        nullHyp = c('ParSep', 'WkSep', 'Sep'),
                        B = 500L, mnBoot = dim(X)[1],
                        calibration = 'elBoot',
                        LSmax = 100L,
                        ELctrl = melt::el_control(maxit = 10000L, maxit_l = 1000L),
                        FVEthres = 0.99,
                        thin = TRUE){

  # Perform checks

  if(!is.array(X) || length(dim(X)) != 3) stop('X must be a 3D array')

  n <- dim(X)[1]; M1 <- dim(X)[2]; M2 <- dim(X)[3]

  if(!is.null(tt1)){
    if(!is.vector(tt1, mode = "numeric")) stop('tt1 and tt2 must be numeric vectors if provided')
    if(length(tt1) != M1) stop('Lengths of tt1 and tt2 must match second and third dimensions of X')
    if(any(diff(tt1) <= 0)) stop('values in tt1 and tt2 must be increasing')
  }
  if(!is.null(tt2)){
    if(!is.vector(tt2, mode = "numeric")) stop('tt1 and tt2 must be numeric vectors if provided')
    if(length(tt2) != M2) stop('Lengths of tt1 and tt2 must match second and third dimensions of X')
    if(any(diff(tt2) <= 0)) stop('values in tt1 and tt2 must be increasing')
  }

  if(mnBoot > n) stop('mnBoot cannot be larger than n')

  if(length(JTest) != length(LTest)){
    stop('LTest and JTest must have the same number of elements')
  }

  if(!(all(is.integer(c(JTest, LTest))) && all(c(JTest, LTest) > 0))){
    stop('LTest and JTest must contain only positive integer values')
  }

  if(any(JTest > M1) || any(LTest > M2)){
    stop('Some values in JTest and LTest are too large')
  }

  if(any(diff(JTest) <0 ) || any(diff(LTest) < 0)){
    stop('Value in JTest and LTest must be nondecreasing')
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

  calibList <- c('elBoot', 'quadForm')
  tmpCal <- calibration %in% calibList
  if(!all(tmpCal)){
    if(!any(tmpCal)){
      warning('All provided values in calibration are invalid - resetting to default')
      calibration <- 'elBoot'
    } else {
      warning('Removing invalid elements of calibration')
      calibration <- calibration[tmpCal]
    }
  }
  calibration <- unique(calibration)

  ELctrlBoot <- ELctrl
  ELctrlBoot@verbose = FALSE

  # Get Marginal Basis Expansion

  JTestMax <- max(JTest)
  LTestMax <- max(LTest)
  MBE <- getMBExp(X = X, tt1 = tt1, tt2 = tt2, useFVE = TRUE)
  J <- ncol(MBE$Psi)
  L <- ncol(MBE$Phi)

  if(all(JTest > J) || all(LTest > L)){
    stop('Value in JTest and LTest are all too large as they all involve negative empirical eigenvalues')
  }
  if(JTestMax > J || LTestMax > L){
    warning('At least one element of JTest/LTest is too large after marginal basis expansions computed - removing infeasible values')
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

  # Per-calibration-method output containers. When only one calibration method is
  # requested (the default), bootPvalList/tStatsBootList are unwrapped to their single
  # element before being returned, so existing single-calibration callers see the same
  # flat structure as always; with more than one, they stay as named lists keyed by
  # calibration method (see `calibration` above)

  bootPvalList <- vector("list", length(calibration))
  names(bootPvalList) <- calibration
  for(cal in calibration) bootPvalList[[cal]] <- bootPval

  # Bootstrap results storage
  if(!thin){
    tStatsBootList <- vector("list", length(calibration))
    names(tStatsBootList) <- calibration
    for(cal in calibration){
      tStatsBootList[[cal]] <- vector(mode = "list", length = numTests)
      names(tStatsBootList[[cal]]) <- nullHyp
    }
  }

  # Get necessary EL Fits

  scrsAug <- getELTestData(MBE$scrs, J, L, JTestMax, LTestMax) # get products of scores
  scrsInd <- getScrsInd(JTestMax, LTestMax) # for indexing columns of scrsAug
  ELRes <- getELFits(scrsAug, nullHyp, JTest, LTest, scrsInd, ELctrl, LSmax)
  ELstats <- ELRes$ELstats
  conv <- ELRes$conv

  # For quadForm calibration, precompute fixed full-sample V0 (and, for Sep, U0)
  # estimates once per (J, L); these feed the Theorem 3 quadratic forms used to
  # calibrate every bootstrap replicate below, replacing bootstrap-level EL
  # re-optimization (see `calibration` above)

  if('quadForm' %in% calibration){
    V0list <- vector("list", q)
    U0list <- vector("list", q)
    lamS0list <- vector("list", q)
    NPvec <- integer(q); NWvec <- integer(q); NSvec <- integer(q)

    for(jl in 1:q){
      Jc <- JTest[jl]; Lc <- LTest[jl]
      NPvec[jl] <- Jc*Jc*Lc*(Lc - 1)/2
      NWvec[jl] <- Lc*Jc*(Jc - 1)/2
      NSvec[jl] <- Jc*Lc

      jlIndsSep <- which((pmax(scrsInd[, 1], scrsInd[, 2]) <= Jc) &
                          (pmax(scrsInd[, 3], scrsInd[, 4]) <= Lc))
      scrsAugCur <- scrsAug[, jlIndsSep, drop = FALSE]
      V0list[[jl]] <- getV0Hat(scrsAugCur)

      if('Sep' %in% nullHyp && isTRUE(conv$Sep[jl])){
        thetaSFull <- tail(ELRes$ELoptInfo$Sep[[jl]]$par, NSvec[jl])
        thetaPars <- getThetaFromThetaS(thetaSFull, Jc, Lc)
        U0list[[jl]] <- getU0Hat(V0list[[jl]], thetaPars$theta01, thetaPars$theta02,
                                 NPvec[jl], NWvec[jl], Jc, Lc)$U0
        lamS0list[[jl]] <- thetaSFull
      }
    }
  }

  if('ParSep' %in% nullHyp){ # Get Test Statistic and Run Bootstrap under ParSep Null

    tStats['ParSep', ] <- ELstats['ParSep', ]

    tStatsBootPS <- vapply(1:B, function(b){
      # Set up outputs as if all likelihoods are zero; columns correspond to calibration
      val <- matrix(NA, nrow = q, ncol = length(calibration), dimnames = list(NULL, calibration))

      if(any(conv$ParSep)){ # Run bootstraps for values (J, L) with non-zero likelihood
        convInd <- which(conv$ParSep)
        JTestMaxBoot <- max(JTest[convInd])
        LTestMaxBoot <- max(LTest[convInd])

        Xboot <- drawBoot(MBE, 'ParSep', mnBoot) # Null bootstrap under ParSep - shared across calibration methods
        MBEboot <- getMBExp(Xboot, tt1 = tt1, tt2 = tt2, J = JTestMaxBoot, L = LTestMaxBoot, useFVE = TRUE)
        MBEboot <- alignMBEBoot(MBE, MBEboot, JTest[convInd], LTest[convInd], tt1, tt2)

        scrsAugBoot <- getELTestData(MBEboot$scrs, J = JTestMaxBoot, L = LTestMaxBoot)
        scrsIndBoot <- getScrsInd(JTestMaxBoot, LTestMaxBoot) # for indexing columns of scrsAugBoot

        if('elBoot' %in% calibration){
          ELResBoot <- getELFits(scrsAugBoot, 'ParSep', JTest[convInd], LTest[convInd], scrsIndBoot, ELctrlBoot, LSmax)
          val[convInd, 'elBoot'] <- ELResBoot$ELstats['ParSep', ]
        }
        if('quadForm' %in% calibration){ # no EL re-optimization, just plug hlam* into the fixed V0 quadratic form
          for(jl in convInd){
            Jc <- JTest[jl]; Lc <- LTest[jl]
            jlIndsBoot <- which((pmax(scrsIndBoot[, 1], scrsIndBoot[, 2]) <= Jc) &
                                 (pmax(scrsIndBoot[, 3], scrsIndBoot[, 4]) <= Lc) &
                                 (scrsIndBoot[, 3] < scrsIndBoot[, 4]))
            hlamBootP <- colMeans(scrsAugBoot[, jlIndsBoot, drop = FALSE])
            qf <- getQFQuantities(hlamBootP, mnBoot,
                                  V0list[[jl]][seq_len(NPvec[jl]), seq_len(NPvec[jl]), drop = FALSE],
                                  NPvec[jl], 0L, NSvec[jl])
            val[jl, 'quadForm'] <- qf["qP"]
          }
        }
      }

      return(as.vector(val))
    }, FUN.VALUE = numeric(q * length(calibration)))
    tStatsBootPS <- matrix(tStatsBootPS, nrow = q * length(calibration), ncol = B) # robust reshape even if q*length(calibration) == 1

    for(ci in seq_along(calibration)){
      cal <- calibration[ci]
      rowsCi <- ((ci - 1) * q + 1):(ci * q)
      tStatsBootPSCal <- matrix(tStatsBootPS[rowsCi, , drop = FALSE], nrow = B, ncol = q, byrow = TRUE, dimnames = list(1:B, dimsList))

      bootPvalList[[cal]]['ParSep', ] <- ifelse(conv$ParSep, sapply(1:q, \(jl) (sum(tStatsBootPSCal[, jl] > tStats['ParSep', jl]) + 1) / (B + 1)), 0)
      if(!thin) tStatsBootList[[cal]]$ParSep <- tStatsBootPSCal
    }
  }

  if('WkSep' %in% nullHyp){ # Run Test for Weak Separability

    tStats['WkSep', ] <- ifelse(conv$ParSep, ifelse(conv$WkSep, ELstats['WkSep', ] - ELstats['ParSep', ], Inf), NA)

    tStatsBootWS <- vapply(1:B, function(b){
      # Set up outputs as if all likelihoods are zero; columns correspond to calibration
      val <- matrix(NA, nrow = q, ncol = length(calibration), dimnames = list(NULL, calibration))

      if(any(conv$WkSep)){ # Run bootstraps for values (J, L) with non-zero likelihood
        convInd <- which(conv$WkSep)
        JTestMaxBoot <- max(JTest[convInd])
        LTestMaxBoot <- max(LTest[convInd])

        Xboot <- drawBoot(MBE, 'WkSep', mnBoot) # Null bootstrap under WkSep - shared across calibration methods
        MBEboot <- getMBExp(Xboot, tt1 = tt1, tt2 = tt2, J = JTestMaxBoot, L = LTestMaxBoot, useFVE = TRUE)
        MBEboot <- alignMBEBoot(MBE, MBEboot, JTest[convInd], LTest[convInd], tt1, tt2)

        scrsAugBoot <- getELTestData(MBEboot$scrs, J = JTestMaxBoot, L = LTestMaxBoot)
        scrsIndBoot <- getScrsInd(JTestMaxBoot, LTestMaxBoot) # for indexing columns of scrsAugBoot

        if('elBoot' %in% calibration){
          ELResBoot <- getELFits(scrsAugBoot, 'WkSep', JTest[convInd], LTest[convInd], scrsIndBoot, ELctrlBoot, LSmax)
          val[convInd, 'elBoot'] <- ifelse(ELResBoot$conv$ParSep, ELResBoot$ELstats['WkSep', ] - ELResBoot$ELstats['ParSep', ], Inf)
        }
        if('quadForm' %in% calibration){
          for(jl in convInd){
            Jc <- JTest[jl]; Lc <- LTest[jl]
            jlIndsBoot <- which((pmax(scrsIndBoot[, 1], scrsIndBoot[, 2]) <= Jc) &
                                 (pmax(scrsIndBoot[, 3], scrsIndBoot[, 4]) <= Lc) &
                                 (((scrsIndBoot[, 3] < scrsIndBoot[, 4])) |
                                    ((scrsIndBoot[, 3] == scrsIndBoot[, 4]) & scrsIndBoot[, 1] < scrsIndBoot[, 2])))
            hlamBootPW <- colMeans(scrsAugBoot[, jlIndsBoot, drop = FALSE])
            NPWc <- NPvec[jl] + NWvec[jl]
            qf <- getQFQuantities(hlamBootPW, mnBoot,
                                  V0list[[jl]][seq_len(NPWc), seq_len(NPWc), drop = FALSE],
                                  NPvec[jl], NWvec[jl], NSvec[jl])
            val[jl, 'quadForm'] <- qf["qPW"] - qf["qP"]
          }
        }
      }

      return(as.vector(val))
    }, FUN.VALUE = numeric(q * length(calibration)))
    tStatsBootWS <- matrix(tStatsBootWS, nrow = q * length(calibration), ncol = B)

    for(ci in seq_along(calibration)){
      cal <- calibration[ci]
      rowsCi <- ((ci - 1) * q + 1):(ci * q)
      tStatsBootWSCal <- matrix(tStatsBootWS[rowsCi, , drop = FALSE], nrow = B, ncol = q, byrow = TRUE, dimnames = list(1:B, dimsList))

      bootPvalList[[cal]]['WkSep', ] <- ifelse(conv$WkSep, sapply(1:q, \(jl) (sum(tStatsBootWSCal[, jl] > tStats['WkSep', jl]) + 1) / (B + 1)), 0)
      if(!thin) tStatsBootList[[cal]]$WkSep <- tStatsBootWSCal
    }
  }

  if('Sep' %in% nullHyp){ # Run Test for Separability

    tStats['Sep', ] <- ifelse(conv$WkSep, ifelse(conv$Sep, ELstats['Sep', ] - ELstats['WkSep', ], Inf), NA)

    tStatsBootS <- vapply(1:B, function(b){
      # Set up outputs as if all likelihoods are zero; columns correspond to calibration
      val <- matrix(NA, nrow = q, ncol = length(calibration), dimnames = list(NULL, calibration))

      if(any(conv$Sep)){ # Run bootstraps for values (J, L) with non-zero likelihood
        convInd <- which(conv$Sep)
        JTestMaxBoot <- max(JTest[convInd])
        LTestMaxBoot <- max(LTest[convInd])
        Xboot <- drawBoot(MBE, 'Sep', mnBoot) # Null bootstrap under Sep - shared across calibration methods
        MBEboot <- getMBExp(Xboot, tt1 = tt1, tt2 = tt2, J = JTestMaxBoot, L = LTestMaxBoot, useFVE = TRUE)
        MBEboot <- alignMBEBoot(MBE, MBEboot, JTest[convInd], LTest[convInd], tt1, tt2)

        scrsAugBoot <- getELTestData(MBEboot$scrs, J = JTestMaxBoot, L = LTestMaxBoot)
        scrsIndBoot <- getScrsInd(JTestMaxBoot, LTestMaxBoot) # for indexing columns of scrsAugBoot

        if('elBoot' %in% calibration){
          ELResBoot <- getELFits(scrsAugBoot, 'Sep', JTest[convInd], LTest[convInd], scrsIndBoot, ELctrlBoot, LSmax)
          val[convInd, 'elBoot'] <- ifelse(ELResBoot$conv$WkSep, ELResBoot$ELstats['Sep', ] - ELResBoot$ELstats['WkSep', ], Inf)
        }
        if('quadForm' %in% calibration){ # no nested outer/inner EL optimization at all
          for(jl in convInd){
            Jc <- JTest[jl]; Lc <- LTest[jl]
            jlIndsBoot <- which((pmax(scrsIndBoot[, 1], scrsIndBoot[, 2]) <= Jc) &
                                 (pmax(scrsIndBoot[, 3], scrsIndBoot[, 4]) <= Lc))
            hlamBootFull <- colMeans(scrsAugBoot[, jlIndsBoot, drop = FALSE])
            qf <- getQFQuantities(hlamBootFull, mnBoot, V0list[[jl]],
                                  NPvec[jl], NWvec[jl], NSvec[jl],
                                  U0hat = U0list[[jl]], lamS0 = lamS0list[[jl]])
            val[jl, 'quadForm'] <- qf["qS"] - qf["qPW"]
          }
        }
      }

      return(as.vector(val))
    }, FUN.VALUE = numeric(q * length(calibration)))
    tStatsBootS <- matrix(tStatsBootS, nrow = q * length(calibration), ncol = B)

    for(ci in seq_along(calibration)){
      cal <- calibration[ci]
      rowsCi <- ((ci - 1) * q + 1):(ci * q)
      tStatsBootSCal <- matrix(tStatsBootS[rowsCi, , drop = FALSE], nrow = B, ncol = q, byrow = TRUE, dimnames = list(1:B, dimsList))

      bootPvalList[[cal]]['Sep', ] <- ifelse(conv$Sep, sapply(1:q, \(jl) (sum(tStatsBootSCal[, jl] > tStats['Sep', jl]) + 1) / (B + 1)), 0)
      if(!thin) tStatsBootList[[cal]]$Sep <- tStatsBootSCal
    }
  }

  # Unwrap to the flat (pre-multi-calibration) structure when only one calibration
  # method was requested, so single-calibration callers see an unchanged return shape
  if(length(calibration) == 1){
    bootPvalOut <- bootPvalList[[1]]
    if(!thin) tStatsBootOut <- tStatsBootList[[1]]
  } else {
    bootPvalOut <- bootPvalList
    if(!thin) tStatsBootOut <- tStatsBootList
  }

  res <- list(tStats = tStats, bootPval = bootPvalOut, ELoptInfo = ELRes$ELoptInfo)
  if('quadForm' %in% calibration){
    res$V0 <- V0list
    if('Sep' %in% nullHyp) res$U0 <- U0list
  }
  if(!thin) res$tStatsBoot <- tStatsBootOut
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
    runEL <- TRUE
    for(jl in 1:q){ # fit EL for each (J, L)
      if(runEL){
        jlInds <- which((pmax(scrsInd[, 1], scrsInd[, 2]) <= JTest[jl]) & # j and k no bigger than JTest[jl]
                          (pmax(scrsInd[, 3], scrsInd[, 4]) <= LTest[jl]) & # l and m no bigger than LTest[jl]
                          (scrsInd[, 3] < scrsInd[, 4])) # l < m keeps only values related to PS
        scrsAugCur <- scrsAug[, jlInds] # extract relevant columns
        ELFitPS <- melt::el_mean(scrsAugCur, rep(0, ncol(scrsAugCur)), control = ELctrl)
        conv$ParSep[jl] <- ELFitPS@optim$convergence
        ELstats['ParSep', jl] <- ifelse(conv$ParSep[jl], ELFitPS@statistic, Inf)
        ELoptInfo$ParSep[[jl]] <- ELFitPS@optim
        runEL <- conv$ParSep[jl] # don't run for larger index values if convergence failed
      } else {
        conv$ParSep[jl] <- FALSE
        ELstats['ParSep', jl] <- Inf
        ELoptInfo$ParSep[[jl]] <- NA
      }
    }
  }

  if('WkSep' %in% fitHyp){
    runEL <- TRUE
    if('Sep' %in% fitHyp) logpWS <- vector("list", q) # store logp if needed later for Sep fits
    for(jl in 1:q){ # fit EL for each (J, L)
      if('ParSep' %in% fitHyp) runEL <- runEL && conv$ParSep[jl] # don't run WkSep if ParSep didn't converge
      if(runEL){
        jlInds <- which((pmax(scrsInd[, 1], scrsInd[, 2]) <= JTest[jl]) & # j and k no bigger than JTest[jl]
                         (pmax(scrsInd[, 3], scrsInd[, 4]) <= LTest[jl]) & # l and m no bigger than LTest[jl]
                         (((scrsInd[, 3] < scrsInd[, 4])) | # all (j, k) are valid if l < m
                            ((scrsInd[, 3] == scrsInd[, 4]) & scrsInd[, 1] < scrsInd[, 2]))) # j < k if l = m
        scrsAugCur <- scrsAug[, jlInds] # extract relevant columns
        ELFitWS <- melt::el_mean(scrsAugCur, rep(0, ncol(scrsAugCur)), control = ELctrl)
        conv$WkSep[jl] <- ELFitWS@optim$convergence
        ELstats['WkSep', jl] <- ifelse(conv$WkSep[jl], ELFitWS@statistic, Inf)
        ELoptInfo$WkSep[[jl]] <- ELFitWS@optim
        if('Sep' %in% fitHyp){
          if(conv$WkSep[jl]){
            logpWS[[jl]] <- ELFitWS@logp
          } else {
            logpWS[[jl]] <- NA
          }
        }
        runEL <- conv$WkSep[jl]
      } else {
        conv$WkSep[jl] <- FALSE
        ELstats['WkSep', jl] <- Inf
        ELoptInfo$WkSep[[jl]] <- NA
        if('Sep' %in% fitHyp) logpWS[[jl]] <- NA
      }
    }
  }

  if('Sep' %in% fitHyp){
    runEL <- TRUE
    for(jl in 1:q){ # fit EL for each (J, L)
      runEL <- runEL && conv$WkSep[jl] # don't run Sep if WkSep didn't converge
      if(runEL){
        jlInds <- which((pmax(scrsInd[, 1], scrsInd[, 2]) <= JTest[jl]) & # j and k no bigger than JTest[jl]
                         (pmax(scrsInd[, 3], scrsInd[, 4]) <= LTest[jl])) # l and m no bigger than LTest[jl]
        scrsAugCur <- scrsAug[, jlInds] # extract relevant columns
        SepInds <- (ncol(scrsAugCur) - JTest[jl] * LTest[jl] + 1):ncol(scrsAugCur) # last columns are for Sep

        w <- expm1(logpWS[[jl]]) + 1
        tmp <- matrix(colSums(scrsAugCur[, SepInds] * w), nrow = JTest[jl], ncol = LTest[jl])
        aStrt <- sum(tmp) # trace of estimated Lambda matrix under weak separability
        gammaStrt <- rowSums(tmp)/aStrt
        betaStrt <- colSums(tmp)/aStrt

        ELFitS <- SepELOpt(scrsAugCur, JTest[jl], LTest[jl], gammaStrt, betaStrt, aStrt,
                        mOtr = ELctrl@maxit, mInr = ELctrl@maxit_l,
                        tolOtr = ELctrl@tol, tolInr = ELctrl@tol_l,
                        LSmax = LSmax, verb = ELctrl@verbose)
        conv$Sep[jl] <- ifelse(length(ELFitS@optim$convergence) == 1, # constrained optimization not performed because initial point has likelihood 0
                               FALSE, ELFitS@optim$convergence$FinalInr)
        ELstats['Sep', jl] <- ifelse(conv$Sep[jl], ELFitS@statistic, Inf)
        ELoptInfo$Sep[[jl]] <- ELFitS@optim
        runEL <- conv$Sep[jl]
      } else {
        conv$Sep[jl] <- FALSE
        ELstats['Sep', jl] <- Inf
        ELoptInfo$Sep[[jl]] <- NA
      }
    }
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





