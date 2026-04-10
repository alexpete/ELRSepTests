#' Outer-layer optimization for empirical likelihood estimates under separability constraint
#'
#' Implements outer-layer optimization (over space of separabe parameters)
#'
#' @param scrsAug augmented marginal basis scores for computing estimating equations
#' @param JTest number of first direction basis elements being tested
#' @param LTest number of second direction basis elements being tested
#' @param gammaStrt JTest starting values for gamma, simplicial vector of normalized first direction eigenvalues
#' @param betaStrt LTest starting values for beta, simplicial vector of normalized second direction eigenvalues
#' @param aStrt starting values for a, the trace of Lambda
#' @param mOtr maximum number of outer iterations (default = 200L)
#' @param mInr maximum number of inner iterations (default = 25L)
#' @param tolOtr convergence tolerance for outer layer (default = 1e-06)
#' @param tolInr convergence tolerance for inner layer (default = 1e-06)
#' @param LSmax maximal number of step size reductions before terminating line search (default = 100L)
#' @param verb logical indicator for printing status (default = FALSE)
#'
#' @return Object of class EL from package 'melt'
#' @export

SepELOpt <- function(scrsAug, JTest, LTest, gammaStrt, betaStrt, aStrt,
                     mOtr = 200L, mInr = 25L, tolOtr = 1e-06, tolInr = 1e-06,
                     LSmax = 100L, verb = FALSE){

  # Initialize parameters
  gammaCur <- gammaStrt; betaCur <- betaStrt; aCur <- aStrt

  # Initial evaluation of EL
  thetaSCur <- aCur*kronecker(betaCur, gammaCur)
  numParsS <- length(thetaSCur)
  elCur <- melt::el_mean(scrsAug, par = c(rep(0, ncol(scrsAug) - numParsS), thetaSCur),
                   control = melt::el_control(maxit_l = mInr, tol_l = tolInr))
  if(!elCur@optim$convergence){ # starting point is bad, likelihood should be zero
    if(verb){
      warning('EL evaluation did not converge at provided starting values.  Consider increasing mInr')
    }
    return(elCur)
  } else {
    SepInd <- (ncol(scrsAug) - numParsS + 1):ncol(scrsAug)

    res <- tryCatch(outer_optimize_rcpp(
        scrsAug = scrsAug,
        betaInit = betaStrt,            # your starting beta
        gammaInit = gammaStrt,          # your starting gamma
        aInit = aStrt,                  # starting a
        LTest = LTest,
        JTest = JTest,
        numParsS = numParsS,
        SepInd = SepInd,               # integer( ) indices (1-based)
        mOtr = mOtr,
        LSmax = LSmax,
        mInr = mInr,
        tolInr = tolInr,
        tolOtr = tolOtr,
        verb = verb
      ),
      error = function(e){
        stop('At some point, the empirical likelihood could not be evaluated.  Consider reducing JTest and/or LTest to simplify the required optimization.')
      }
    )

    elNew <- res$elNew
    LSFail <- res$LSFail
    convGrad <- res$convGrad; convStep <- res$convStep; maxGrad <- res$maxGrad
    niter <- res$niter
    elNew@optim$par <- rep(0, ncol(scrsAug))
    elNew@optim$par[SepInd] <- res$thetaS
    elNew@optim$iterations <- niter # number of outer iterations
    elNew@optim$convergence <- list('LSFail' = LSFail,
                                    'Obj' = res$convObj,
                                    'Grad' = convGrad,
                                    'Step' = convStep,
                                    'maxGrad' = maxGrad,
                                    'FinalInr' = elNew@optim$convergence,  # did final el_eval converge?
                                    'objTrace' = res$objTrace,
                                    'gradTrace' = res$gradTrace)
    elNew@df <- as.integer(LTest + JTest - 1)

    return(elNew)
  }
}
