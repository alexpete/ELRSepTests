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
  SepInd <- (ncol(scrsAug) - numParsS + 1):ncol(scrsAug)

  makeSepELReturn <- function(el, thetaS, niter = 0L, LSFail = NA,
                              converged = FALSE, convObj = FALSE,
                              convGrad = FALSE, convStep = FALSE,
                              maxGrad = NA_real_, objTrace = NULL,
                              gradTrace = NULL, nEvalFail = 0L,
                              nInnerFail = 0L, outerAttempted = FALSE,
                              status = "outer_not_run", message = NULL){
    finalInr <- el@optim$convergence

    if(is.null(objTrace)){
      objTrace <- el@logl / nrow(scrsAug)
    }
    if(is.null(gradTrace)){
      gradTrace <- NA_real_
    }

    el@optim$par <- rep(0, ncol(scrsAug))
    el@optim$par[SepInd] <- thetaS
    el@optim$iterations <- as.integer(niter)
    el@optim$convergence <- list(
      LSFail = LSFail,
      Obj = convObj,
      Grad = convGrad,
      Step = convStep,
      maxGrad = maxGrad,
      FinalInr = finalInr,
      objTrace = objTrace,
      gradTrace = gradTrace,
      converged = converged,
      outerAttempted = outerAttempted,
      status = status,
      message = message,
      nEvalFail = as.integer(nEvalFail),
      nInnerFail = as.integer(nInnerFail)
    )
    el@df <- as.integer(LTest + JTest - 1)
    el
  }

  elCur <- melt::el_mean(scrsAug,
                         par = c(rep(0, ncol(scrsAug) - numParsS), thetaSCur),
                         control = melt::el_control(maxit_l = mInr, tol_l = tolInr))

  if(!elCur@optim$convergence){
    if(verb){
      warning('EL evaluation did not converge at provided starting values. Outer optimization was not run. Consider increasing mInr.',
              call. = FALSE)
    }

    return(makeSepELReturn(
      el = elCur,
      thetaS = thetaSCur,
      niter = 0L,
      LSFail = NA,
      converged = FALSE,
      convObj = FALSE,
      convGrad = FALSE,
      convStep = FALSE,
      maxGrad = NA_real_,
      objTrace = elCur@logl / nrow(scrsAug),
      gradTrace = NA_real_,
      nEvalFail = 0L,
      nInnerFail = 0L,
      outerAttempted = FALSE,
      status = "initial_inner_failed",
      message = "Initial inner EL optimization failed; outer optimization was not run."
    ))
  }

  res <- outer_optimize_rcpp(
    scrsAug = scrsAug,
    betaInit = betaStrt,
    gammaInit = gammaStrt,
    aInit = aStrt,
    LTest = LTest,
    JTest = JTest,
    numParsS = numParsS,
    SepInd = SepInd,
    mOtr = mOtr,
    LSmax = LSmax,
    mInr = mInr,
    tolInr = tolInr,
    tolOtr = tolOtr,
    verb = verb
  )

  status <- if(isTRUE(res$converged)){
    "converged"
  } else if(isTRUE(res$LSFail)){
    "outer_line_search_failed"
  } else if(res$niter >= mOtr){
    "outer_iteration_limit"
  } else {
    "outer_not_converged"
  }

  message <- switch(status,
                    converged = "Outer optimization converged.",
                    outer_line_search_failed = "Outer optimization stopped after line search failed; returning best available iterate.",
                    outer_iteration_limit = "Outer optimization reached the maximum number of iterations; returning best available iterate.",
                    outer_not_converged = "Outer optimization stopped without satisfying convergence criteria; returning best available iterate.")

  elNew <- makeSepELReturn(
    el = res$elNew,
    thetaS = res$thetaS,
    niter = res$niter,
    LSFail = res$LSFail,
    converged = res$converged,
    convObj = res$convObj,
    convGrad = res$convGrad,
    convStep = res$convStep,
    maxGrad = res$maxGrad,
    objTrace = res$objTrace,
    gradTrace = res$gradTrace,
    nEvalFail = res$nEvalFail,
    nInnerFail = res$nInnerFail,
    outerAttempted = TRUE,
    status = status,
    message = message
  )

  if(!isTRUE(res$converged)){
    warning(message, call. = FALSE)
  }
  if((res$nEvalFail > 0L || res$nInnerFail > 0L) && verb){
    warning("Line search rejected candidate steps: ",
            res$nEvalFail, " failed EL evaluations and ",
            res$nInnerFail, " inner convergence failures.",
            call. = FALSE)
  }

  return(elNew)
}
