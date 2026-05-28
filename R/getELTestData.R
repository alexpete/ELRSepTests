#' Compute and arrange products of marginal tensor basis scores for use in ELR tests
#'
#' @param scrs n-by-JL matrix of scores in marginal tensor product basis
#' @param J number of first direction basis elements used in the expansion
#' @param L number of second direction basis elements used in the expansion
#' @param JTest integer indicating how many first direction scores to consider
#'              when forming products. Can be at most `J`, with `J` being the default
#' @param LTest integer indicating how many second direction scores to consider
#'              when forming products. Can be at most `L`, with `L` being the default
#'
#'
#' @return matrix of augmented data with columns organized as follows:
#'          cross-products for different direction 2 bases - `LTest * (LTest - 1) * JTest^2 / 2` columns for ParSep
#'          cross-products for different direction 1 bases but same direction 2 basis - `JTest * (JTest - 1) * LTest / 2` columns for WkSep
#'          squares - `LTest * JTest` columns for Sep
#'

getELTestData <- function(scrs, J, L, JTest = J, LTest = L) {

  # Needs to be a matrix
  if (!is.matrix(scrs)) {
    stop("scores must be a matrix")
  }

  # Needs to be numeric
  if (!is.numeric(scrs)) {
    stop("scores must be numeric")
  }

  # Needs to be finite
  if (any(!is.finite(scrs))) {
    stop("scores contain NA, NaN, or Inf values")
  }

  # Making sure that J and L have the right properties
  if (!is.numeric(J) || length(J) != 1 || J <= 0 || J != as.integer(J)) {
    stop("J must be a positive integer scalar")
  }

  if (!is.numeric(L) || length(L) != 1 || L <= 0 || L != as.integer(L)) {
    stop("L must be a positive integer scalar")
  }

  if (!is.numeric(JTest) || length(JTest) != 1 || JTest <= 0 || JTest != as.integer(JTest)) {
    stop("JTest must be a positive integer scalar")
  }

  if (!is.numeric(LTest) || length(LTest) != 1 || LTest <= 0 || LTest != as.integer(LTest)) {
    stop("LTest must be a positive integer scalar")
  }

  if(J*L != ncol(scrs)){
    stop('Provided values for J and L do not correspond to the number of columns in scrs')
  }

  if(JTest > J || LTest > L){
    stop('JTest and LTest can be at most J and L, respectively')
  }

  # Might update later
  if (JTest < 1 || LTest < 1) {
    stop("JTest and LTest must be at least 1")
  }

  # Make sure each score has at least 1 row
  if (nrow(scrs) < 1) {
    stop("scrs must have at least one row")
  }

  n <- nrow(scrs)
  if(JTest < J || LTest < L){ # need to extract only columns related to tested scores
    idx <- sort(sapply(1:LTest, \(l) (l - 1)*J + 1:JTest))
    scrs <- scrs[, idx]
  }

  # Get all distinct products of columns of scrs
  A <- getIndSets(JTest, LTest)

  scrsAug <- scrs[, A$A1, drop = FALSE] * scrs[, A$A2, drop = FALSE]

  return(scrsAug = scrsAug)

}

## Helper function to get correctly ordered indices for products of scores

getIndSets <- function(J, L){

  # column index of S corresponding to pair (j, l)
  q <- function(j, l) j + (l - 1L) * J

  # sizes of the three blocks
  nP <- J * J * L * (L - 1L) / 2L
  nW <- L * J * (J - 1L) / 2L
  nS <- J * L

  A1P <- integer(nP)
  A2P <- integer(nP)
  h <- 1L
  for (l in 1:(L - 1L)) {
    for (m in (l + 1L):L) {
      rng <- h:(h + J * J - 1L)
      A1P[rng] <- q(rep(1:J, each = J), l)   # j varies slow, k varies fast
      A2P[rng] <- q(rep(1:J, times = J), m)
      h <- h + J * J
    }
  }

  A1W <- integer(nW)
  A2W <- integer(nW)
  h <- 1L
  for (l in 1:L) {
    for (j in 1:(J - 1L)) {
      kseq <- (j + 1L):J
      len <- length(kseq)
      rng <- h:(h + len - 1L)
      A1W[rng] <- q(j, l)
      A2W[rng] <- q(kseq, l)
      h <- h + len
    }
  }

  A1S <- integer(nS)
  A2S <- integer(nS)
  h <- 1L
  for (l in 1:L) {
    rng <- h:(h + J - 1L)
    cols <- q(1:J, l)
    A1S[rng] <- cols
    A2S[rng] <- cols
    h <- h + J
  }

  A1 <- c(A1P, A1W, A1S)
  A2 <- c(A2P, A2W, A2S)

  return(list('A1' = A1, 'A2' = A2))
}
