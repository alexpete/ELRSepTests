# Internal helpers implementing the quadratic-form bootstrap calibration
# (Theorem 3 of the manuscript), used by ELRSepTests(calibration = 'quadForm').
# None of these are exported; see ELRSepTests.R for the public entry point.

#' Sample covariance of score-product vectors (estimate of V0 in Theorem 3)
#'
#' Internal helper. Computes the plug-in estimate of \eqn{V_0 = \mathrm{Var}(\chi_1)}
#' from Theorem 3 of the manuscript, i.e. the sample covariance (with divisor n,
#' not n-1) of the rows of an augmented score-product matrix such as that
#' returned by \code{getELTestData}.
#'
#' @param scrsAugCur n-by-N matrix of score products (columns ordered as
#'                    P-block, then W-block, then S-block, as produced by
#'                    \code{getELTestData}/\code{getIndSets})
#'
#' @return N-by-N estimate of V0
#' @keywords internal

getV0Hat <- function(scrsAugCur) {
  n <- nrow(scrsAugCur)
  Xc <- sweep(scrsAugCur, 2, colMeans(scrsAugCur), FUN = "-")
  crossprod(Xc) / n
}

#' Recover separable parameters from a fitted Kronecker vector
#'
#' Internal helper. \code{SepELOpt}/\code{outer_optimize_rcpp} parameterize the
#' fitted separable value as \code{thetaS = a * kronecker(beta, gamma)}, with
#' both \code{gamma} and \code{beta} constrained to the simplex (they sum to
#' one), and \code{a} an unconstrained positive scale. The manuscript's
#' parameter spaces \eqn{\Theta_1} (unconstrained scale, positive entries) and
#' \eqn{\Theta_2} (simplex) instead correspond to \code{theta1 = a * gamma} and
#' \code{theta2 = beta}. Because both \code{gamma} and \code{beta} are
#' simplicial, \code{theta1} and \code{theta2} can be recovered exactly (no
#' re-optimization, no SVD) from row/column sums of the fitted J-by-L matrix:
#' \code{rowSums} give \code{a * gamma} directly and \code{colSums / sum} give
#' \code{beta} directly.
#'
#' @param thetaS length JL fitted Kronecker vector (e.g. the tail of
#'               \code{ELFitS@optim$par} corresponding to the separable block)
#' @param J number of first-direction indices
#' @param L number of second-direction indices
#'
#' @return list with elements \code{theta01} (= a*gamma, length J, manuscript's
#'         \eqn{\theta_1}) and \code{theta02} (= beta, length L, manuscript's
#'         \eqn{\theta_2}, sums to one)
#' @keywords internal

getThetaFromThetaS <- function(thetaS, J, L) {
  Lam <- matrix(thetaS, nrow = J, ncol = L)
  theta01 <- rowSums(Lam)
  theta02 <- colSums(Lam) / sum(Lam)
  list(theta01 = theta01, theta02 = theta02)
}

#' Analytic Moore-Penrose pseudoinverse for a rank-deficient matrix with known
#' one-dimensional null space
#'
#' Internal helper. If \code{M0} is symmetric positive semi-definite with null
#' space spanned \emph{exactly} by \code{u0}, then, writing \code{un} for the
#' unit vector in the direction of \code{u0},
#' \deqn{M_0^+ = (M_0 + u_nu_n^T)^{-1} - u_nu_n^T.}
#' This holds because \code{M0 + un \%*\% t(un)} is positive definite (the null
#' direction gets eigenvalue 1, all other eigenvalues/eigenvectors of \code{M0}
#' are untouched since \code{un} is orthogonal to them), and inversion and
#' pseudoinversion agree on the eigenspaces they share. Verified against
#' \code{MASS::ginv} to match to numerical precision, and satisfies all four
#' Moore-Penrose conditions.
#'
#' This is used in place of a generic SVD-based pseudoinverse (e.g.
#' \code{MASS::ginv}) because, in context, \code{M0} is estimated and its true
#' null space direction \code{u0 = c(theta02, -theta01)} is known analytically;
#' shifting along this known direction avoids relying on numerically detecting
#' a near-zero singular value of an estimated matrix.
#'
#' @param M0 symmetric positive semi-definite matrix with a one-dimensional
#'           null space
#' @param u0 vector spanning the null space of \code{M0} (need not be
#'           normalized)
#'
#' @return Moore-Penrose pseudoinverse of \code{M0}
#' @keywords internal

getM0Plus <- function(M0, u0) {
  un <- u0 / sqrt(sum(u0^2))
  P <- tcrossprod(un)
  solve(M0 + P) - P
}

#' Tangent-space matrix D0 for the separable parameterization
#'
#' Internal helper implementing eq. (S1.11)-type definition
#' \code{D0 = [I_L %x% theta1 | theta2 %x% I_J]}, an NS-by-(L+J) matrix, where
#' NS = J*L. Its null space (as a map from R^(L+J)) is spanned by
#' \code{c(theta2, -theta1)}, reflecting the scale invariance of the raw
#' Kronecker parameterization.
#'
#' @param theta01 length J vector, manuscript's \eqn{\theta_1} (unconstrained
#'                scale, positive entries)
#' @param theta02 length L vector, manuscript's \eqn{\theta_2} (sums to one)
#' @param J number of first-direction indices
#' @param L number of second-direction indices
#'
#' @return NS-by-(L+J) matrix D0, NS = J*L
#' @keywords internal

getD0 <- function(theta01, theta02, J, L) {
  cbind(kronecker(diag(L), theta01), kronecker(theta02, diag(J)))
}

#' Corrected U0 matrix linking the separable-block EL estimator to the full
#' score-product mean vector
#'
#' Internal helper computing the (corrected) matrix \code{U0} from Theorem 3(c)
#' of the manuscript, satisfying
#' \deqn{\hat\lambda - \hat\lambda_{0,S} = U_0(\hat\lambda - \lambda_0) + o_p(n^{-1/2}),}
#' where \eqn{\hat\lambda_{0,S} = (0, 0, \hat\lambda_0^S)}. The version
#' currently displayed in the supplement, \code{U0 = rbind(G0, C0)}, omits an
#' identity term on the separable block. Direct substitution of Lemma
#' (hlamS_0 - lamS_0 = C0 (hlam - lambda0)) into the definition of
#' \code{hlam - hlam_{0,S}} shows the separable block of U0 must instead be
#' \code{G0S - C0}, where \code{G0S = [0 | I_NS]} selects the separable
#' coordinates. This has been verified numerically (see accompanying tests) to
#' exactly reproduce the direct quadratic-form expansion of
#' \code{2*ellS(0,0,hlamS_0)} (eq. ellSExp in the supplement), which the
#' uncorrected \code{U0 = rbind(G0, C0)} does not.
#'
#' Note: this correction only affects the \emph{separable}-null calibration; it
#' has no effect on \code{QP} or \code{QW}, and does not affect the observed
#' test statistics computed elsewhere in the package (which have always used
#' direct EL optimization via \code{SepELOpt}, not this asymptotic
#' representation).
#'
#' @param V0hat N-by-N estimate of V0 (see \code{getV0Hat}), N = NP+NW+NS,
#'              with rows/columns ordered P-block, W-block, S-block
#' @param theta01 length J vector, manuscript's \eqn{\theta_1} (see
#'                \code{getThetaFromThetaS})
#' @param theta02 length L vector, manuscript's \eqn{\theta_2}
#' @param NP number of partial-separability indices (may be 0)
#' @param NW number of weak-separability indices (may be 0)
#' @param J number of first-direction indices
#' @param L number of second-direction indices
#'
#' @return list with elements \code{U0} (N-by-N, corrected), \code{C0}
#'         (NS-by-N), \code{D0} (NS-by-(L+J)), and \code{M0plus}
#'         ((L+J)-by-(L+J), the analytic pseudoinverse from \code{getM0Plus})
#' @keywords internal

getU0Hat <- function(V0hat, theta01, theta02, NP, NW, J, L) {

  N  <- nrow(V0hat)
  NS <- J * L
  idxS <- (NP + NW + 1):N

  V0inv  <- solve(V0hat)
  V0SS   <- V0inv[idxS, idxS, drop = FALSE]
  V0Srow <- V0inv[idxS, , drop = FALSE]

  D0 <- getD0(theta01, theta02, J, L)
  u0 <- c(theta02, -theta01)

  M0     <- crossprod(D0, V0SS %*% D0)
  M0plus <- getM0Plus(M0, u0)

  C0 <- D0 %*% M0plus %*% t(D0) %*% V0Srow   # NS x N

  G0  <- cbind(diag(NP + NW), matrix(0, NP + NW, NS))
  G0S <- cbind(matrix(0, NS, NP + NW), diag(NS))

  U0 <- rbind(G0, G0S - C0)

  list(U0 = U0, C0 = C0, D0 = D0, M0plus = M0plus)
}

#' Plug-in quadratic-form analogs of the leading quantities in Theorem 3
#'
#' Internal helper. Given a score-product mean vector \code{hlamVec} (either
#' the observed \code{hlam} or a bootstrap \code{hlam*}) and \emph{fixed}
#' full-sample estimates \code{V0hat} (and, if computing the separable-level
#' quantity, \code{U0hat} and the separable center \code{lamS0}), returns
#' \code{nObs} times the quadratic forms from Theorem 3 that stand in for
#' \code{2*ellP(0)}, \code{2*ellW(0,0)}, and \code{2*ellS(0,0,hlamS_0)}.
#'
#' @param hlamVec length N vector (P-block, then W-block, then S-block)
#' @param nObs sample size multiplying the quadratic forms (n for the observed
#'             statistic, or the bootstrap sample size mnBoot for a bootstrap
#'             replicate)
#' @param V0hat N-by-N fixed estimate of V0 from the full sample
#' @param NP number of partial-separability indices (may be 0)
#' @param NW number of weak-separability indices (may be 0)
#' @param NS number of separability indices (= J*L)
#' @param U0hat optional N-by-N corrected U0 matrix (see \code{getU0Hat}); only
#'              needed to compute the separable-level quantity \code{qS}
#' @param lamS0 optional length NS vector giving the separable center (e.g. the
#'              full-sample EL-constrained estimate \code{hlamS_0}); only
#'              needed together with \code{U0hat}
#'
#' @return named numeric vector with elements \code{qP}, \code{qPW}, and
#'         \code{qS} (the latter \code{NA} unless \code{U0hat}/\code{lamS0}
#'         are supplied)
#' @keywords internal

getQFQuantities <- function(hlamVec, nObs, V0hat, NP, NW, NS,
                            U0hat = NULL, lamS0 = NULL) {

  idxP  <- if (NP > 0) seq_len(NP) else integer(0)
  idxW  <- if (NW > 0) (NP + 1):(NP + NW) else integer(0)
  idxPW <- c(idxP, idxW)

  qP <- if (NP > 0) {
    nObs * as.numeric(crossprod(hlamVec[idxP],
                                 solve(V0hat[idxP, idxP, drop = FALSE], hlamVec[idxP])))
  } else {
    0
  }

  qPW <- nObs * as.numeric(crossprod(hlamVec[idxPW],
                                     solve(V0hat[idxPW, idxPW, drop = FALSE], hlamVec[idxPW])))

  qS <- NA_real_
  if (!is.null(U0hat)) {
    lam0 <- c(rep(0, NP + NW), lamS0)
    d  <- hlamVec - lam0
    Ud <- as.numeric(U0hat %*% d)
    qS <- nObs * as.numeric(crossprod(Ud, solve(V0hat, Ud)))
  }

  c(qP = qP, qPW = qPW, qS = qS)
}
