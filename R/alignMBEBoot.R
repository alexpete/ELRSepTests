#' Align a bootstrap marginal basis expansion to the original basis
#'
#' Internal helper for bootstrap calibration in \code{ELRSepTests}. The function
#' aligns the leading bootstrap marginal eigenfunctions to the corresponding
#' original-sample marginal eigenfunctions, allowing for permutation and sign
#' indeterminacy. It then returns a modified copy of \code{MBEboot} containing
#' only the aligned columns needed for the requested test dimensions, with the
#' score matrix permuted and sign-adjusted consistently.
#'
#' Alignment is performed separately for the two marginal bases \code{Psi} and
#' \code{Phi}. For each margin, the function compares the leading \code{K}
#' original eigenfunctions to a small candidate set of leading bootstrap
#' eigenfunctions using absolute inner products. It then chooses the ordered set
#' of \code{K} bootstrap eigenfunctions with maximal total absolute alignment.
#' Finally, signs are chosen so that the aligned bootstrap eigenfunctions have
#' nonnegative inner products with their matched original counterparts.
#'
#' This helper is intentionally dependency-free. It uses exhaustive search over
#' ordered candidate matches, so it is intended for small testing dimensions such
#' as \code{max(JTest)} and \code{max(LTest)} no larger than about 5 or 6. The
#' number of bootstrap candidate eigenfunctions considered is controlled by
#' \code{buffer}; if \code{buffer = NULL}, the candidate set size is
#' \code{min(ncol(Boot), max(K + 2L, 2L * K))} for each margin.
#'
#' @param MBE List returned by \code{getMBExp} for the original sample. Must
#'   contain elements \code{Psi}, \code{Phi}, and \code{scrs}.
#' @param MBEboot List returned by \code{getMBExp} for a bootstrap sample. Must
#'   contain elements \code{Psi}, \code{Phi}, and \code{scrs}.
#' @param JTest Integer or integer vector giving the first-margin testing
#'   dimensions. The function keeps \code{max(JTest)} aligned columns of
#'   \code{Psi}.
#' @param LTest Integer or integer vector giving the second-margin testing
#'   dimensions. The function keeps \code{max(LTest)} aligned columns of
#'   \code{Phi}.
#' @param buffer Optional nonnegative integer giving the number of extra
#'   bootstrap eigenfunctions, beyond \code{K}, to consider as alignment
#'   candidates in each margin. If \code{NULL}, a default candidate-set size of
#'   \code{max(K + 2L, 2L * K)} is used, truncated at the number of available
#'   bootstrap eigenfunctions.
#'
#' @return A modified copy of \code{MBEboot} with aligned and truncated
#'   \code{Psi}, \code{Phi}, and \code{scrs}. The returned object also contains
#'   an \code{alignInfo} element with the selected permutations, signs,
#'   alignments, and candidate-set sizes for each margin.
#'
#' @keywords internal
alignMBEBoot <- function(MBE, MBEboot, JTest, LTest, buffer = NULL) {
  J <- max(JTest)
  L <- max(LTest)

  if (!is.matrix(MBE$Psi) || !is.matrix(MBEboot$Psi)) {
    stop("MBE$Psi and MBEboot$Psi must be matrices.")
  }
  if (!is.matrix(MBE$Phi) || !is.matrix(MBEboot$Phi)) {
    stop("MBE$Phi and MBEboot$Phi must be matrices.")
  }
  if (!is.matrix(MBEboot$scrs)) {
    stop("MBEboot$scrs must be a matrix.")
  }
  if (nrow(MBE$Psi) != nrow(MBEboot$Psi)) {
    stop("MBE$Psi and MBEboot$Psi must have the same number of rows.")
  }
  if (nrow(MBE$Phi) != nrow(MBEboot$Phi)) {
    stop("MBE$Phi and MBEboot$Phi must have the same number of rows.")
  }
  if (!is.null(buffer)) {
    if (length(buffer) != 1L || buffer < 0 || buffer != as.integer(buffer)) {
      stop("buffer must be NULL or a single nonnegative integer.")
    }
    buffer <- as.integer(buffer)
  }

  allPerms <- function(x, k) {
    if (k == 1L) {
      return(matrix(x, ncol = 1L))
    }

    do.call(rbind, lapply(seq_along(x), function(i) {
      cbind(x[i], allPerms(x[-i], k - 1L))
    }))
  }

  alignOne <- function(Boot, Orig, K, buffer) {
    if (K > ncol(Orig)) {
      stop("Requested alignment dimension exceeds number of original eigenfunctions.")
    }
    if (K > ncol(Boot)) {
      stop("Requested alignment dimension exceeds number of bootstrap eigenfunctions.")
    }

    if (is.null(buffer)) {
      nCand <- min(ncol(Boot), max(K + 2L, 2L * K))
    } else {
      nCand <- min(ncol(Boot), K + buffer)
    }

    BootCand <- Boot[, seq_len(nCand), drop = FALSE]
    OrigK <- Orig[, seq_len(K), drop = FALSE]

    ## Alignment scores: rows index original eigenfunctions and columns index
    ## candidate bootstrap eigenfunctions.
    A <- abs(crossprod(OrigK, BootCand))

    ## Exhaustively search ordered selections of K distinct bootstrap candidate
    ## eigenfunctions. This is feasible because K is small in the testing use case.
    perms <- allPerms(seq_len(nCand), K)
    scores <- apply(perms, 1L, function(p) {
      sum(A[cbind(seq_len(K), p)])
    })
    perm <- perms[which.max(scores), ]

    Aligned <- BootCand[, perm, drop = FALSE]

    signs <- sign(diag(crossprod(OrigK, Aligned)))
    signs[signs == 0] <- 1
    Aligned <- sweep(Aligned, 2L, signs, `*`)

    list(
      eigfuns = Aligned,
      perm = perm,
      signs = signs,
      alignment = diag(crossprod(OrigK, Aligned)),
      nCand = nCand
    )
  }

  aPsi <- alignOne(MBEboot$Psi, MBE$Psi, J, buffer)
  aPhi <- alignOne(MBEboot$Phi, MBE$Phi, L, buffer)

  qBoot <- ncol(MBEboot$Psi)

  if (ncol(MBEboot$scrs) < qBoot * ncol(MBEboot$Phi)) {
    stop("MBEboot$scrs has fewer columns than expected from MBEboot$Psi and MBEboot$Phi.")
  }

  ## Score columns are assumed to be ordered with the first-margin index moving
  ## fastest: column (j, l) is (l - 1) * qBoot + j.
  scrsNew <- matrix(NA_real_, nrow = nrow(MBEboot$scrs), ncol = J * L)

  ctr <- 1L
  for (l in seq_len(L)) {
    for (j in seq_len(J)) {
      oldCol <- (aPhi$perm[l] - 1L) * qBoot + aPsi$perm[j]
      scrsNew[, ctr] <- MBEboot$scrs[, oldCol] *
        aPsi$signs[j] * aPhi$signs[l]
      ctr <- ctr + 1L
    }
  }

  colnames(scrsNew) <- NULL

  MBEbootAligned <- MBEboot
  MBEbootAligned$Psi <- aPsi$eigfuns
  MBEbootAligned$Phi <- aPhi$eigfuns
  MBEbootAligned$scrs <- scrsNew
  MBEbootAligned$alignInfo <- list(
    Psi = aPsi,
    Phi = aPhi
  )

  MBEbootAligned
}
