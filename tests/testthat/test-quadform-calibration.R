test_that("getM0Plus matches SVD-based pseudoinverse and satisfies MP conditions", {
  skip_if_not_installed("MASS")
  set.seed(101)

  J <- 3; L <- 4
  theta01 <- abs(rnorm(J)); theta01 <- theta01 / sum(theta01)
  theta02 <- abs(rnorm(L)); theta02 <- theta02 / sum(theta02)

  D0 <- ELRSepTests:::getD0(theta01, theta02, J, L)
  u0 <- c(theta02, -theta01)
  expect_lt(max(abs(D0 %*% u0)), 1e-10)

  NS <- J * L
  Bmat <- matrix(rnorm(NS * NS), NS, NS)
  V0SS <- Bmat %*% t(Bmat) + NS * diag(NS)

  M0 <- crossprod(D0, V0SS %*% D0)
  Mp_analytic <- ELRSepTests:::getM0Plus(M0, u0)
  Mp_svd <- MASS::ginv(M0)

  expect_equal(Mp_analytic, Mp_svd, tolerance = 1e-8)

  # Moore-Penrose conditions
  expect_equal(M0 %*% Mp_analytic %*% M0, M0, tolerance = 1e-6)
  expect_equal(Mp_analytic %*% M0 %*% Mp_analytic, Mp_analytic, tolerance = 1e-6)
  expect_equal(M0 %*% Mp_analytic, t(M0 %*% Mp_analytic), tolerance = 1e-6)
  expect_equal(Mp_analytic %*% M0, t(Mp_analytic %*% M0), tolerance = 1e-6)
})

test_that("corrected U0 reproduces the direct quadratic-form expansion of 2*ellS(0,0,hlamS_0)", {
  set.seed(202)

  NP <- 3; NW <- 2; J <- 2; L <- 3
  NS <- J * L
  N <- NP + NW + NS

  Amat <- matrix(rnorm(N * N), N, N)
  V0 <- Amat %*% t(Amat) + N * diag(N)

  theta01 <- abs(rnorm(J)); theta01 <- theta01 / sum(theta01)
  theta02 <- abs(rnorm(L)); theta02 <- theta02 / sum(theta02)
  lamS0 <- kronecker(theta02, theta01)  # matches thetaS = a*kronecker(beta,gamma) with a folded into theta01

  U0info <- ELRSepTests:::getU0Hat(V0, theta01, theta02, NP, NW, J, L)
  U0 <- U0info$U0
  C0 <- U0info$C0

  V0inv <- solve(V0)

  # Direct construction: pick a random d = hlam - lambda0, define
  # hlamS_0 - lamS_0 := C0 %*% d (per Lemma S.9), then hlam - hlam_{0,S}
  # by its own definition, and compare the two quadratic forms.
  d <- rnorm(N)
  deltaS0 <- as.numeric(C0 %*% d)
  idxP <- 1:NP; idxW <- (NP + 1):(NP + NW); idxS <- (NP + NW + 1):N
  diff_direct <- c(d[idxP], d[idxW], d[idxS] - deltaS0)

  Q_direct <- as.numeric(t(diff_direct) %*% V0inv %*% diff_direct)
  Q_via_U0 <- as.numeric(t(d) %*% t(U0) %*% V0inv %*% U0 %*% d)

  expect_equal(Q_via_U0, Q_direct, tolerance = 1e-8)

  # And confirm the *uncorrected* U0 = rbind(G0, C0) does NOT match in general
  G0  <- cbind(diag(NP + NW), matrix(0, NP + NW, NS))
  U0_wrong <- rbind(G0, C0)
  Q_via_U0_wrong <- as.numeric(t(d) %*% t(U0_wrong) %*% V0inv %*% U0_wrong %*% d)
  expect_false(isTRUE(all.equal(Q_via_U0_wrong, Q_direct, tolerance = 1e-8)))
})

test_that("getThetaFromThetaS exactly recovers theta01 (=a*gamma) and theta02 (=beta)", {
  set.seed(303)
  J <- 4; L <- 3
  gamma <- abs(rnorm(J)); gamma <- gamma / sum(gamma)
  beta  <- abs(rnorm(L)); beta  <- beta / sum(beta)
  a <- 2.7

  thetaS <- a * kronecker(beta, gamma)
  out <- ELRSepTests:::getThetaFromThetaS(thetaS, J, L)

  expect_equal(out$theta01, a * gamma, tolerance = 1e-10)
  expect_equal(out$theta02, beta, tolerance = 1e-10)
})

test_that("ELRSepTests(calibration = 'quadForm') runs end-to-end and returns valid p-values", {
  skip_on_cran()
  set.seed(404)

  n <- 40; M1 <- 6; M2 <- 8
  X <- array(rnorm(n * M1 * M2), dim = c(n, M1, M2))

  res <- ELRSepTests(X, JTest = 2L, LTest = 2L, nullHyp = c('ParSep', 'WkSep', 'Sep'),
                     B = 30L, calibration = 'quadForm', thin = TRUE)

  expect_true(all(is.na(res$bootPval) | (res$bootPval >= 0 & res$bootPval <= 1)))
  expect_true(is.list(res$V0))
  expect_true(is.list(res$U0))
})

test_that("ELRSepTests observed tStats agree exactly between calibration = 'elBoot' and 'quadForm'", {
  skip_on_cran()
  set.seed(505)

  n <- 40; M1 <- 6; M2 <- 8
  X <- array(rnorm(n * M1 * M2), dim = c(n, M1, M2))

  res_el <- ELRSepTests(X, JTest = 2L, LTest = 2L, B = 10L, calibration = 'elBoot')
  res_qf <- ELRSepTests(X, JTest = 2L, LTest = 2L, B = 10L, calibration = 'quadForm')

  expect_equal(res_el$tStats, res_qf$tStats)
  expect_null(res_el$V0)
  expect_true(is.list(res_qf$V0))
})

test_that("calibration = c('elBoot','quadForm') matches single-calibration runs exactly (same bootstrap draw) and shape-shifts output", {
  skip_on_cran()

  n <- 40; M1 <- 6; M2 <- 8

  set.seed(606)
  X <- array(rnorm(n * M1 * M2), dim = c(n, M1, M2))

  set.seed(707)
  res_both <- ELRSepTests(X, JTest = 2L, LTest = 2L, B = 15L,
                          calibration = c('elBoot', 'quadForm'), thin = FALSE)

  set.seed(707)
  res_el_only <- ELRSepTests(X, JTest = 2L, LTest = 2L, B = 15L,
                             calibration = 'elBoot', thin = FALSE)

  set.seed(707)
  res_qf_only <- ELRSepTests(X, JTest = 2L, LTest = 2L, B = 15L,
                             calibration = 'quadForm', thin = FALSE)

  # observed statistics unaffected by calibration
  expect_equal(res_both$tStats, res_el_only$tStats)
  expect_equal(res_both$tStats, res_qf_only$tStats)

  # output shape-shifts to a named list keyed by calibration method
  expect_true(is.list(res_both$bootPval) && setequal(names(res_both$bootPval), c('elBoot', 'quadForm')))
  expect_true(is.list(res_both$tStatsBoot) && setequal(names(res_both$tStatsBoot), c('elBoot', 'quadForm')))

  # sharing the same bootstrap draw means each component matches its single-calibration
  # counterpart exactly, given the same seed
  expect_equal(res_both$bootPval$elBoot, res_el_only$bootPval)
  expect_equal(res_both$bootPval$quadForm, res_qf_only$bootPval)
  expect_equal(res_both$tStatsBoot$elBoot, res_el_only$tStatsBoot)
  expect_equal(res_both$tStatsBoot$quadForm, res_qf_only$tStatsBoot)

  # single-calibration runs keep the flat (non-list) structure
  expect_true(is.matrix(res_el_only$bootPval))
  expect_true(is.matrix(res_qf_only$bootPval))
})
