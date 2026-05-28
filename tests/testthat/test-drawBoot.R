library(testthat)

# 1. MBExp must be a list
test_that("MBExp must be a list", {

  expect_error(drawBoot("not a list", "ParSep", 3))
})


# 2. Required components exist
test_that("MBExp must contain required components", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4)
    # Lambda missing
  )

  expect_error(drawBoot(MBExp, "ParSep", 3))
})


# 3. Non-finite checks
test_that("MBExp rejects NA/NaN/Inf values", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(4)
  )

  MBExp$Psi[1,1] <- NA
  expect_error(drawBoot(MBExp, "ParSep", 3))

  MBExp$Psi[1,1] <- 1
  MBExp$Phi[1,1] <- Inf
  expect_error(drawBoot(MBExp, "ParSep", 3))

  MBExp$Phi[1,1] <- 1
  MBExp$scrs[1,1] <- NaN
  expect_error(drawBoot(MBExp, "ParSep", 3))
})


# 4. Dimension checks (Psi / Phi)
test_that("Psi and Phi must have valid dimensions", {

  MBExp <- list(
    Psi = matrix(numeric(0), 0, 0),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(4)
    )

  expect_error(drawBoot(MBExp, "ParSep", 3))
  })


# 5. scrs dimension consistency
test_that("scrs must match Psi x Phi structure", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(25), 5, 5),
    Lambda = diag(4)
  )

  expect_error(drawBoot(MBExp, "ParSep", 3))
})


# 6. mnBoot checks (integer, finite, bounds)
test_that("mnBoot validation works", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(4)
  )

  expect_error(drawBoot(MBExp, "ParSep", 2.5))  # fractional
  expect_error(drawBoot(MBExp, "ParSep", "3"))  # character
  expect_error(drawBoot(MBExp, "ParSep", c(3,4))) # vector
  expect_error(drawBoot(MBExp, "ParSep", -1))    # negative
  expect_error(drawBoot(MBExp, "ParSep", 0))     # zero
  expect_error(drawBoot(MBExp, "ParSep", 1000))  # too large
  expect_error(drawBoot(MBExp, "ParSep", NA))  # NA value
})


# 7. nullType checks
test_that("nullType must be valid", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(4)
  )

  expect_error(drawBoot(MBExp, 123, 3))
  expect_error(drawBoot(MBExp, c("ParSep","WkSep"), 3))
  expect_error(drawBoot(MBExp, "invalid", 3))
})


# 8. Lambda checks
test_that("Lambda must be matrix with positive diagonal", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(4)
  )

  MBExp$Lambda <- as.data.frame(MBExp$Lambda)
  expect_error(drawBoot(MBExp, "ParSep", 3))
})

test_that("Lambda diagonal must be strictly positive", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(c(1, 0, 1, 1))
  )

  expect_error(drawBoot(MBExp, "ParSep", 3))
})

# 9. Output shape tests
test_that("ParSep returns correct array shape", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(4)
  )

  out <- drawBoot(MBExp, "ParSep", 3)

  expect_equal(dim(out), c(3, 3, 4))
})

test_that("WkSep returns correct array shape", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(4)
  )

  out <- drawBoot(MBExp, "WkSep", 3)

  expect_equal(dim(out), c(3, 3, 4))
})

test_that("Sep returns correct array shape", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(4)
  )

  out <- drawBoot(MBExp, "Sep", 3)

  expect_equal(dim(out), c(3, 3, 4))
})

# 10. Basic output validity
test_that("output is numeric and finite array", {

  MBExp <- list(
    Psi = matrix(rnorm(6), 3, 2),
    Phi = matrix(rnorm(8), 4, 2),
    scrs = matrix(rnorm(20), 5, 4),
    Lambda = diag(4)
  )

  out <- drawBoot(MBExp, "ParSep", 3)

  expect_true(is.array(out))
  expect_true(is.numeric(out))
  expect_true(all(is.finite(out)))
})



