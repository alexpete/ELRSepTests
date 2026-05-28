test_that("getELTestData runs on valid input", {

  scrs <- matrix(1:12, nrow = 2)  # 2 rows, 6 cols
  J <- 2
  L <- 3

  out <- getELTestData(scrs, J, L)

  expect_true(is.matrix(out))
  expect_true(nrow(out) == 2)
  expect_true(ncol(out) > 0)
})

# Fix this
test_that("output dimensions are consistent with theory", {

  scrs <- matrix(1:27, nrow = 3)
  J <- 3
  L <- 3

  out <- getELTestData(scrs, J, L)

  expected_cols <- (L*(L-1)/2)*J^2 +
    L*(J*(J-1)/2) +
    L*J

  expect_equal(ncol(out), expected_cols)
})

# Where the JTests and LTests are 2
test_that("output dimensions are consistent with theory", {

  scrs <- matrix(1:27, nrow = 3)
  J <- 3
  L <- 3
  JTest <- 2
  LTest <- 2

  out <- getELTestData(scrs, J, L, JTest = 2, LTest = 2)

  expected_cols <- (LTest*(LTest-1)/2)*JTest^2 +
    LTest*(JTest*(JTest-1)/2) +
    LTest*JTest

  expect_equal(ncol(out), expected_cols)
})

# Wrong Matrix shape
test_that("errors when J*L does not match ncol(scrs)", {

  scrs <- matrix(1:10, nrow = 2)
  J <- 2
  L <- 3  # 2*3 = 6 ≠ 5 cols

  expect_error(getELTestData(scrs, J, L))
})

# Non-numeric input
test_that("errors when scrs is not numeric", {

  scrs <- matrix(letters[1:6], nrow = 2)
  J <- 2
  L <- 3

  expect_error(getELTestData(scrs, J, L))
})

# NA/Inf handling
test_that("errors when scrs contains NA", {

  scrs <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2)
  J <- 2
  L <- 3

  expect_error(getELTestData(scrs, J, L))
})

# JTest/LTest too large
test_that("errors when JTest or LTest exceed limits", {

  scrs <- matrix(1:12, nrow = 2)
  J <- 2
  L <- 3

  expect_error(getELTestData(scrs, J, L, JTest = 5))
})

# Invalid J/L types
test_that("errors when J is not integer", {

  scrs <- matrix(1:12, nrow = 2)

  expect_error(getELTestData(scrs, 2.5, 3))
})

test_that("output is deterministic and stable", {

  scrs <- matrix(1:12, nrow = 2)
  J <- 2
  L <- 3

  out1 <- getELTestData(scrs, J, L)
  out2 <- getELTestData(scrs, J, L)

  expect_equal(out1, out2)
})


























