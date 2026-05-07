test_that("getTrapzVec() computes correct weights", {
  expect_equal(getTrapzVec(c(0,1,2,3)), c(0.5,1,1,0.5))
})

test_that("getTrapzVec() handles uneven spacing", {
  expect_equal(getTrapzVec(c(0,1,3,4)), c(0.5,1.5,1.5,0.5))
})

test_that("getTrapzVec() integrates constant function", {
  expect_equal(sum(getTrapzVec(c(0,1,2,3))), 3)
})

test_that("errors if tt is a matrix", {
  tt <- matrix(c(0, 1, 2, 3), ncol = 2)

  expect_error(getTrapzVec(tt))
})

test_that("errors if tt is a data frame", {
  tt <- data.frame(x = c(0, 1, 2, 3))

  expect_error(getTrapzVec(tt))
})

test_that("errors if tt is a scalar", {
  expect_error(getTrapzVec(5))
})

test_that("errors if tt is character", {
  expect_error(getTrapzVec(c("a", "b", "c")))
})

test_that("errors if tt contains NA", {
  expect_error(getTrapzVec(c(0, 1, NA, 2)))
})

test_that("errors if tt is not increasing", {
  expect_error(getTrapzVec(c(0, 2, 1, 3)))
})

test_that("errors if tt contains Inf", {
  expect_error(getTrapzVec(c(0, 1, 2, Inf)))
})

test_that("errors if tt contains -Inf", {
  expect_error(getTrapzVec(c(-Inf, 0, 1, 2)))
})

