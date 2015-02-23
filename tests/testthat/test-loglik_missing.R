context("Log Likelihood Calculation")

rows = 100
cols = 10
set.seed(2015)
mat = matrix(sample(c(0, 1), rows * cols, TRUE), rows, cols)
q = 2 * mat - 1
mis = matrix(sample(c(FALSE, TRUE), rows * cols, TRUE), rows, cols)
theta = matrix(rnorm(rows * cols), rows, cols)

test_that("no missing data", {
  expect_equal(log_like_Bernoulli(x = mat, theta = theta),
               sum(log(inv.logit.mat(q * theta))))
  expect_equal(log_like_Bernoulli(q = q, theta = theta),
               sum(log(inv.logit.mat(q * theta))))
})

is.na(mat[mis]) <- TRUE
q = 2 * mat - 1
q[mis] = 0

test_that("missing data", {
  expect_equal(log_like_Bernoulli(x = mat, theta = theta),
               sum(log(inv.logit.mat(q * theta))[q!=0]))
  expect_equal(log_like_Bernoulli(q = q, theta = theta),
               sum(log(inv.logit.mat(q * theta))[q!=0]))
})
