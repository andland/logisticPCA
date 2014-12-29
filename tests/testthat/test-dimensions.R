context("Dimensions")

# construct a low rank matrix in the logit scale
rows = 100
cols = 10
k = 1
set.seed(1)
mat_logit = outer(rnorm(rows), rnorm(cols))

mat = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit)) * 1.0

lpca = logisticPCA(mat, M = 4, k = k, main_effects = FALSE)
lsvd = logisticSVD(mat, k = k, main_effects = FALSE, max_iters = 50)

pred1 = predict(lpca, mat)
fit1a = fitted(lpca, type = "link")
fit1b = fitted(lpca, type = "response")

pred2 = predict(lsvd, mat)
fit2a = fitted(lsvd, type = "link")
fit2b = fitted(lsvd, type = "response")

test_that("k = 1 LPCA", {
  expect_equal(dim(lpca$U), c(cols, k))
  expect_equal(dim(lpca$PCs), c(rows, k))
  expect_equal(length(lpca$mu), cols)
  
  expect_equal(dim(pred1), c(rows, k))
  expect_equal(dim(fit1a), c(rows, cols))
  expect_equal(dim(fit1b), c(rows, cols))
})

test_that("k = 1 LSVD", {
  expect_equal(dim(lsvd$B), c(cols, k))
  expect_equal(dim(lsvd$A), c(rows, k))
  expect_equal(length(lsvd$mu), cols)
  
  expect_equal(dim(pred2), c(rows, k))
  expect_equal(dim(fit2a), c(rows, cols))
  expect_equal(dim(fit2b), c(rows, cols))
})


rm(lsvd, lpca, pred1, pred2, fit1a, fit1b, fit2a, fit2b)

k = 2
lpca = logisticPCA(mat, M = 4, k = k, main_effects = FALSE)
lsvd = logisticSVD(mat, k = k, main_effects = FALSE, max_iters = 50)

pred1 = predict(lpca, mat)
fit1a = fitted(lpca, type = "link")
fit1b = fitted(lpca, type = "response")

pred2 = predict(lsvd, mat)
fit2a = fitted(lsvd, type = "link")
fit2b = fitted(lsvd, type = "response")

test_that("k = 2 LPCA", {
  expect_equal(dim(lpca$U), c(cols, k))
  expect_equal(dim(lpca$PCs), c(rows, k))
  expect_equal(length(lpca$mu), cols)
  
  expect_equal(dim(pred1), c(rows, k))
  expect_equal(dim(fit1a), c(rows, cols))
  expect_equal(dim(fit1b), c(rows, cols))
})

test_that("k = 2 LSVD", {
  expect_equal(dim(lsvd$B), c(cols, k))
  expect_equal(dim(lsvd$A), c(rows, k))
  expect_equal(length(lsvd$mu), cols)
  
  expect_equal(dim(pred2), c(rows, k))
  expect_equal(dim(fit2a), c(rows, cols))
  expect_equal(dim(fit2b), c(rows, cols))
})
