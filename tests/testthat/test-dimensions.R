context("Dimensions")

# construct a low rank matrix in the logit scale
rows = 100
cols = 10
k = 1
set.seed(1)
mat_logit = outer(rnorm(rows), rnorm(cols))

mat = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit)) * 1.0

lpca = logisticPCA(mat, k = k, m = 4, main_effects = FALSE)
lsvd = logisticSVD(mat, k = k, main_effects = FALSE, conv_criteria = 1e-4, partial_decomp = FALSE)
clpca = convexLogisticPCA(mat, m = 4, k = k, main_effects = FALSE)

lpca_part = logisticPCA(mat, k = k, m = 4, main_effects = FALSE, partial_decomp = TRUE)
lsvd_part = logisticSVD(mat, k = k, main_effects = FALSE, conv_criteria = 1e-4, partial_decomp = TRUE)
clpca_part = convexLogisticPCA(mat, m = 4, k = k, main_effects = FALSE, partial_decomp = TRUE)

pred1 = predict(lpca, mat)
pred1l = predict(lpca, mat, type = "link")
pred1r = predict(lpca, mat, type = "response")
fit1l = fitted(lpca, type = "link")
fit1r = fitted(lpca, type = "response")

pred2 = predict(lsvd, mat)
pred2l = predict(lsvd, mat, type = "link")
pred2r = predict(lsvd, mat, type = "response")
fit2l = fitted(lsvd, type = "link")
fit2r = fitted(lsvd, type = "response")

pred3 = predict(clpca, mat)
pred3l = predict(clpca, mat, type = "link")
pred3r = predict(clpca, mat, type = "response")

test_that("correct classes", {
  expect_is(lpca, "lpca")
  expect_is(lsvd, "lsvd")
  expect_is(clpca, "clpca")

  expect_is(pred1, "matrix")
  expect_is(pred1l, "matrix")
  expect_is(pred1r, "matrix")
  expect_is(fit1l, "matrix")
  expect_is(fit1r, "matrix")

  expect_is(pred2, "matrix")
  expect_is(pred2l, "matrix")
  expect_is(pred2r, "matrix")
  expect_is(fit2l, "matrix")
  expect_is(fit2r, "matrix")

  expect_is(pred3, "matrix")
  expect_is(pred3l, "matrix")
  expect_is(pred3r, "matrix")
})

test_that("k = 1 LPCA", {
  expect_equal(dim(lpca$U), c(cols, 1))
  expect_equal(dim(lpca$PCs), c(rows, 1))
  expect_equal(length(lpca$mu), cols)

  expect_equal(dim(pred1), c(rows, 1))
  expect_equal(dim(pred1l), c(rows, cols))
  expect_equal(dim(pred1r), c(rows, cols))
  expect_equal(dim(fit1l), c(rows, cols))
  expect_equal(dim(fit1r), c(rows, cols))
})

test_that("k = 1 LSVD", {
  expect_equal(dim(lsvd$B), c(cols, 1))
  expect_equal(dim(lsvd$A), c(rows, 1))
  expect_equal(length(lsvd$mu), cols)

  expect_equal(dim(pred2), c(rows, 1))
  expect_equal(dim(pred2l), c(rows, cols))
  expect_equal(dim(pred2r), c(rows, cols))
  expect_equal(dim(fit2l), c(rows, cols))
  expect_equal(dim(fit2r), c(rows, cols))
})

test_that("k = 1 CLPCA", {
  expect_equal(dim(clpca$U), c(cols, 1))
  expect_equal(dim(clpca$H), c(cols, cols))
  expect_equal(dim(clpca$PCs), c(rows, 1))
  expect_equal(length(clpca$mu), cols)

  expect_equal(dim(pred3), c(rows, 1))
  expect_equal(dim(pred3l), c(rows, cols))
  expect_equal(dim(pred3r), c(rows, cols))
})


rm(lsvd, lpca, clpca, pred1, pred1l, pred1r, pred2, pred2l, pred2r,
   pred3, pred3l, pred3r, fit1l, fit1r, fit2l, fit2r)

k = 2
lpca = logisticPCA(mat, m = 4, k = k, main_effects = FALSE)
lsvd = logisticSVD(mat, k = k, main_effects = FALSE, conv_criteria = 1e-4, partial_decomp = FALSE)
clpca = convexLogisticPCA(mat, m = 4, k = k, main_effects = FALSE)

pred1 = predict(lpca, mat)
pred1l = predict(lpca, mat, type = "link")
pred1r = predict(lpca, mat, type = "response")
fit1l = fitted(lpca, type = "link")
fit1r = fitted(lpca, type = "response")

pred2 = predict(lsvd, mat)
pred2l = predict(lsvd, mat, type = "link")
pred2r = predict(lsvd, mat, type = "response")
fit2l = fitted(lsvd, type = "link")
fit2r = fitted(lsvd, type = "response")

pred3 = predict(clpca, mat)
pred3l = predict(clpca, mat, type = "link")
pred3r = predict(clpca, mat, type = "response")

test_that("k = 2 LPCA", {
  expect_equal(dim(lpca$U), c(cols, 2))
  expect_equal(dim(lpca$PCs), c(rows, 2))
  expect_equal(length(lpca$mu), cols)

  expect_equal(dim(pred1), c(rows, 2))
  expect_equal(dim(pred1l), c(rows, cols))
  expect_equal(dim(pred1r), c(rows, cols))
  expect_equal(dim(fit1l), c(rows, cols))
  expect_equal(dim(fit1r), c(rows, cols))
})

test_that("k = 2 LSVD", {
  expect_equal(dim(lsvd$B), c(cols, 2))
  expect_equal(dim(lsvd$A), c(rows, 2))
  expect_equal(length(lsvd$mu), cols)

  expect_equal(dim(pred2), c(rows, 2))
  expect_equal(dim(pred2l), c(rows, cols))
  expect_equal(dim(pred2r), c(rows, cols))
  expect_equal(dim(fit2l), c(rows, cols))
  expect_equal(dim(fit2r), c(rows, cols))
})

test_that("k = 2 CLPCA", {
  expect_equal(dim(clpca$U), c(cols, 2))
  expect_equal(dim(clpca$H), c(cols, cols))
  expect_equal(dim(clpca$PCs), c(rows, 2))
  expect_equal(length(clpca$mu), cols)

  expect_equal(dim(pred3), c(rows, 2))
  expect_equal(dim(pred3l), c(rows, cols))
  expect_equal(dim(pred3r), c(rows, cols))
})

is_less_than_equal <- function (expected, label = NULL, ...)
{
  if (is.null(label)) {
    label <- testthat:::find_expr("expected")
  }
  else if (!is.character(label) || length(label) != 1) {
    label <- deparse(label)
  }
  function(actual) {
    diff <- expected - actual
    expectation(diff >= 0, paste0("not less than ", label,
                                  ". Difference: ", format(diff)), paste0("is less than ",
                                                                          label))
  }
}

test_that("response between 0 and 1", {
  expect_gte(min(pred1r), 0)
  expect_gte(min(fit1r), 0)
  
  expect_lte(max(pred1r), 1)
  expect_lte(max(fit1r), 1)
  
  expect_gte(min(pred2r), 0)
  expect_gte(min(fit2r), 0)
  
  expect_lte(max(pred2r), 1)
  expect_lte(max(fit2r), 1)
  
  expect_gte(min(pred3r), 0)
  expect_lte(max(pred3r), 1)
})
