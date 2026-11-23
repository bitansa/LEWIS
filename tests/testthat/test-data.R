test_that("the data set tcga_prostate has the documented form", {
  X <- LEWIS::tcga_prostate
  expect_true(is.matrix(X))
  expect_true(is.double(X))
  expect_identical(dim(X), c(462L, 8L))
  expect_identical(colnames(X), c("AR", "FOXA1", "HOXB13", "KLK3", "KLK2",
                                  "TMPRSS2", "NKX3-1", "FKBP5"))
  expect_identical(anyDuplicated(rownames(X)), 0L)
  expect_true(all(startsWith(rownames(X), "TCGA-")))
  expect_true(all(is.finite(X)))
  expect_true(all(abs(colMeans(X)) < 1e-8))
  expect_true(all(abs(apply(X, 2L, sd) - 1) < 1e-8))
})

test_that("lewis() runs on tcga_prostate", {
  fit <- lewis(LEWIS::tcga_prostate, n_iter = 1000, center = FALSE, ess_min = 0,
               seed = 1)
  expect_s3_class(fit, "lewis")
  expect_identical(colnames(fit$D_hat), colnames(LEWIS::tcga_prostate))
})
