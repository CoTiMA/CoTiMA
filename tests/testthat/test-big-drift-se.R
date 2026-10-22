test_that("fixed covariance SEs do not prevent drift analyses", {
  errors <- matrix(.1, nrow = 3, ncol = 12)
  colnames(errors) <- c("V1toV1", "V2toV1", "V1toV2", "V2toV2",
    paste0("covariance", 1:8))
  errors[, 11] <- 0
  expect_true(.ctmaBiGCheckDriftSE(errors, 2))
  errors[2, 3] <- 0
  expect_error(.ctmaBiGCheckDriftSE(errors, 2), "DRIFT.*V1toV2")
  expect_true(.ctmaBiGCheckDriftSE(errors, 2, excludeStudies = 2))
  for (value in c(NA_real_, Inf, -.1)) {
    errors[2, 3] <- value
    expect_error(.ctmaBiGCheckDriftSE(errors, 2), "strictly positive DRIFT")
  }
})
