test_that("ctmaOTL supports arbitrary drift dimensions and exact grid indices", {
  drift <- matrix(c(
    -.7, .2, 0,
    -.1, -.5, .1,
    .05, 0, -.4
  ), 3, 3, byrow = TRUE)
  grid <- seq(0, 5, by = .25)
  result <- ctmaOTL(driftMat = drift, timeRange = grid)

  expected_values <- vapply(
    grid,
    function(time) OpenMx::expm(drift * time)[1, 2],
    numeric(1)
  )
  expect_equal(result$optimalCrossLag[1, 2],
               grid[which.max(abs(expected_values))])
  expect_equal(dim(result$optimalCrossLag), c(3L, 3L))
})

test_that("ctmaOTL produces draw-wise uncertainty intervals", {
  drift_draws <- list(
    matrix(c(-.8, .20, -.10, -.6), 2, 2, byrow = TRUE),
    matrix(c(-.7, .25, -.08, -.5), 2, 2, byrow = TRUE),
    matrix(c(-.9, .15, -.12, -.7), 2, 2, byrow = TRUE)
  )
  result <- ctmaOTL(
    driftMat = drift_draws,
    timeRange = seq(0, 5, by = .25),
    CI = TRUE,
    confidence = .90
  )

  expect_equal(result$nDriftDraws, 3L)
  expect_equal(dim(result$optimalCrossLagCI), c(2L, 2L, 3L))
  expect_equal(dimnames(result$optimalCrossLagCI)[[3]],
               c("lower", "median", "upper"))
  expect_equal(result$confidence, .90)
})

test_that("ctmaOTL needs drift uncertainty to calculate intervals", {
  drift <- matrix(c(-.8, .2, -.1, -.6), 2, 2, byrow = TRUE)
  expect_error(
    ctmaOTL(driftMat = drift, timeRange = 0:5, CI = TRUE),
    "multiple drift matrices"
  )
})
