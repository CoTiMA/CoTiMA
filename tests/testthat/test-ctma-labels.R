test_that("custom drift and diffusion specifications retain labels and fixed values", {
  labels <- ctmaLabels(
    n.latent = 2,
    drift = c("V1toV1", 0, -0.15, "V2toV2"),
    diff = c("V1withV1", 0, 0, "V2withV2")
  )

  expect_identical(
    labels$driftParams,
    c("V1toV1", "0", "-0.15", "V2toV2")
  )
  expect_identical(
    labels$diffParams,
    c("V1withV1", "0", "0", "V2withV2")
  )
  expect_identical(labels$diffNames, c("V1withV1", "V2withV2"))
  expect_false(exists("diffParamsParams", inherits = FALSE))
})

test_that("custom drift and diffusion dimensions are checked", {
  expect_error(
    ctmaLabels(n.latent = 2, drift = c("a", "b", "c")),
    "exactly 4 entries",
    fixed = TRUE
  )
  expect_error(
    ctmaLabels(n.latent = 2, diff = c("a", "b", "c")),
    "exactly 4 entries",
    fixed = TRUE
  )
})
