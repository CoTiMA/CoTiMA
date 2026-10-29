test_that("removed ctsem vb argument is handled compatibly", {
  args <- .ctmaSanitizeCtFitArgs(list(
    optimize = FALSE,
    priors = FALSE,
    vb = FALSE
  ))

  expect_false("vb" %in% names(args))
  expect_identical(args$optimize, FALSE)
  expect_identical(args$priors, FALSE)
  expect_error(
    .ctmaSanitizeCtFitArgs(list(vb = TRUE)),
    "Variational Bayes is no longer supported by ctsem"
  )
})
