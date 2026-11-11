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

test_that("Stan is selected only for ctsem APIs accepting backend", {
  old <- function(optimize) optimize
  new <- function(optimize, backend = "julia") list(optimize, backend)
  args <- list(optimize = TRUE, vb = FALSE)
  expect_true(.ctmaStanCall(old, args))
  expect_identical(.ctmaStanCall(new, args), list(TRUE, "stan"))
})
