test_that("editable model text reconstructs matrices and constraints", {
  model <- ctsem::ctModel(
    type = "ct", n.latent = 2, n.manifest = 2,
    LAMBDA = diag(2),
    DRIFT = matrix(c("a", "b", "c", "d"), 2, 2)
  )
  model$pars$custom_constraint <- FALSE
  text <- .ctmaModelText(model)
  expect_match(text, "^m <- ctsem::ctModel")
  expect_match(text, "matrix\\(")
  expect_match(text, "nrow = 2, ncol = 2")
  expect_false(grepl("structure(", text, fixed = TRUE))
  expect_false(grepl("[0-9]L\\b", text))
  target <- new.env(parent = globalenv())
  eval(parse(text = text), target)
  expect_equal(target$m$matrices, model$matrices)
  expect_equal(target$m$pars, model$pars)
})
