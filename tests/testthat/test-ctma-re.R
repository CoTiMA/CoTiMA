test_that("CoTiMA-RE covariance aliases are normalized", {
  expect_identical(CoTiMA:::.ctma_covariance_name("ure"), "independent")
  expect_identical(CoTiMA:::.ctma_covariance_name("ere"), "equicorrelated")
  expect_identical(CoTiMA:::.ctma_covariance_name("equal"), "equicorrelated")
  expect_identical(CoTiMA:::.ctma_covariance_name("cre"), "full")
  expect_error(CoTiMA:::.ctma_covariance_name("unknown"),
    "Unknown covariance structure")
})

test_that("automatic matrix blocks count correlation parameters correctly", {
  spec <- data.frame(
    re = 1:8,
    matrix = c(rep("DRIFT", 4), rep("DIFFUSION", 3), "T0VAR"),
    row = c(1, 2, 1, 2, 1, 2, 2, 2),
    col = c(1, 1, 2, 2, 1, 1, 2, 1),
    parameter = paste0("p", 1:8),
    stringsAsFactors = FALSE
  )

  result <- CoTiMA:::.ctma_re_blocks(
    spec = spec,
    covariance = "block",
    covarianceBlocks = list(
      DRIFT = "equicorrelated",
      DIFFUSION = "independent",
      T0VAR = "full"
    ),
    nlatent = 2,
    verbose = FALSE
  )

  expect_identical(result$blocks$block,
    c("DRIFT", "DIFFUSION", "T0VAR"))
  expect_identical(result$blocks$size, c(4L, 3L, 1L))
  expect_identical(result$blocks$structure,
    c("equicorrelated", "independent", "independent"))
  code <- CoTiMA:::.ctma_stan_block_code(result$blocks)
  expect_true(code$has_corr)
  expect_identical(code$ncorr, 1L)
  expect_match(code$par, "studyre_rho_b1", fixed = TRUE)
})

test_that("custom blocks require exactly the selected cells", {
  spec <- data.frame(
    re = 1L, matrix = "DRIFT", row = 1L, col = 1L,
    parameter = "drift_1_1", stringsAsFactors = FALSE
  )
  assignment <- matrix(NA_character_, 2, 2)
  assignment[1, 1] <- "dynamic"
  assignment[2, 2] <- "not_selected"

  expect_error(
    CoTiMA:::.ctma_re_blocks(
      spec = spec,
      covariance = "block",
      reBlocks = list(DRIFT = assignment),
      blockStructures = c(dynamic = "independent"),
      nlatent = 2,
      verbose = FALSE
    ),
    "not selected"
  )
})

test_that("MAP and conditional ML set opposite prior switches", {
  standata <- list(priors = 1L, studyre_use_hyperpriors = 1L)
  map_data <- CoTiMA:::.ctma_optimizer_data(standata, "map")
  ml_data <- CoTiMA:::.ctma_optimizer_data(standata, "ml")

  expect_identical(map_data$priors, 1L)
  expect_identical(map_data$studyre_use_hyperpriors, 1L)
  expect_identical(ml_data$priors, 0L)
  expect_identical(ml_data$studyre_use_hyperpriors, 0L)
  expect_error(CoTiMA:::.ctma_optimizer_data(list(), "map"),
    "expected 'priors' switch", fixed = TRUE)
})

test_that("optimizer output receives a leading draw dimension", {
  result <- CoTiMA:::.ctma_add_draw_dimension(list(
    scalar = 2,
    vector = c(1, 2),
    matrix = matrix(1:4, 2, 2)
  ))

  expect_identical(dim(result$scalar), c(1L, 1L))
  expect_identical(dim(result$vector), c(1L, 2L))
  expect_identical(dim(result$matrix), c(1L, 2L, 2L))
})

