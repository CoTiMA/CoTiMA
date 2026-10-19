make_par_map_fit <- function() {
  matsetup <- matrix(
    c(1, 1, 1, 1, 0, 0),
    nrow = 1L,
    dimnames = list(
      "DRIFT",
      c("row", "col", "param", "transform", "indvarying", "tipred")
    )
  )
  data <- list(
    matsetup = matsetup,
    TIPREDEFFECTsetup = matrix(1L, 1L, 1L),
    nparams = 1L,
    nindvarying = 0L,
    nindvaryingoffdiagonals = 0L,
    intoverpop = 1L,
    nsubjects = 5L,
    ntipredeffects = 1L,
    tipredsdata = matrix(c(0, 1), ncol = 1L,
      dimnames = list(NULL, "TI1")),
    choleskymats = 0L
  )
  list(
    ctstanmodelbase = list(
      pars = data.frame(
        matrix = "DRIFT", row = 1L, col = 1L,
        param = "V1toV1", transform = "param", indvarying = FALSE
      ),
      TIpredNames = "TI1"
    ),
    data = data,
    stanfit = list(rawposterior = cbind(
      rawpopmeans = c(1, 2),
      tipredeffectparams = c(1, 1)
    ))
  )
}

test_that("ctmaParMap accepts a direct ctsem fit", {
  map <- ctmaParMap(make_par_map_fit())
  expect_equal(map$rawpopmeans$rawpopmeans_col, 1L)
  expect_equal(map$rawpopmeans$rawposterior_col, 1L)
  expect_equal(map$tipredeffectparams$effect_label, "TI1_TO_V1toV1")
  expect_equal(map$tipredeffectparams$rawposterior_col, 2L)
  expect_identical(attr(map, "ctsem_fit_source"), "ctsem fit")
})

test_that("ctmaReMoments reconstructs studies and uses equal-weight moments", {
  object <- structure(
    list(studyFitList = make_par_map_fit()),
    class = "CoTiMAFit"
  )
  result <- suppressMessages(ctmaReMoments(
    object,
    include_matrices = FALSE,
    include_individual = FALSE
  ))
  expect_s3_class(result, "CoTiMAReMoments")
  expect_s3_class(result, "CoTiMAFit")
  expect_equal(result$transformedStudyDraws[, 1L, 1L], c(1, 2))
  expect_equal(result$transformedStudyDraws[, 1L, 2L], c(2, 3))
  population <- subset(result$distributionParameters,
    distribution_parameter == "mean")
  heterogeneity <- subset(result$distributionParameters,
    distribution_parameter == "between_study_sd")
  expect_equal(population$mean, 2)
  expect_equal(heterogeneity$mean, sqrt(.5))
  expect_equal(summary(result)$numberOfStudies, 2L)
})

test_that("ctmaReMoments supports inverse-variance study weighting", {
  fit <- make_par_map_fit()
  fit$stanfit$rawposterior[, 2L] <- c(0, 2)
  object <- structure(list(studyFitList = fit), class = "CoTiMAFit")
  result <- ctmaReMoments(
    object,
    SEweighting = TRUE,
    include_matrices = FALSE,
    include_individual = FALSE
  )
  expect_identical(result$weighting, "inverse_variance")
  expect_equal(unname(colSums(result$studyWeights)), 1)
  expect_gt(result$studyWeights[1L, 1L], result$studyWeights[2L, 1L])
  expect_lt(result$effectiveNumberOfStudies[1L], 2)
  population <- subset(result$distributionParameters,
    distribution_parameter == "mean")
  expect_equal(population$mean, 1.6)
  expect_identical(result$covarianceWeighting, "equal study weights")
})

test_that("ctmaParMap extracts one pooled fit from a CoTiMA object", {
  object <- structure(
    list(studyFitList = make_par_map_fit()),
    class = "CoTiMAFit"
  )
  map <- ctmaParMap(object)
  expect_equal(map$rawpopmeans$parameter, "V1toV1")
  expect_identical(attr(map, "ctsem_fit_source"), "CoTiMA fit")
})

test_that("ctmaParMap rejects ambiguous separate-study CoTiMA fits", {
  object <- structure(
    list(studyFitList = list(make_par_map_fit(), make_par_map_fit())),
    class = "CoTiMAFit"
  )
  expect_error(
    ctmaParMap(object),
    "contains 2 separately fitted ctsem models",
    fixed = TRUE
  )
})
