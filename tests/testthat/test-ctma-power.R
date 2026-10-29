.power_fixture <- function(n = 2L, waves = 2L, start = 0, rho = 0) {
  a <- diag(-seq(.3, .5, length.out = n), n)
  a[1, 2] <- .15; a[2, 1] <- -.08
  if (n > 2) { a[2, 3] <- .03; a[3, 1] <- -.04 }
  initial <- matrix(rho, n, n); diag(initial) <- 1
  q <- -(a %*% initial + initial %*% t(a))
  table <- data.frame(matrix = rep(c("DRIFT", "DIFFUSIONcov", "T0cov"), each = n^2),
    row = rep(rep(seq_len(n), n), 3), col = rep(rep(seq_len(n), each = n), 3),
    Mean = c(a, q, initial), Sd = .05, check.names = FALSE)
  labels <- as.vector(outer(seq_len(n), seq_len(n), function(i, j) paste0("V", j, "toV", i)))
  list(init = list(n.latent = n, n.studies = 2, n.manifest = n, activeDirectory = tempdir(),
      ctModel = list(manifesttype = rep(0, n)),
      parameterNames = list(DRIFT = labels),
      statisticsList = list(allSampleSizes = c(100, 200), allTpoints = rep(waves, 2),
        allDeltas = rep(c(1, 2) + start, each = waves - 1)),
      studyList = list(list(originalStudyNo = 1), list(originalStudyNo = 2))),
    pooled = list(summary = list(parmatrices = table),
      modelResults = list(DRIFT = a, DIFFUSION = q, T0VAR = initial)))
}

.power_design <- function(x, range = c(0, 2, 1), probabilities = c(.8, .95))
  .ctmaPowerDesign(x$init, range, probabilities, 100, .05, tempdir(), TRUE, NULL, 4)

test_that("direct Gaussian regression matches lavaan for asymmetric 2D and 3D systems", {
  for (n in c(2L, 3L)) for (rho in c(0, .3)) {
    x <- .power_fixture(n, rho = rho)
    design <- .power_design(x)
    model <- .ctmaPowerRegressionModel(n)
    expect_true(model$unrestricted)
    for (lag in c(.2, 1, 5)) {
      correlation <- .ctmaPowerImplied(.ctmaPowerMatrices(x$pooled, n), lag)
      direct <- .ctmaPowerDirect(correlation, design$targets, 250)
      fitted <- .ctmaPowerLavaan(correlation, model, design$targets, 250, tempdir(), character(), character(), 1, TRUE)
      expect_equal(direct$estimate, fitted$estimate, tolerance = 1e-5)
      expect_equal(direct$se, fitted$se, tolerance = 1e-5)
      expect_equal(direct$pvalue, fitted$pvalue, tolerance = 1e-5)
      expect_lt(max(abs(.ctmaPowerF2(direct) - .ctmaPowerF2(fitted))), 1e-7)
    }
  }
})

test_that("two-sided probabilities have the correct null limit and match both t tails", {
  for (n in c(10, 100, 1000)) for (f2 in c(0, .01, .1)) {
    critical <- stats::qt(.975, n - 3)
    expected <- stats::pt(-critical, n - 3, ncp = sqrt(n * f2)) +
      stats::pt(critical, n - 3, ncp = sqrt(n * f2), lower.tail = FALSE)
    expect_equal(.ctmaPowerProbability(n, f2, 2), expected, tolerance = 1e-8)
  }
  expect_equal(.ctmaPowerProbability(100, 0, 2), .05)
  expect_equal(.ctmaPowerProbability(100, 0, 2, useMBESS = TRUE), .025)
  expect_true(is.na(.ctmaPowerProbability(3, .1, 2)))
})

test_that("integer sample-size search handles zero, strong and weak effects", {
  expect_identical(.ctmaPowerRequiredN(0, 2, .8), Inf)
  expect_equal(.ctmaPowerRequiredN(100, 2, .8), 4)
  for (f2 in c(1e-8, .0001, .01, .25, 1)) {
    N <- .ctmaPowerRequiredN(f2, 2, .8)
    expect_gte(.ctmaPowerProbability(N, f2, 2), .8)
    if (N > 4) expect_lt(.ctmaPowerProbability(N - 1, f2, 2), .8)
  }
  for (f2 in c(.01, .04, .25)) expect_equal(.ctmaPowerRequiredN(f2, 2, .8, TRUE),
    MBESS::ss.power.reg.coef(Cohen.f2 = f2, p = 2, desired.power = .8)$Necessary.Sample.Size)
})

test_that("matrix extraction uses coordinates, lower triangles and original time units", {
  x <- .power_fixture(3L)
  original <- .ctmaPowerMatrices(x$pooled, 3L)
  tab <- x$pooled$summary$parmatrices
  tab <- tab[rev(seq_len(nrow(tab))), ]
  tab <- tab[tab$matrix == "DRIFT" | tab$row >= tab$col, ]
  x$pooled$summary$parmatrices <- tab
  x$pooled$modelResults$DRIFT[,] <- 999
  expect_equal(.ctmaPowerMatrices(x$pooled, 3L)$DRIFT, original$DRIFT)
  expect_equal(.ctmaPowerMatrices(x$pooled, 3L)$DIFFUSION, original$DIFFUSION)
  x$pooled$CoTiMAStanctArgs <- list(scaleTime = 1/12)
  x$pooled$summary$parmatrices$Mean[x$pooled$summary$parmatrices$matrix != "T0cov"] <-
    x$pooled$summary$parmatrices$Mean[x$pooled$summary$parmatrices$matrix != "T0cov"] * 12
  restored <- .ctmaPowerMatrices(x$pooled, 3L)
  expect_equal(.ctmaPowerImplied(restored, 2), .ctmaPowerImplied(original, 2))
})

test_that("implied covariances support singular drift and keep 3D dimensions", {
  matrices <- list(DRIFT = matrix(0, 3, 3), DIFFUSION = diag(3), T0VAR = diag(3))
  implied <- .ctmaPowerImplied(matrices, 2)
  expect_equal(dim(implied), c(6L, 6L))
  expect_equal(unname(diag(implied[1:3, 4:6])), rep(1/sqrt(3), 3))
})

test_that("eligibility detects restrictions, indicator types and missing metadata", {
  x <- .power_fixture()
  full <- .ctmaPowerRegressionModel(2)
  expect_true(full$unrestricted)
  restricted <- .ctmaPowerRegressionModel(2, c(full$syntax, "V1T1 ~~ 0*V2T1"))
  expect_false(restricted$unrestricted)
  fixed <- .ctmaPowerRegressionModel(2, c(full$syntax, "V1T1 ~ 0*V2T0"))
  expect_false(fixed$unrestricted)
  expect_error(.ctmaPowerRegressionModel(2, c(full$syntax, "V1T1 ~ common*V1T0", "V2T1 ~ common*V2T0")), "Shared|constraints")
  expect_true(.ctmaPowerGaussian(x$init, x$pooled)$verified)
  x$init$ctModel$manifesttype <- c(0, 1)
  expect_error(.ctmaPowerGaussian(x$init, x$pooled), "binary or ordinal")
  x$init$ctModel <- NULL
  expect_false(.ctmaPowerGaussian(x$init, x$pooled)$verified)
  x <- .power_fixture(); x$init$n.manifest <- 3
  expect_error(.ctmaPowerGaussian(x$init, x$pooled), "directly observed")
  x <- .power_fixture(); x$init$ctModel$MANIFESTVAR <- diag(.1, 2)
  expect_error(.ctmaPowerGaussian(x$init, x$pooled), "measurement-error")
  x <- .power_fixture(); x$init$ctModel$LAMBDA <- diag(2, 2)
  expect_error(.ctmaPowerGaussian(x$init, x$pooled), "measurement loadings")
})

test_that("complete power analysis preserves two-wave data, bounds, labels and signed effects", {
  x <- .power_fixture(3L, 2L, 5)
  design <- .power_design(x, c(5, 7, 1))
  expect_equal(design$times, 5:7)
  analysis <- suppressMessages(.ctmaPowerAnalyse(x$init, x$pooled, design, FALSE, "auto", NULL, character(), character(), TRUE, 4))
  expect_identical(analysis$route, "covariance")
  expect_true(grepl("No lavaan", analysis$note))
  expect_equal(rownames(analysis$required)[1:3], as.character(5:7))
  expect_length(analysis$posthoc, 6L)
  expect_equal(dim(analysis$posthoc[[1]]), c(4L, 4L))
  expect_true(all(is.finite(analysis$posthoc[[1]][1:2, 3])))
  for (i in seq_along(analysis$details)) {
    result <- analysis$details[[i]]
    expect_equal(result$pvalue, 2 * stats::pnorm(abs(result$estimate/result$se), lower.tail = FALSE))
    expect_true(any(result$estimate < 0))
    expect_equal(result$estimate, unname(analysis$coefficients[i, ]))
  }
})

test_that("fallback is explicit and fixed targets are not assigned fabricated power", {
  x <- .power_fixture(); design <- .power_design(x, c(1, 1, 1), .8)
  syntax <- c(.ctmaPowerRegressionModel(2)$syntax, "V1T1 ~ 0*V2T0")
  expect_message(a <- .ctmaPowerAnalyse(x$init, x$pooled, design, FALSE, "auto", syntax, character(), character(), TRUE, 4), "Lavaan regression fits used")
  expect_identical(a$route, "lavaan")
  fixed <- which(design$targets$parameter == "V2toV1")
  expect_true(is.na(a$details[[1]]$f2[fixed]))
  expect_true(is.na(a$required[1, fixed]))
  expect_error(.ctmaPowerAnalyse(x$init, x$pooled, design, FALSE, "covariance", syntax, character(), character(), TRUE, 4), "unavailable")
})

test_that("minimum and tied-lag reporting selects an actual minimizing lag", {
  tied <- .ctmaPowerOptima(c(100, 200, 100), c(1, 2, 3))
  expect_equal(tied$first, 1)
  expect_equal(tied$all, c(1, 3))
  expect_equal(tied$minimum, 100)
  expect_true(is.na(.ctmaPowerOptima(c(NA, NA), 1:2)$first))
  # Real zero-effect system gives Inf, not arbitrary optimal lags or warnings.
  x <- .power_fixture(); x$pooled$summary$parmatrices$Mean[x$pooled$summary$parmatrices$matrix == "DRIFT" &
    x$pooled$summary$parmatrices$row != x$pooled$summary$parmatrices$col] <- 0
  x$pooled$summary$parmatrices$Mean[x$pooled$summary$parmatrices$matrix == "DIFFUSIONcov" &
    x$pooled$summary$parmatrices$row != x$pooled$summary$parmatrices$col] <- 0
  a <- suppressMessages(.ctmaPowerAnalyse(x$init, x$pooled, .power_design(x), FALSE, "auto", NULL, character(), character(), TRUE, 4))
  expect_true(all(is.infinite(a$required["Min N", 1:4])))
  expect_true(all(is.na(a$required["Opt. Lag", 1:4])))
  expect_true(all(a$posthoc[[1]][1:2, 3] == .05))
})

test_that("reduced-fit caches are reused across power levels and validate the model", {
  x <- .power_fixture(); directory <- withr::local_tempdir()
  design <- .power_design(x, c(1, 1, 1), c(.8, .95))
  design$directory <- directory
  saved <- suppressMessages(.ctmaPowerAnalyse(x$init, x$pooled, design, FALSE,
    "auto", NULL, character(), "reduced", FALSE, 4))
  expect_identical(saved$route, "lavaan")
  expect_length(list.files(directory, pattern = "rds$"), 2)
  loaded <- suppressMessages(.ctmaPowerAnalyse(x$init, x$pooled, design, FALSE,
    "auto", NULL, "reduced", character(), FALSE, 4))
  expect_equal(loaded$required, saved$required)
  x$pooled$summary$parmatrices$Mean[1] <- -.7
  expect_error(suppressMessages(.ctmaPowerAnalyse(x$init, x$pooled, design, FALSE,
    "auto", NULL, "reduced", character(), FALSE, 4)), "does not match")
})

test_that("unsupported inputs fail clearly and missing powers stay missing", {
  x <- .power_fixture()
  expect_error(.power_design(x, c(3, 1, 1)), "timeRange")
  expect_error(.power_design(x, c(0, 0, 1)), "positive lag")
  expect_error(.power_design(x, probabilities = c(.8, NA)), "probabilities")
  x$init$statisticsList$allDeltas <- 1
  expect_error(.power_design(x), "mapping is ambiguous")
  x <- .power_fixture()
  design <- .power_design(x, c(.5, .5, .1), .8)
  result <- suppressMessages(.ctmaPowerAnalyse(x$init, x$pooled, design, FALSE,
    "auto", NULL, character(), character(), FALSE, 4))
  expect_true(all(is.na(result$posthoc[[1]][, 3])))
})

test_that("the public function loads cached fits and retains summary route notes", {
  x <- .power_fixture()
  directory <- withr::local_tempdir()
  saveRDS(x$pooled, file.path(directory, "pooled.rds"))
  expect_message(result <- ctmaPower(x$init, activeDirectory = directory,
    statisticalPower = c(.8, .95), timeRange = c(0, 2, 1), loadAllInvFit = "pooled"), "Direct covariance-matrix")
  expect_s3_class(result, "CoTiMAFit")
  expect_identical(result$power_method, "covariance")
  expect_match(result$summary$calculation_note, "unrestricted Gaussian")
  expect_equal(dim(result$summary$estimates[["Estimates of Model with all Effects Invariant"]]), c(3L, 12L))
  expect_equal(rownames(result$summary$estimates[["Estimates of Model with all Effects Invariant"]]),
    c("Fixed Effects Drift", "Fixed Effects Diffusion", "Fixed Effects T0Var"))
  expect_equal(result$summary$estimates[["Estimates of Model with all Effects Invariant"]][1, 1], -.3)
  expect_true("CoTiMAStanctArgs" %in% names(result$argumentList))
  expect_match(colnames(result$summary$estimates[["Power (post hoc) for Drift Effects"]][[1]])[3], "\u03b1", fixed = TRUE)
  expect_equal(result$regression_details$parameter, rep(.power_design(x)$targets$parameter, 2))
  plot_function <- get("ctmaPlot", envir = environment(ctmaPower))
  plot_environment <- new.env(parent = environment(plot_function))
  saved <- character()
  plot_environment$.ctmaPlotSave <- function(path) { saved <<- c(saved, path); invisible(path) }
  environment(plot_function) <- plot_environment
  grDevices::pdf(file = NULL); device <- grDevices::dev.cur()
  on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device))
  expect_no_error(plot_function(result, plotAutoEffects = FALSE, plotCrossEffects = FALSE))
  expect_length(saved, 2L)
  expect_true(all(vapply(.power_design(x)$targets$parameter,
    function(name) any(grepl(name, saved, fixed = TRUE)), logical(1))))
})

test_that("the public orchestration delegates fitting once with its controls", {
  x <- .power_fixture()
  environment <- new.env(parent = environment(ctmaPower))
  calls <- list()
  environment$ctmaAllInvFit <- function(...) { calls[[length(calls) + 1L]] <<- list(...); x$pooled }
  f <- ctmaPower; environment(f) <- environment
  suppressMessages(result <- f(x$init, activeDirectory = tempdir(), statisticalPower = c(.8, .95),
    timeRange = c(1, 1, 1), coresToUse = 1, finishsamples = 123, optimize = TRUE, priors = FALSE))
  expect_length(calls, 1L)
  expect_identical(calls[[1]]$finishsamples, 123)
  expect_false(calls[[1]]$priors)
  expect_identical(result$power_method, "covariance")
})
