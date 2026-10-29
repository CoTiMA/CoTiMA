test_that("zero-sum observed rows are not classified as all missing", {
  data <- data.frame(
    first = c(-1, NA, 0),
    second = c(1, NA, 0)
  )

  expect_identical(
    CoTiMA:::.ctmaRowsAllMissing(data, c("first", "second")),
    c(FALSE, TRUE, FALSE)
  )
})

test_that("within-time-A detects missing T0 values rather than zero sums", {
  data <- data.frame(
    first_T0 = c(-1, NA, 0),
    second_T0 = c(1, 1, 0)
  )

  expect_identical(
    CoTiMA:::.ctmaRowsAnyMissing(data, c("first_T0", "second_T0")),
    c(FALSE, TRUE, FALSE)
  )
})

test_that("time-point matching distinguishes T1 from T10 and T11", {
  expect_identical(
    CoTiMA:::.ctmaColumnsAtTime(c("X_T1", "X_T10", "X_T11"), 1),
    "X_T1"
  )

  time_points <- 12L
  process_names <- unlist(lapply(0:11, function(i) paste0(c("x", "y"), i)))
  process_data <- as.data.frame(
    matrix(seq_len(2 * time_points * 2), nrow = 2, byrow = TRUE)
  )
  names(process_data) <- process_names
  time_data <- as.data.frame(
    rbind(c(0, 0.0001, 2:11), c(0, 0.0001, 2:11))
  )
  names(time_data) <- paste0("time_input_", 0:11)

  shaped <- suppressMessages(CoTiMA::ctmaShapeRawData(
    dataFrame = cbind(process_data, time_data),
    inputDataFrameFormat = "wide",
    inputTimeFormat = "time",
    n.manifest = 2,
    Tpoints = time_points,
    orderInputVariablesNames = "time",
    targetInputVariablesNames = process_names,
    targetTimeVariablesNames = names(time_data),
    outputDataFrameFormat = "wide",
    outputVariablesNames = c("X", "Y"),
    outputTimeFormat = "time",
    minTolDelta = 0.0002
  ))

  expect_true(all(is.na(shaped[, c("X_T1", "Y_T1")])))
  expect_false(anyNA(shaped[, c("X_T10", "Y_T10", "X_T11", "Y_T11")]))
})
