# Validate the DRIFT standard errors used in bias and heterogeneity analyses.
# DIFFUSION and T0VAR may contain fixed cells with legitimate zero SEs.
.ctmaBiGCheckDriftSE <- function(all_SE, n.latent, excludeStudies = integer()) {
  columns <- seq_len(n.latent^2)
  studies <- setdiff(seq_len(nrow(all_SE)), excludeStudies)
  drift_se <- all_SE[studies, columns, drop = FALSE]
  invalid <- which(!is.finite(drift_se) | drift_se <= 0, arr.ind = TRUE)
  if (!nrow(invalid)) return(invisible(TRUE))
  parameters <- colnames(all_SE)[columns]
  if (is.null(parameters)) {
    parameters <- paste0("DRIFT[", rep(seq_len(n.latent), each = n.latent),
      ",", rep(seq_len(n.latent), n.latent), "]")
  }
  details <- data.frame(
    Study = studies[invalid[, 1]],
    Matrix = "DRIFT",
    Parameter = parameters[invalid[, 2]],
    SE = drift_se[invalid],
    check.names = FALSE
  )
  stop(
    "ctmaBiG requires finite, strictly positive DRIFT standard errors.\n",
    paste(capture.output(print(details, row.names = FALSE)), collapse = "\n"),
    "\nInspect these parameters; exclude a study only when substantively justified.",
    call. = FALSE
  )
}
