# Restore the requested observation families after a model-construction branch
# has replaced the original ctsem model (for example, the augmented-state
# random-intercept representation).
.ctma_restore_manifesttype <- function(model, binaries) {
  if (is.null(binaries)) return(model)
  nmanifest <- as.integer(model$n.manifest)[1L]
  if (!is.finite(nmanifest) || nmanifest < 1L || length(binaries) < nmanifest)
    stop("binaries must contain at least one entry per manifest variable.",
      call. = FALSE)
  model$manifesttype <- as.numeric(binaries[seq_len(nmanifest)])
  model
}

# Validate location parameters for binary manifest variables.
#
# This check must use the final ctsem model. In particular, ctmaFit() can
# replace its initial model with a process-based random-intercept model.
.ctma_check_binary_intercepts <- function(model, binaries,
    process_intercepts = FALSE, inform = TRUE) {
  if (is.null(binaries)) return(invisible(model))
  if (is.null(model$pars) || !is.data.frame(model$pars)) stop(
    "The ctsem model has no parameter table for binary-intercept validation.",
    call. = FALSE
  )

  nmanifest <- as.integer(model$n.manifest)[1L]
  nlatent <- as.integer(model$n.latent)[1L]
  if (!is.finite(nmanifest) || nmanifest < 1L ||
      !is.finite(nlatent) || nlatent < 1L) stop(
    "The ctsem model does not report valid manifest and latent dimensions.",
    call. = FALSE
  )
  if (length(binaries) < nmanifest) stop(
    "binaries must contain at least one entry per manifest variable.",
    call. = FALSE
  )
  binary_manifest <- which(as.integer(binaries[seq_len(nmanifest)]) == 1L)
  if (!length(binary_manifest) || isTRUE(process_intercepts))
    return(invisible(model))

  pars <- model$pars
  free_by_row <- function(matrix_name, nr) {
    answer <- rep(FALSE, nr)
    selected <- pars[!is.na(pars$matrix) & pars$matrix == matrix_name,
      , drop = FALSE]
    if (!nrow(selected)) return(answer)
    for (i in seq_len(nrow(selected))) {
      r <- as.integer(selected$row[i])
      if (!is.na(r) && r >= 1L && r <= nr)
        answer[r] <- answer[r] || !is.na(selected$param[i])
    }
    answer
  }

  manifest_mean_free <- free_by_row("MANIFESTMEANS", nmanifest)
  cint_free <- free_by_row("CINT", nlatent)

  # Map latent CINTs to manifest variables through active loadings instead of
  # assuming that manifest and latent variable indices are interchangeable.
  lambda_active <- matrix(FALSE, nmanifest, nlatent)
  lambda <- pars[!is.na(pars$matrix) & pars$matrix == "LAMBDA",
    , drop = FALSE]
  if (nrow(lambda)) for (i in seq_len(nrow(lambda))) {
    r <- as.integer(lambda$row[i]); cc <- as.integer(lambda$col[i])
    if (is.na(r) || is.na(cc) || r < 1L || r > nmanifest ||
        cc < 1L || cc > nlatent) next
    is_free <- !is.na(lambda$param[i])
    is_fixed_nonzero <- !is.na(lambda$value[i]) && lambda$value[i] != 0
    lambda_active[r, cc] <- is_free || is_fixed_nonzero
  }
  cint_support <- vapply(seq_len(nmanifest), function(j)
    any(lambda_active[j, ] & cint_free), logical(1L))

  has_cint <- cint_support[binary_manifest]
  has_manifest_mean <- manifest_mean_free[binary_manifest]
  labels <- model$manifestNames
  if (is.null(labels) || length(labels) < nmanifest)
    labels <- paste0("manifest ", seq_len(nmanifest))
  describe <- function(index) paste(labels[index], collapse = ", ")

  missing <- !has_cint & !has_manifest_mean
  if (any(missing)) stop(
    "Binary variable(s) ", describe(binary_manifest[missing]),
    " have neither a freely estimated MANIFESTMEANS parameter nor a freely ",
    "estimated CINT connected through LAMBDA. Prefer freely estimated CINTs ",
    "with MANIFESTMEANS fixed to zero, or explicitly free MANIFESTMEANS.",
    call. = FALSE
  )

  manifest_only <- !has_cint & has_manifest_mean
  if (isTRUE(inform) && any(manifest_only)) message(
    "Binary variable(s) ", describe(binary_manifest[manifest_only]),
    " use freely estimated MANIFESTMEANS while their connected latent CINTs ",
    "are fixed."
  )

  both <- has_cint & has_manifest_mean
  if (isTRUE(inform) && any(both)) warning(
    "Binary variable(s) ", describe(binary_manifest[both]),
    " have both freely estimated MANIFESTMEANS and connected CINT parameters. ",
    "Check location identification.",
    call. = FALSE
  )
  invisible(model)
}
