.ctmaModelText <- function(model, object_name = "m") {
  if (!is.list(model) || !is.data.frame(model$pars) ||
      !is.list(model$matrices)) {
    stop("'model' must be a ctsem model with 'pars' and 'matrices'.", call. = FALSE)
  }
  if (!is.character(object_name) || length(object_name) != 1L ||
      is.na(object_name) || !nzchar(object_name) ||
      !identical(make.names(object_name), object_name)) {
    stop("'object_name' must be one valid R object name.", call. = FALSE)
  }

  matrices <- model$matrices
  required <- c(
    "LAMBDA", "T0VAR", "T0MEANS", "MANIFESTMEANS", "MANIFESTVAR",
    "DRIFT", "CINT", "DIFFUSION"
  )
  missing <- setdiff(required, names(matrices))
  if (length(missing)) {
    stop(
      "The ctsem model is missing matrix specifications: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }

  args <- list(
    type = if (isTRUE(model$continuoustime)) "ct" else "dt",
    n.manifest = as.integer(model$n.manifest),
    n.latent = as.integer(model$n.latent),
    Tpoints = model$Tpoints,
    manifestNames = model$manifestNames,
    manifesttype = model$manifesttype,
    latentNames = model$latentNames,
    id = model$subjectIDname,
    time = model$timeName,
    LAMBDA = matrices$LAMBDA,
    T0VAR = matrices$T0VAR,
    T0MEANS = matrices$T0MEANS,
    MANIFESTMEANS = matrices$MANIFESTMEANS,
    MANIFESTVAR = matrices$MANIFESTVAR,
    DRIFT = matrices$DRIFT,
    CINT = matrices$CINT,
    DIFFUSION = matrices$DIFFUSION,
    n.TDpred = as.integer(model$n.TDpred),
    TDpredNames = model$TDpredNames,
    n.TIpred = as.integer(model$n.TIpred),
    TIpredNames = model$TIpredNames
  )

  optional <- c(
    "Tpoints", "manifestNames", "manifesttype", "latentNames", "id", "time",
    "TDpredNames", "TIpredNames"
  )
  for (name in optional) {
    if (is.null(args[[name]]) || !length(args[[name]])) args[[name]] <- NULL
  }
  if (is.null(args$TDpredNames)) args$n.TDpred <- 0L
  if (is.null(args$TIpredNames)) args$n.TIpred <- 0L

  if (length(model$ncategories)) args$ncategories <- model$ncategories
  if (length(model$censormin)) args$censormin <- model$censormin
  if (length(model$censormax)) args$censormax <- model$censormax

  # Emit ordinary constructors so users can edit the model specification.
  as_code <- function(value) {
    if (is.matrix(value)) {
      parts <- list(
        quote(matrix), as_code(as.vector(value)),
        nrow = as.numeric(nrow(value)), ncol = as.numeric(ncol(value))
      )
      if (!is.null(dimnames(value))) {
        parts$dimnames <- as_code(dimnames(value))
      }
      return(as.call(parts))
    }
    if (is.data.frame(value)) {
      columns <- lapply(value, as_code)
      options <- list(check.names = FALSE, stringsAsFactors = FALSE)
      if (!identical(rownames(value), as.character(seq_len(nrow(value))))) {
        options$row.names <- as_code(rownames(value))
      }
      return(as.call(c(list(quote(data.frame)), columns, options)))
    }
    if (is.list(value)) {
      return(as.call(c(list(quote(list)), lapply(value, as_code))))
    }
    if (is.factor(value)) {
      return(call("factor", as_code(as.character(value)),
        levels = as_code(levels(value)), ordered = is.ordered(value)))
    }
    if (is.integer(value)) value <- as.numeric(value)
    parse(text = paste(capture.output(dput(value)), collapse = "\n"))[[1L]]
  }
  args <- lapply(args, as_code)
  model_call <- as.call(c(list(quote(ctsem::ctModel)), args))
  assignment <- call("<-", as.name(object_name), model_call)
  model_text <- paste(
    deparse(
      assignment,
      width.cutoff = 120L,
      control = c("keepNA", "keepInteger", "niceNames")
    ),
    collapse = "\n"
  )
  pars_assignment <- call(
    "<-",
    call("$", as.name(object_name), as.name("pars")),
    as_code(model$pars)
  )
  pars_text <- paste(
    deparse(
      pars_assignment,
      width.cutoff = 120L,
      control = c("keepNA", "keepInteger", "niceNames")
    ),
    collapse = "\n"
  )
  paste(
    model_text,
    "",
    "# Preserve CoTiMA's parameter table and its post-construction constraints.",
    pars_text,
    sep = "\n"
  )
}
