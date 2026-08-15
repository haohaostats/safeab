.column_value <- function(data, column, label) {
  if (!is.character(column) || length(column) != 1L || !column %in% names(data)) {
    stop(label, " must name a column in data.", call. = FALSE)
  }
  data[[column]]
}

.validate_binomial <- function(events, total, endpoint) {
  if (!is.numeric(events) || !is.numeric(total) ||
      any(!is.finite(events)) || any(!is.finite(total))) {
    stop(endpoint, " event and total columns must be finite numeric values.", call. = FALSE)
  }
  if (any(events < 0) || any(total <= 0) || any(events > total) ||
      any(events != floor(events)) || any(total != floor(total))) {
    stop(endpoint, " requires integer values satisfying 0 <= events <= total.", call. = FALSE)
  }
}

#' Prepare aggregated dose-response data for SAFE-AB
#'
#' @param data A data frame with one row per study-dose combination.
#' @param target Value of the study column identifying the target study.
#' @param study,dose Column names for the study identifier and numerical dose.
#' @param toxicity_events,toxicity_total Column names for toxicity event counts
#'   and evaluated totals. Both may be `NULL` when toxicity is unavailable.
#' @param efficacy_events,efficacy_total Column names for efficacy event counts
#'   and evaluated totals. Both may be `NULL` when efficacy is unavailable.
#' @param dose_transform Dose transformation. The default, `"rank"`, maps the
#'   ordered doses within each study to equally spaced positions on `[0, 1]`.
#'   Alternatives are pooled range standardization, pooled log-dose range
#'   standardization, or no transformation. With `"none"`, dose values must
#'   already lie in `[0, 1]`.
#' @return A validated `safeab_data` object.
#' @export
safeab_data <- function(
    data,
    target,
    study = "study",
    dose = "dose",
    toxicity_events = "toxicity",
    toxicity_total = "n_toxicity",
    efficacy_events = "efficacy",
    efficacy_total = "n_efficacy",
    dose_transform = c("rank", "range", "log_range", "none")) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("data must be a non-empty data frame.", call. = FALSE)
  }
  dose_transform <- match.arg(dose_transform)
  study_value <- as.character(.column_value(data, study, "study"))
  dose_value <- as.numeric(.column_value(data, dose, "dose"))
  if (anyNA(study_value) || any(!nzchar(study_value))) {
    stop("Study identifiers cannot be missing or empty.", call. = FALSE)
  }
  if (any(!is.finite(dose_value))) stop("Dose values must be finite.", call. = FALSE)
  target <- as.character(target)
  if (length(target) != 1L || !target %in% study_value) {
    stop("target must identify exactly one study present in data.", call. = FALSE)
  }
  key <- paste(study_value, dose_value, sep = "\r")
  if (anyDuplicated(key)) {
    stop("Each study-dose combination must occupy exactly one row.", call. = FALSE)
  }

  transform_input <- dose_value
  if (dose_transform == "rank") {
    x <- numeric(length(dose_value))
    for (study_id in unique(study_value)) {
      index <- which(study_value == study_id)
      if (length(index) < 2L) {
        stop(
          "Within-study rank standardization requires at least two doses in every study; ",
          study_id, " has only one.", call. = FALSE
        )
      }
      ordered <- index[order(dose_value[index])]
      x[ordered] <- (seq_along(ordered) - 1) / (length(ordered) - 1)
    }
    transform_info <- list(type = "rank")
  } else if (dose_transform == "log_range") {
    if (any(dose_value <= 0)) stop("log_range requires strictly positive doses.", call. = FALSE)
    transform_input <- log(dose_value)
  }
  if (dose_transform == "none") {
    x <- transform_input
    if (any(x < 0 | x > 1)) stop("With dose_transform='none', doses must lie in [0, 1].", call. = FALSE)
    transform_info <- list(type = "none", minimum = 0, maximum = 1)
  } else if (dose_transform != "rank") {
    limits <- range(transform_input)
    if (diff(limits) <= 0) stop("At least two distinct pooled doses are required.", call. = FALSE)
    x <- (transform_input - limits[1L]) / diff(limits)
    transform_info <- list(type = dose_transform, minimum = limits[1L], maximum = limits[2L])
  }

  out <- data.frame(study = study_value, dose = dose_value, x = x, stringsAsFactors = FALSE)
  has_tox <- !is.null(toxicity_events) || !is.null(toxicity_total)
  has_eff <- !is.null(efficacy_events) || !is.null(efficacy_total)
  if (xor(is.null(toxicity_events), is.null(toxicity_total))) {
    stop("Provide both toxicity columns or set both to NULL.", call. = FALSE)
  }
  if (xor(is.null(efficacy_events), is.null(efficacy_total))) {
    stop("Provide both efficacy columns or set both to NULL.", call. = FALSE)
  }
  if (!has_tox && !has_eff) stop("At least one endpoint must be supplied.", call. = FALSE)
  if (has_tox) {
    out$toxicity <- as.numeric(.column_value(data, toxicity_events, "toxicity_events"))
    out$n_toxicity <- as.numeric(.column_value(data, toxicity_total, "toxicity_total"))
    .validate_binomial(out$toxicity, out$n_toxicity, "Toxicity")
  }
  if (has_eff) {
    out$efficacy <- as.numeric(.column_value(data, efficacy_events, "efficacy_events"))
    out$n_efficacy <- as.numeric(.column_value(data, efficacy_total, "efficacy_total"))
    .validate_binomial(out$efficacy, out$n_efficacy, "Efficacy")
  }
  out <- out[order(out$study, out$x), , drop = FALSE]
  rownames(out) <- NULL
  structure(list(
    data = out,
    target = target,
    donors = setdiff(unique(out$study), target),
    endpoints = c(if (has_tox) "toxicity", if (has_eff) "efficacy"),
    dose_transform = transform_info
  ), class = "safeab_data")
}

#' @export
print.safeab_data <- function(x, ...) {
  cat("SAFE-AB data\n")
  cat("  Target:", x$target, "\n")
  cat("  Donors:", if (length(x$donors)) paste(x$donors, collapse = ", ") else "none", "\n")
  cat("  Endpoints:", paste(x$endpoints, collapse = ", "), "\n")
  cat("  Study-dose rows:", nrow(x$data), "\n")
  invisible(x)
}
