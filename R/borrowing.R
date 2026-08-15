.endpoint_data <- function(data, endpoint) {
  if (endpoint == "toxicity") {
    data.frame(
      study = data$study, dose = data$dose, x = data$x,
      y = data$toxicity, n = data$n_toxicity,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      study = data$study, dose = data$dose, x = data$x,
      y = data$efficacy, n = data$n_efficacy,
      stringsAsFactors = FALSE
    )
  }
}

.conflict_weights <- function(target_fit, donor_fits, donor_data, endpoint, control) {
  if (nrow(donor_data) == 0L) {
    return(data.frame(
      endpoint = character(), donor = character(), dose = numeric(), x = numeric(),
      events = numeric(), total = numeric(), conflict = numeric(),
      raw_weight = numeric(), donor_multiplier = numeric(),
      total_multiplier = numeric(), weight = numeric(), borrowed_n = numeric()
    ))
  }
  parts <- lapply(names(donor_fits), function(donor) {
    donor_rows <- donor_data[donor_data$study == donor, , drop = FALSE]
    target_moments <- .linear_predictor_moments(target_fit, donor_rows$x)
    donor_moments <- .linear_predictor_moments(donor_fits[[donor]], donor_rows$x)
    conflict <- abs(target_moments$mean - donor_moments$mean) /
      sqrt(target_moments$variance + donor_moments$variance +
        control$conflict_scale[[endpoint]]^2)
    data.frame(
      endpoint = endpoint,
      donor = donor,
      dose = donor_rows$dose,
      x = donor_rows$x,
      events = donor_rows$y,
      total = donor_rows$n,
      conflict = conflict,
      raw_weight = exp(-0.5 * conflict^2),
      stringsAsFactors = FALSE
    )
  })
  weights <- do.call(rbind, parts)
  target_n <- sum(target_fit$data$n)
  donor_n <- tapply(weights$total, weights$donor, sum)
  donor_multiplier <- pmin(1, control$donor_cap * target_n / donor_n)
  names(donor_multiplier) <- names(donor_n)
  weights$donor_multiplier <- unname(donor_multiplier[weights$donor])
  borrowed_after_donor_cap <- sum(
    weights$total * weights$raw_weight * weights$donor_multiplier
  )
  total_multiplier <- if (borrowed_after_donor_cap > 0) {
    min(1, control$total_cap * target_n / borrowed_after_donor_cap)
  } else {
    1
  }
  weights$total_multiplier <- total_multiplier
  weights$weight <- weights$raw_weight * weights$donor_multiplier * total_multiplier
  weights$borrowed_n <- weights$weight * weights$total
  rownames(weights) <- NULL
  weights
}

.fit_borrowed_endpoint <- function(endpoint_data, target, donors, endpoint, prior, control) {
  studies <- c(target, donors)
  local_fits <- stats::setNames(lapply(studies, function(study) {
    rows <- endpoint_data[endpoint_data$study == study, , drop = FALSE]
    .fit_endpoint(rows, endpoint, prior, control)
  }), studies)
  target_rows <- endpoint_data[endpoint_data$study == target, , drop = FALSE]
  donor_rows <- endpoint_data[endpoint_data$study %in% donors, , drop = FALSE]
  weights <- .conflict_weights(
    local_fits[[target]], local_fits[donors], donor_rows, endpoint, control
  )
  if (nrow(weights) == 0L) {
    borrowed_fit <- local_fits[[target]]
  } else {
    borrowed_data <- data.frame(
      study = weights$donor, dose = weights$dose, x = weights$x,
      y = weights$events, n = weights$total,
      stringsAsFactors = FALSE
    )
    combined <- rbind(target_rows, borrowed_data)
    cell_weight <- c(rep(1, nrow(target_rows)), weights$weight)
    borrowed_fit <- .fit_endpoint(combined, endpoint, prior, control, cell_weight)
  }
  list(
    local_target = local_fits[[target]],
    local_donors = local_fits[donors],
    borrowed = borrowed_fit,
    weights = weights
  )
}

#' Summarize selective borrowing
#'
#' @param object A fitted `safeab_fit` object.
#' @param by Summarize by endpoint and donor, or return dose-level weights.
#' @return A data frame describing conflicts, weights, and borrowed effective
#'   sample sizes.
#' @export
borrowing_summary <- function(object, by = c("donor", "dose")) {
  if (!inherits(object, "safeab_fit")) stop("object must be a safeab_fit.", call. = FALSE)
  by <- match.arg(by)
  weights <- do.call(rbind, object$weights)
  if (is.null(weights) || nrow(weights) == 0L) return(weights)
  rownames(weights) <- NULL
  if (by == "dose") return(weights)
  keys <- interaction(weights$endpoint, weights$donor, drop = TRUE)
  pieces <- lapply(split(weights, keys), function(rows) {
    data.frame(
      endpoint = rows$endpoint[1L],
      donor = rows$donor[1L],
      available_n = sum(rows$total),
      borrowed_n = sum(rows$borrowed_n),
      relative_weight = sum(rows$borrowed_n) / sum(rows$total),
      mean_conflict = stats::weighted.mean(rows$conflict, rows$total),
      max_conflict = max(rows$conflict),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, pieces)
  rownames(result) <- NULL
  result[order(result$endpoint, result$donor), , drop = FALSE]
}

#' Borrowed effective sample size
#'
#' @param object A fitted `safeab_fit` object.
#' @param endpoint Optional endpoint name.
#' @return A named numeric vector of borrowed effective sample sizes.
#' @export
effective_sample_size <- function(object, endpoint = NULL) {
  if (!inherits(object, "safeab_fit")) stop("object must be a safeab_fit.", call. = FALSE)
  if (!is.null(endpoint)) endpoint <- match.arg(endpoint, object$endpoints)
  selected <- if (is.null(endpoint)) object$endpoints else endpoint
  stats::setNames(vapply(selected, function(name) {
    weights <- object$weights[[name]]
    if (is.null(weights) || nrow(weights) == 0L) 0 else sum(weights$borrowed_n)
  }, numeric(1)), selected)
}

#' Inspect target-donor commensurability
#'
#' @param object A fitted `safeab_fit` object.
#' @return Dose-level standardized conflicts and raw compatibility weights.
#' @export
commensurability <- function(object) {
  details <- borrowing_summary(object, by = "dose")
  if (is.null(details) || nrow(details) == 0L) return(details)
  details[, c("endpoint", "donor", "dose", "x", "conflict", "raw_weight")]
}
