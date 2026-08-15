#' Specify priors for SAFE-AB endpoint models
#'
#' @param toxicity_mean,toxicity_sd Numeric vectors of length two for the
#'   monotone toxicity model.
#' @param efficacy_mean,efficacy_sd Numeric vectors of length three for the
#'   quadratic efficacy model.
#' @return A `safeab_prior` object.
#' @export
safeab_prior <- function(
    toxicity_mean = c(stats::qlogis(0.10), 0),
    toxicity_sd = c(1.5, 0.75),
    efficacy_mean = c(stats::qlogis(0.15), 0, 0),
    efficacy_sd = c(1.5, 2, 2)) {
  check_prior <- function(mean, sd, n, endpoint) {
    if (!is.numeric(mean) || length(mean) != n || any(!is.finite(mean))) {
      stop(endpoint, " prior mean must contain ", n, " finite values.", call. = FALSE)
    }
    if (!is.numeric(sd) || length(sd) != n || any(!is.finite(sd)) || any(sd <= 0)) {
      stop(endpoint, " prior sd must contain ", n, " positive values.", call. = FALSE)
    }
  }
  check_prior(toxicity_mean, toxicity_sd, 2L, "Toxicity")
  check_prior(efficacy_mean, efficacy_sd, 3L, "Efficacy")
  structure(list(
    toxicity = list(mean = toxicity_mean, sd = toxicity_sd),
    efficacy = list(mean = efficacy_mean, sd = efficacy_sd)
  ), class = "safeab_prior")
}

#' Control SAFE-AB fitting and borrowing
#'
#' @param conflict_scale Positive endpoint-specific conflict scales.
#' @param donor_cap Maximum borrowed effective sample size from each donor,
#'   expressed as a fraction of the target endpoint sample size.
#' @param total_cap Maximum total borrowed effective sample size, expressed as
#'   a fraction of the target endpoint sample size.
#' @param quadrature_points Number of adaptive Gauss-Hermite quadrature points.
#' @return A `safeab_control` object.
#' @export
safeab_control <- function(
    conflict_scale = c(toxicity = 0.50, efficacy = 0.50),
    donor_cap = 0.25,
    total_cap = 1,
    quadrature_points = 7L) {
  if (is.null(names(conflict_scale)) && length(conflict_scale) == 2L) {
    names(conflict_scale) <- c("toxicity", "efficacy")
  }
  required <- c("toxicity", "efficacy")
  if (!all(required %in% names(conflict_scale)) ||
      any(!is.finite(conflict_scale[required])) || any(conflict_scale[required] <= 0)) {
    stop("conflict_scale must contain positive toxicity and efficacy values.", call. = FALSE)
  }
  caps <- c(donor_cap, total_cap)
  if (any(!is.finite(caps)) || any(caps < 0)) {
    stop("BESS caps must be non-negative finite values.", call. = FALSE)
  }
  quadrature_points <- as.integer(quadrature_points)
  if (length(quadrature_points) != 1L || is.na(quadrature_points) || quadrature_points < 3L) {
    stop("quadrature_points must be an integer of at least 3.", call. = FALSE)
  }
  structure(list(
    conflict_scale = conflict_scale[required],
    donor_cap = donor_cap,
    total_cap = total_cap,
    quadrature_points = quadrature_points
  ), class = "safeab_control")
}

.safeab_starts <- function(endpoint) {
  if (endpoint == "toxicity") {
    rbind(
      c(stats::qlogis(0.10), 0), c(stats::qlogis(0.05), -0.5),
      c(stats::qlogis(0.20), 0.5), c(stats::qlogis(0.05), 0.5),
      c(stats::qlogis(0.20), -0.5)
    )
  } else {
    rbind(
      c(stats::qlogis(0.15), 0, 0), c(stats::qlogis(0.05), -1, 0),
      c(stats::qlogis(0.30), 1, 0), c(stats::qlogis(0.15), 0, -1),
      c(stats::qlogis(0.15), 0, 1)
    )
  }
}
