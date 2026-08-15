#' Construct the safety-admissible dose set
#'
#' @param object A fitted `safeab_fit` containing toxicity.
#' @param doses Candidate original or standardized dose values. Target-study
#'   observed doses are used by default.
#' @param toxicity_limit Maximum acceptable toxicity probability.
#' @param borrowed_cutoff Maximum posterior exceedance probability under the
#'   selectively borrowed model.
#' @param local_cutoff Maximum posterior exceedance probability under the
#'   target-only local model.
#' @param local_veto Apply the target-only safety veto.
#' @param dose_scale Scale used by `doses`.
#' @return A dose-level safety table.
#' @export
admissible_doses <- function(
    object,
    doses = NULL,
    toxicity_limit = 0.30,
    borrowed_cutoff = 0.10,
    local_cutoff = 0.20,
    local_veto = TRUE,
    dose_scale = c("dose", "standardized")) {
  if (!inherits(object, "safeab_fit") || !"toxicity" %in% object$endpoints) {
    stop("object must be a safeab_fit containing toxicity.", call. = FALSE)
  }
  cutoffs <- c(toxicity_limit, borrowed_cutoff, local_cutoff)
  if (any(!is.finite(cutoffs)) || any(cutoffs <= 0 | cutoffs >= 1)) {
    stop("Toxicity limit and posterior cutoffs must lie in (0, 1).", call. = FALSE)
  }
  dose_scale <- match.arg(dose_scale)
  grid <- .as_prediction_grid(object, doses, dose_scale)
  borrowed <- .posterior_summary(
    object$endpoint_fits$toxicity, grid$x, threshold = toxicity_limit
  )
  local <- .posterior_summary(
    object$local_fits$toxicity, grid$x, threshold = toxicity_limit
  )
  admissible <- borrowed$probability_exceeding <= borrowed_cutoff
  if (isTRUE(local_veto)) {
    admissible <- admissible & local$probability_exceeding <= local_cutoff
  }
  data.frame(
    dose = grid$dose,
    x = grid$x,
    toxicity = borrowed$estimate,
    toxicity_lower = borrowed$lower,
    toxicity_upper = borrowed$upper,
    borrowed_exceedance = borrowed$probability_exceeding,
    local_exceedance = local$probability_exceeding,
    local_veto = isTRUE(local_veto) & local$probability_exceeding > local_cutoff,
    admissible = admissible
  )
}

#' Recommend a dose under SAFE-AB
#'
#' @param object A fitted `safeab_fit` containing toxicity and efficacy.
#' @param doses Candidate doses. Target-study observed doses are used by default.
#' @param toxicity_limit,borrowed_cutoff,local_cutoff,local_veto Safety-rule
#'   arguments passed to [admissible_doses()].
#' @param efficacy_margin Select the lowest admissible dose whose posterior mean
#'   efficacy is within this absolute margin of the best admissible dose.
#' @param utility Optional function receiving the completed dose table and
#'   returning one finite utility value per row. The admissible dose with maximum
#'   utility is selected.
#' @param dose_scale Scale used by `doses`.
#' @return A `safeab_recommendation` object with the recommendation and dose table.
#' @export
recommend_dose <- function(
    object,
    doses = NULL,
    toxicity_limit = 0.30,
    borrowed_cutoff = 0.10,
    local_cutoff = 0.20,
    local_veto = TRUE,
    efficacy_margin = 0.05,
    utility = NULL,
    dose_scale = c("dose", "standardized")) {
  if (!inherits(object, "safeab_fit") ||
      !all(c("toxicity", "efficacy") %in% object$endpoints)) {
    stop("Dose recommendation requires fitted toxicity and efficacy endpoints.", call. = FALSE)
  }
  if (!is.numeric(efficacy_margin) || length(efficacy_margin) != 1L ||
      !is.finite(efficacy_margin) || efficacy_margin < 0) {
    stop("efficacy_margin must be one non-negative finite value.", call. = FALSE)
  }
  dose_scale <- match.arg(dose_scale)
  safety <- admissible_doses(
    object, doses, toxicity_limit, borrowed_cutoff, local_cutoff,
    local_veto, dose_scale
  )
  efficacy <- .posterior_summary(object$endpoint_fits$efficacy, safety$x)
  safety$efficacy <- efficacy$estimate
  safety$efficacy_lower <- efficacy$lower
  safety$efficacy_upper <- efficacy$upper
  candidate <- which(safety$admissible)
  selected <- integer()
  if (length(candidate)) {
    if (is.null(utility)) {
      best <- max(safety$efficacy[candidate])
      eligible <- candidate[safety$efficacy[candidate] >= best - efficacy_margin]
      selected <- eligible[which.min(safety$x[eligible])]
      safety$utility <- NA_real_
    } else {
      if (!is.function(utility)) stop("utility must be a function or NULL.", call. = FALSE)
      value <- utility(safety)
      if (!is.numeric(value) || length(value) != nrow(safety) || any(!is.finite(value))) {
        stop("utility must return one finite numeric value per candidate dose.", call. = FALSE)
      }
      safety$utility <- value
      selected <- candidate[which.max(value[candidate])]
    }
  } else {
    safety$utility <- NA_real_
  }
  result <- list(
    recommended_dose = if (length(selected)) safety$dose[selected] else NA_real_,
    recommended_x = if (length(selected)) safety$x[selected] else NA_real_,
    selected_row = if (length(selected)) selected else NA_integer_,
    dose_table = safety,
    settings = list(
      toxicity_limit = toxicity_limit,
      borrowed_cutoff = borrowed_cutoff,
      local_cutoff = local_cutoff,
      local_veto = local_veto,
      efficacy_margin = efficacy_margin
    )
  )
  structure(result, class = "safeab_recommendation")
}

#' @export
print.safeab_recommendation <- function(x, ...) {
  if (is.na(x$recommended_x)) {
    cat("SAFE-AB recommendation: no admissible dose\n")
  } else if (is.na(x$recommended_dose)) {
    cat("SAFE-AB recommendation (standardized dose):", format(x$recommended_x), "\n")
  } else {
    cat("SAFE-AB recommended dose:", format(x$recommended_dose), "\n")
  }
  invisible(x)
}
