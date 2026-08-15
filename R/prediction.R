.as_prediction_grid <- function(object, newdata, dose_scale) {
  dose_scale <- match.arg(dose_scale, c("dose", "standardized"))
  if (missing(newdata) || is.null(newdata)) {
    target_rows <- object$data$data[object$data$data$study == object$target, , drop = FALSE]
    return(data.frame(dose = target_rows$dose, x = target_rows$x))
  }
  values <- as.numeric(newdata)
  if (any(!is.finite(values))) stop("newdata must contain finite doses.", call. = FALSE)
  if (dose_scale == "standardized") {
    x <- values
    dose <- rep(NA_real_, length(values))
  } else {
    dose <- values
    transform <- object$data$dose_transform
    if (transform$type == "rank") {
      target_rows <- object$data$data[
        object$data$data$study == object$target, c("dose", "x"), drop = FALSE
      ]
      matched <- match(values, target_rows$dose)
      if (anyNA(matched)) {
        stop(
          "With within-study rank standardization, original-dose predictions ",
          "must use doses observed in the target study; use ",
          "dose_scale='standardized' for an arbitrary grid.", call. = FALSE
        )
      }
      x <- target_rows$x[matched]
    } else {
      transformed <- if (transform$type == "log_range") {
        if (any(values <= 0)) stop("Log-dose prediction requires positive doses.", call. = FALSE)
        log(values)
      } else values
      x <- if (transform$type == "none") transformed else
        (transformed - transform$minimum) / (transform$maximum - transform$minimum)
    }
  }
  if (any(x < 0 | x > 1)) warning("Some prediction doses lie outside the fitted standardized range.", call. = FALSE)
  data.frame(dose = dose, x = x)
}

#' Posterior dose-response predictions
#'
#' @param object A fitted `safeab_fit` object.
#' @param newdata Numeric dose values. The target-study observed doses are used
#'   by default.
#' @param endpoint Endpoint to predict. If omitted, all fitted endpoints are used.
#' @param model Return borrowed, target-only local, or both posterior predictions.
#' @param dose_scale Whether `newdata` contains original or standardized doses.
#' @param probs Two posterior interval probabilities.
#' @param ... Unused.
#' @return A data frame of posterior probabilities and intervals.
#' @export
predict.safeab_fit <- function(
    object,
    newdata = NULL,
    endpoint = NULL,
    model = c("borrowed", "local", "both"),
    dose_scale = c("dose", "standardized"),
    probs = c(0.025, 0.975),
    ...) {
  model <- match.arg(model)
  dose_scale <- match.arg(dose_scale)
  if (is.null(endpoint)) endpoint <- object$endpoints
  endpoint <- unique(as.character(endpoint))
  if (length(setdiff(endpoint, object$endpoints))) {
    stop("Requested endpoint was not fitted.", call. = FALSE)
  }
  if (!is.numeric(probs) || length(probs) != 2L || any(probs <= 0 | probs >= 1) || probs[1L] >= probs[2L]) {
    stop("probs must contain two increasing probabilities in (0, 1).", call. = FALSE)
  }
  grid <- .as_prediction_grid(object, newdata, dose_scale)
  models <- if (model == "both") c("local", "borrowed") else model
  pieces <- list()
  index <- 1L
  for (endpoint_name in endpoint) {
    for (model_name in models) {
      fit <- if (model_name == "local") object$local_fits[[endpoint_name]] else
        object$endpoint_fits[[endpoint_name]]
      result <- .posterior_summary(fit, grid$x, probs)
      result$dose <- grid$dose
      result$endpoint <- endpoint_name
      result$model <- model_name
      pieces[[index]] <- result[, c("endpoint", "model", "dose", "x", "estimate", "lower", "upper")]
      index <- index + 1L
    }
  }
  result <- do.call(rbind, pieces)
  rownames(result) <- NULL
  result
}
