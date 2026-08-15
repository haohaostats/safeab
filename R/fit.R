#' Fit a SAFE-AB model
#'
#' @param data A `safeab_data` object.
#' @param donors Optional character vector selecting donor studies. By default,
#'   all non-target studies in `data` are used.
#' @param endpoints Endpoints to fit. The default fits every available endpoint.
#' @param prior A `safeab_prior` object.
#' @param control A `safeab_control` object.
#' @return A fitted object of class `safeab_fit`.
#' @export
fit_safeab <- function(
    data,
    donors = NULL,
    endpoints = NULL,
    prior = safeab_prior(),
    control = safeab_control()) {
  if (!inherits(data, "safeab_data")) {
    stop("data must be created by safeab_data().", call. = FALSE)
  }
  if (!inherits(prior, "safeab_prior")) stop("prior must be a safeab_prior.", call. = FALSE)
  if (!inherits(control, "safeab_control")) stop("control must be a safeab_control.", call. = FALSE)
  if (is.null(donors)) donors <- data$donors
  donors <- unique(as.character(donors))
  unknown <- setdiff(donors, data$donors)
  if (length(unknown)) {
    stop("Unknown or target donor identifiers: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (is.null(endpoints)) endpoints <- data$endpoints
  endpoints <- unique(as.character(endpoints))
  unknown_endpoints <- setdiff(endpoints, data$endpoints)
  if (!length(endpoints) || length(unknown_endpoints)) {
    stop("endpoints must select available toxicity and/or efficacy outcomes.", call. = FALSE)
  }
  fits <- stats::setNames(lapply(endpoints, function(endpoint) {
    .fit_borrowed_endpoint(
      .endpoint_data(data$data, endpoint), data$target, donors,
      endpoint, prior, control
    )
  }), endpoints)
  structure(list(
    call = match.call(),
    data = data,
    target = data$target,
    donors = donors,
    endpoints = endpoints,
    endpoint_fits = lapply(fits, `[[`, "borrowed"),
    local_fits = lapply(fits, `[[`, "local_target"),
    donor_fits = lapply(fits, `[[`, "local_donors"),
    weights = lapply(fits, `[[`, "weights"),
    prior = prior,
    control = control
  ), class = "safeab_fit")
}
