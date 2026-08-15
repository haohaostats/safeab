#' Fit a SAFE-AB model
#'
#' @param data A `safeab_data` object.
#' @param donors Optional character vector used for every fitted endpoint, or a
#'   named list defining endpoint-specific donor sets, for example
#'   `list(toxicity = c("A", "B"), efficacy = "A")`. By default, all
#'   non-target studies in `data` are used for every endpoint.
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
  if (is.null(endpoints)) endpoints <- data$endpoints
  endpoints <- unique(as.character(endpoints))
  unknown_endpoints <- setdiff(endpoints, data$endpoints)
  if (!length(endpoints) || length(unknown_endpoints)) {
    stop("endpoints must select available toxicity and/or efficacy outcomes.", call. = FALSE)
  }
  if (is.null(donors)) {
    donors_by_endpoint <- stats::setNames(
      rep(list(data$donors), length(endpoints)), endpoints
    )
  } else if (is.list(donors)) {
    if (is.null(names(donors)) || any(!nzchar(names(donors))) || anyDuplicated(names(donors))) {
      stop("A donor list must have unique endpoint names.", call. = FALSE)
    }
    valid_values <- vapply(donors, function(value) {
      is.null(value) || is.character(value)
    }, logical(1))
    if (!all(valid_values)) {
      stop("Every donor-list element must be a character vector or NULL.", call. = FALSE)
    }
    unknown_names <- setdiff(names(donors), data$endpoints)
    missing_names <- setdiff(endpoints, names(donors))
    if (length(unknown_names) || length(missing_names)) {
      stop(
        "A donor list must contain every fitted endpoint name and only available endpoint names. Missing: ",
        if (length(missing_names)) paste(missing_names, collapse = ", ") else "none",
        "; unknown: ",
        if (length(unknown_names)) paste(unknown_names, collapse = ", ") else "none",
        ".", call. = FALSE
      )
    }
    donors_by_endpoint <- lapply(donors[endpoints], function(value) {
      if (is.null(value)) character() else unique(as.character(value))
    })
  } else {
    common_donors <- unique(as.character(donors))
    donors_by_endpoint <- stats::setNames(
      rep(list(common_donors), length(endpoints)), endpoints
    )
  }
  donor_union <- unique(as.character(unlist(donors_by_endpoint, use.names = FALSE)))
  unknown <- setdiff(donor_union, data$donors)
  if (length(unknown)) {
    stop("Unknown or target donor identifiers: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  fits <- stats::setNames(lapply(endpoints, function(endpoint) {
    .fit_borrowed_endpoint(
      .endpoint_data(data$data, endpoint), data$target, donors_by_endpoint[[endpoint]],
      endpoint, prior, control
    )
  }), endpoints)
  structure(list(
    call = match.call(),
    data = data,
    target = data$target,
    donors = donor_union,
    donors_by_endpoint = donors_by_endpoint,
    endpoints = endpoints,
    endpoint_fits = lapply(fits, `[[`, "borrowed"),
    local_fits = lapply(fits, `[[`, "local_target"),
    donor_fits = lapply(fits, `[[`, "local_donors"),
    weights = lapply(fits, `[[`, "weights"),
    prior = prior,
    control = control
  ), class = "safeab_fit")
}
