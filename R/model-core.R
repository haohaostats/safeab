.log1pexp <- function(x) {
  out <- numeric(length(x))
  high <- x > 0
  out[high] <- x[high] + log1p(exp(-x[high]))
  out[!high] <- log1p(exp(x[!high]))
  out
}

.model_eta <- function(theta, x, endpoint) {
  if (endpoint == "toxicity") {
    theta[1L] + exp(theta[2L]) * x
  } else {
    z <- x - 0.5
    theta[1L] + theta[2L] * z + theta[3L] * z^2
  }
}

.model_design <- function(theta, x, endpoint) {
  if (endpoint == "toxicity") {
    cbind(1, exp(theta[2L]) * x)
  } else {
    z <- x - 0.5
    cbind(1, z, z^2)
  }
}

.make_logposterior <- function(data, endpoint, prior, cell_weight = NULL) {
  if (is.null(cell_weight)) cell_weight <- rep(1, nrow(data))
  x <- data$x
  y <- data$y
  n <- data$n
  w <- as.numeric(cell_weight)
  mu <- prior$mean
  prior_sd <- prior$sd
  fn <- function(theta) {
    eta <- .model_eta(theta, x, endpoint)
    sum(w * (y * eta - n * .log1pexp(eta))) -
      0.5 * sum(((theta - mu) / prior_sd)^2)
  }
  gr <- function(theta) {
    eta <- .model_eta(theta, x, endpoint)
    p <- stats::plogis(eta)
    design <- .model_design(theta, x, endpoint)
    as.numeric(crossprod(design, w * (y - n * p))) -
      (theta - mu) / prior_sd^2
  }
  he <- function(theta) {
    eta <- .model_eta(theta, x, endpoint)
    p <- stats::plogis(eta)
    design <- .model_design(theta, x, endpoint)
    hessian <- -crossprod(design, design * (w * n * p * (1 - p)))
    if (endpoint == "toxicity") {
      hessian[2L, 2L] <- hessian[2L, 2L] +
        sum(w * (y - n * p) * exp(theta[2L]) * x)
    }
    hessian - diag(1 / prior_sd^2, nrow = length(theta))
  }
  list(fn = fn, gr = gr, he = he)
}

.posterior_mode <- function(functions, starts) {
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    tryCatch(stats::optim(
      starts[i, ],
      fn = function(theta) -functions$fn(theta),
      gr = function(theta) -functions$gr(theta),
      method = "BFGS",
      control = list(maxit = 1000, reltol = 1e-10)
    ), error = function(e) NULL)
  })
  valid <- vapply(fits, function(fit) {
    !is.null(fit) && is.finite(fit$value) && fit$convergence == 0
  }, logical(1))
  if (!any(valid)) stop("Posterior optimization failed from every start.", call. = FALSE)
  fits <- fits[valid]
  fits[[which.min(vapply(fits, function(fit) fit$value, numeric(1)))]]$par
}

.fit_endpoint <- function(data, endpoint, prior, control, cell_weight = NULL) {
  functions <- .make_logposterior(data, endpoint, prior[[endpoint]], cell_weight)
  start <- .posterior_mode(functions, .safeab_starts(endpoint))
  quadrature <- aghq::aghq(
    ff = functions,
    k = control$quadrature_points,
    startingvalue = start,
    control = aghq::default_control(method = "BFGS")
  )
  nodes <- aghq::get_nodesandweights(quadrature)
  parameter_names <- paste0("theta", seq_along(prior[[endpoint]]$mean))
  mass <- nodes$weights * exp(nodes$logpost_normalized)
  mass <- mass / sum(mass)
  structure(list(
    endpoint = endpoint,
    data = data,
    cell_weight = if (is.null(cell_weight)) rep(1, nrow(data)) else cell_weight,
    theta = as.matrix(nodes[, parameter_names, drop = FALSE]),
    mass = mass,
    mode = aghq::get_mode(quadrature),
    quadrature = quadrature
  ), class = "safeab_endpoint_fit")
}

.posterior_grid <- function(fit, x) {
  eta <- vapply(seq_len(nrow(fit$theta)), function(i) {
    .model_eta(fit$theta[i, ], x, fit$endpoint)
  }, numeric(length(x)))
  if (length(x) == 1L) eta <- matrix(eta, nrow = 1L)
  list(eta = eta, probability = stats::plogis(eta), mass = fit$mass)
}

.weighted_quantile <- function(x, weights, probabilities) {
  ordering <- order(x)
  x <- x[ordering]
  weights <- weights[ordering] / sum(weights)
  cumulative <- cumsum(weights)
  vapply(probabilities, function(probability) {
    x[which(cumulative >= probability)[1L]]
  }, numeric(1))
}

.posterior_summary <- function(fit, x, probs = c(0.025, 0.975), threshold = NULL) {
  grid <- .posterior_grid(fit, x)
  mean_probability <- as.numeric(grid$probability %*% grid$mass)
  intervals <- t(vapply(seq_along(x), function(i) {
    .weighted_quantile(grid$probability[i, ], grid$mass, probs)
  }, numeric(length(probs))))
  result <- data.frame(
    x = x,
    estimate = mean_probability,
    lower = intervals[, 1L],
    upper = intervals[, 2L]
  )
  if (!is.null(threshold)) {
    result$probability_exceeding <- as.numeric((grid$probability > threshold) %*% grid$mass)
  }
  result
}

.linear_predictor_moments <- function(fit, x) {
  grid <- .posterior_grid(fit, x)
  mean <- as.numeric(grid$eta %*% grid$mass)
  second <- as.numeric((grid$eta^2) %*% grid$mass)
  list(mean = mean, variance = pmax(0, second - mean^2))
}
