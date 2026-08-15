#' @export
print.safeab_fit <- function(x, ...) {
  cat("SAFE-AB fit\n")
  cat("  Target:", x$target, "\n")
  cat("  Donors:", if (length(x$donors)) paste(x$donors, collapse = ", ") else "none", "\n")
  cat("  Endpoints:", paste(x$endpoints, collapse = ", "), "\n")
  bess <- effective_sample_size(x)
  cat("  Borrowed effective sample size:",
      paste(paste0(names(bess), "=", format(round(bess, 2), nsmall = 2)), collapse = ", "), "\n")
  invisible(x)
}

#' Summarize a SAFE-AB fit
#'
#' @param object A fitted `safeab_fit` object.
#' @param ... Unused.
#' @return A `summary.safeab_fit` object.
#' @export
summary.safeab_fit <- function(object, ...) {
  target_data <- object$data$data[object$data$data$study == object$target, , drop = FALSE]
  target_n <- vapply(object$endpoints, function(endpoint) {
    if (endpoint == "toxicity") sum(target_data$n_toxicity) else sum(target_data$n_efficacy)
  }, numeric(1))
  bess <- effective_sample_size(object)
  result <- list(
    call = object$call,
    target = object$target,
    donors = object$donors,
    endpoints = object$endpoints,
    endpoint_summary = data.frame(
      endpoint = object$endpoints,
      target_n = unname(target_n),
      borrowed_n = unname(bess),
      relative_bess = unname(bess / target_n)
    ),
    donor_summary = borrowing_summary(object)
  )
  class(result) <- "summary.safeab_fit"
  result
}

#' @export
print.summary.safeab_fit <- function(x, ...) {
  cat("SAFE-AB model summary\n")
  cat("Target:", x$target, "\n")
  cat("Donors:", if (length(x$donors)) paste(x$donors, collapse = ", ") else "none", "\n\n")
  print(x$endpoint_summary, row.names = FALSE, digits = 3)
  if (!is.null(x$donor_summary) && nrow(x$donor_summary)) {
    cat("\nBorrowing by donor:\n")
    print(x$donor_summary, row.names = FALSE, digits = 3)
  }
  invisible(x)
}

#' Plot posterior SAFE-AB dose-response curves
#'
#' @param x A fitted `safeab_fit` object.
#' @param endpoint Endpoints to plot.
#' @param n_grid Number of standardized dose-grid points.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#' @return The fitted object, invisibly.
#' @export
plot.safeab_fit <- function(x, endpoint = x$endpoints, n_grid = 101L, ...) {
  endpoint <- intersect(as.character(endpoint), x$endpoints)
  if (!length(endpoint)) stop("No fitted endpoint selected.", call. = FALSE)
  n_grid <- as.integer(n_grid)
  if (is.na(n_grid) || n_grid < 10L) stop("n_grid must be at least 10.", call. = FALSE)
  old <- graphics::par(mfrow = c(1, length(endpoint)))
  on.exit(graphics::par(old), add = TRUE)
  grid <- seq(0, 1, length.out = n_grid)
  target <- x$data$data[x$data$data$study == x$target, , drop = FALSE]
  for (endpoint_name in endpoint) {
    predictions <- stats::predict(
      x, grid, endpoint_name, model = "both", dose_scale = "standardized"
    )
    local <- predictions[predictions$model == "local", ]
    borrowed <- predictions[predictions$model == "borrowed", ]
    graphics::plot(
      borrowed$x, borrowed$estimate, type = "n", ylim = c(0, 1),
      xlab = "Standardized dose", ylab = paste("Posterior", endpoint_name, "probability"),
      ...
    )
    graphics::polygon(
      c(borrowed$x, rev(borrowed$x)), c(borrowed$lower, rev(borrowed$upper)),
      border = NA, col = grDevices::adjustcolor("#2C7FB8", alpha.f = 0.16)
    )
    graphics::lines(local$x, local$estimate, lty = 2, lwd = 1.5, col = "#666666")
    graphics::lines(borrowed$x, borrowed$estimate, lwd = 2, col = "#2C7FB8")
    observed <- if (endpoint_name == "toxicity") {
      target$toxicity / target$n_toxicity
    } else {
      target$efficacy / target$n_efficacy
    }
    graphics::points(target$x, observed, pch = 21, bg = "white", col = "black")
    graphics::legend(
      "topleft", c("SAFE-AB", "LOCAL", "Target observed"),
      col = c("#2C7FB8", "#666666", "black"),
      lty = c(1, 2, NA), pch = c(NA, NA, 21), pt.bg = "white", bty = "n"
    )
  }
  invisible(x)
}
