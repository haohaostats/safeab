.safeab_plot_colors <- function() {
  c(
    toxicity = "#D97706",
    efficacy = "#2563A6",
    admissible = "#0F766E",
    excluded = "#B94E78",
    recommendation = "#C43D3D",
    neutral = "#6B7280"
  )
}

.plot_probability_interval <- function(x, estimate, lower, upper, color, pch, offset) {
  graphics::segments(x, lower, x, upper, col = color, lwd = 1.7)
  graphics::segments(x - offset, lower, x + offset, lower, col = color, lwd = 1.3)
  graphics::segments(x - offset, upper, x + offset, upper, col = color, lwd = 1.3)
  graphics::points(x, estimate, pch = pch, bg = "white", col = color, cex = 1.15, lwd = 1.5)
}

#' Plot selective borrowing diagnostics
#'
#' @param object A fitted `safeab_fit` object.
#' @param metric Display the borrowed fraction of available donor information or
#'   the borrowed effective sample size.
#' @param endpoint Optional endpoint subset.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#' @return The plotted borrowing summary, invisibly.
#' @export
plot_borrowing <- function(
    object,
    metric = c("relative_weight", "borrowed_n"),
    endpoint = NULL,
    ...) {
  if (!inherits(object, "safeab_fit")) stop("object must be a safeab_fit.", call. = FALSE)
  metric <- match.arg(metric)
  result <- borrowing_summary(object)
  if (is.null(result) || !nrow(result)) {
    stop("No donor borrowing is available to plot.", call. = FALSE)
  }
  if (!is.null(endpoint)) {
    endpoint <- match.arg(endpoint, object$endpoints, several.ok = TRUE)
    result <- result[result$endpoint %in% endpoint, , drop = FALSE]
  }
  value <- result[[metric]]
  ordering <- order(value)
  result <- result[ordering, , drop = FALSE]
  value <- value[ordering]
  labels <- paste(
    tools::toTitleCase(gsub("_", " ", result$donor, fixed = TRUE)),
    tools::toTitleCase(result$endpoint), sep = "  -  "
  )
  colors <- .safeab_plot_colors()[result$endpoint]
  shapes <- ifelse(result$endpoint == "toxicity", 21, 22)
  y <- seq_len(nrow(result))
  upper <- max(value) * 1.18
  if (!is.finite(upper) || upper <= 0) upper <- 1
  graphics::plot(
    value, y, type = "n", xlim = c(0, upper), ylim = c(0.5, nrow(result) + 0.5),
    yaxt = "n", ylab = "",
    xlab = if (metric == "relative_weight") {
      "Borrowed fraction of available donor information"
    } else {
      "Borrowed effective sample size"
    },
    main = "Selective borrowing by donor and endpoint", bty = "n", ...
  )
  graphics::axis(2, at = y, labels = labels, tick = FALSE, las = 1)
  graphics::segments(0, y, value, y, col = grDevices::adjustcolor(colors, 0.45), lwd = 2)
  for (i in seq_len(nrow(result))) {
    graphics::points(value[i], y[i], pch = shapes[i], bg = colors[i], col = colors[i], cex = 1.2)
    label <- if (metric == "relative_weight") {
      paste0("BESS ", format(round(result$borrowed_n[i], 1), nsmall = 1))
    } else {
      paste0(format(round(100 * result$relative_weight[i], 1), nsmall = 1), "%")
    }
    graphics::text(value[i] + 0.025 * upper, y[i], label, pos = 4, cex = 0.78)
  }
  graphics::legend(
    "bottomright", c("Toxicity", "Efficacy"), pch = c(21, 22),
    pt.bg = .safeab_plot_colors()[c("toxicity", "efficacy")],
    col = .safeab_plot_colors()[c("toxicity", "efficacy")], bty = "n"
  )
  invisible(result)
}

#' Plot a SAFE-AB dose recommendation
#'
#' @param x A `safeab_recommendation` object.
#' @param show_intervals Show 95 percent posterior credible intervals.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#' @return The recommendation object, invisibly.
#' @export
plot.safeab_recommendation <- function(x, show_intervals = TRUE, ...) {
  table <- x$dose_table
  use_original_dose <- all(is.finite(table$dose))
  dose <- if (use_original_dose) table$dose else table$x
  axis_label <- if (use_original_dose) "Candidate dose" else "Standardized dose"
  dose_span <- diff(range(dose))
  if (!is.finite(dose_span) || dose_span <= 0) dose_span <- 1
  offset <- 0.012 * dose_span
  colors <- .safeab_plot_colors()
  upper_values <- c(table$toxicity_upper, table$efficacy_upper, x$settings$toxicity_limit)
  y_max <- min(1, max(upper_values, na.rm = TRUE) * 1.10)
  y_max <- max(y_max, 0.55)
  graphics::plot(
    dose, table$toxicity, type = "n", ylim = c(0, y_max), xaxt = "n",
    xlab = axis_label, ylab = "Posterior event probability",
    main = "Safety-constrained dose recommendation", ...
  )
  graphics::axis(1, at = dose, labels = format(dose, trim = TRUE))
  graphics::abline(h = x$settings$toxicity_limit, col = colors["neutral"], lty = 3, lwd = 1.4)
  graphics::text(
    max(dose) - 0.20 * dose_span, x$settings$toxicity_limit + 0.02 * y_max,
    "Toxicity limit", col = colors["neutral"], cex = 0.78
  )
  if (isTRUE(show_intervals)) {
    .plot_probability_interval(
      dose - offset, table$toxicity, table$toxicity_lower, table$toxicity_upper,
      colors["toxicity"], 21, offset * 0.55
    )
    .plot_probability_interval(
      dose + offset, table$efficacy, table$efficacy_lower, table$efficacy_upper,
      colors["efficacy"], 22, offset * 0.55
    )
  } else {
    graphics::points(dose - offset, table$toxicity, pch = 21, bg = "white", col = colors["toxicity"])
    graphics::points(dose + offset, table$efficacy, pch = 22, bg = "white", col = colors["efficacy"])
  }
  status_y <- 0.025 * y_max
  graphics::points(
    dose, rep(status_y, length(dose)),
    pch = ifelse(table$admissible, 16, 4),
    col = ifelse(table$admissible, colors["admissible"], colors["excluded"]),
    cex = 1.05, lwd = 1.6
  )
  selected <- x$selected_row
  if (!is.na(selected)) {
    direction <- if (dose[selected] <= mean(range(dose))) 1 else -1
    label_x <- dose[selected] + direction * 0.28 * dose_span
    label_y <- 0.12 * y_max
    graphics::arrows(
      label_x - direction * 0.025 * dose_span, label_y,
      dose[selected] + direction * 0.012 * dose_span, status_y * 1.25,
      length = 0.10, angle = 24, code = 2,
      col = colors["recommendation"], lwd = 1.8
    )
    graphics::text(
      label_x, label_y, "Recommended",
      pos = if (direction > 0) 4 else 2,
      col = colors["recommendation"], font = 2, cex = 0.82
    )
  }
  graphics::legend(
    "topleft", c("Toxicity", "Efficacy", "Admissible", "Excluded"),
    pch = c(21, 22, 16, 4), pt.bg = c("white", "white", NA, NA),
    col = colors[c("toxicity", "efficacy", "admissible", "excluded")],
    bty = "n", ncol = 2
  )
  invisible(x)
}
