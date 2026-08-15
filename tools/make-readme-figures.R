.libPaths(c(normalizePath(file.path("..", ".Rlib")), .libPaths()))

library(safeab)

dir.create(file.path("man", "figures"), recursive = TRUE, showWarnings = FALSE)

raw <- example_safeab_data()
dat <- safeab_data(raw, target = "target")
fit <- fit_safeab(dat)

ink <- "#263238"
muted <- "#6B7280"
blue <- "#2563A6"
blue_fill <- grDevices::adjustcolor(blue, alpha.f = 0.16)
orange <- "#D97706"
green <- "#0F766E"
pink <- "#B94E78"
red <- "#C43D3D"

open_png <- function(filename, width = 1600, height = 720) {
  grDevices::png(
    file.path("man", "figures", filename),
    width = width, height = height, res = 180, type = "cairo-png",
    bg = "white"
  )
}

plot_interval <- function(x, estimate, lower, upper, color, pch) {
  graphics::segments(x, lower, x, upper, col = color, lwd = 1.7)
  graphics::segments(x - 0.07, lower, x + 0.07, lower, col = color, lwd = 1.4)
  graphics::segments(x - 0.07, upper, x + 0.07, upper, col = color, lwd = 1.4)
  graphics::points(x, estimate, pch = pch, bg = "white", col = color, cex = 1.25, lwd = 1.6)
}

open_png("posterior-dose-response.png")
graphics::par(
  mfrow = c(1, 2), mar = c(4.1, 4.2, 2.5, 0.8), oma = c(0, 0, 1.5, 0),
  family = "sans", las = 1, bty = "l", col.axis = ink, col.lab = ink,
  cex.axis = 0.86, cex.lab = 0.92
)
target <- dat$data[dat$data$study == dat$target, , drop = FALSE]
grid <- seq(0, 1, length.out = 151)
for (endpoint in c("toxicity", "efficacy")) {
  pred <- predict(fit, grid, endpoint = endpoint, model = "both", dose_scale = "standardized")
  local <- pred[pred$model == "local", ]
  borrowed <- pred[pred$model == "borrowed", ]
  observed <- if (endpoint == "toxicity") {
    target$toxicity / target$n_toxicity
  } else {
    target$efficacy / target$n_efficacy
  }
  graphics::plot(
    borrowed$x, borrowed$estimate, type = "n", ylim = c(0, 0.72),
    xlab = "Standardized dose", ylab = "Posterior event probability",
    main = if (endpoint == "toxicity") "A  Toxicity" else "B  Efficacy"
  )
  graphics::polygon(
    c(borrowed$x, rev(borrowed$x)), c(borrowed$lower, rev(borrowed$upper)),
    border = NA, col = blue_fill
  )
  graphics::lines(local$x, local$estimate, col = muted, lwd = 1.8, lty = 2)
  graphics::lines(borrowed$x, borrowed$estimate, col = blue, lwd = 2.4)
  graphics::points(target$x, observed, pch = 21, bg = "white", col = ink, cex = 1.05, lwd = 1.3)
  if (endpoint == "toxicity") {
    graphics::legend(
      "topleft", c("SAFE-AB", "95% credible interval", "LOCAL", "Target observed"),
      col = c(blue, blue, muted, ink), lty = c(1, NA, 2, NA), lwd = c(2.4, NA, 1.8, NA),
      pch = c(NA, 15, NA, 21), pt.bg = c(NA, blue_fill, NA, "white"),
      pt.cex = c(NA, 1.4, NA, 1), bty = "n", cex = 0.74
    )
  }
}
graphics::mtext("Posterior dose-response estimates", outer = TRUE, font = 2, cex = 1.05)
grDevices::dev.off()

borrow <- borrowing_summary(fit)
borrow$label <- paste(
  gsub("_", " ", sub("donor_", "Donor ", borrow$donor, fixed = TRUE)),
  tools::toTitleCase(borrow$endpoint), sep = "  -  "
)
borrow <- borrow[order(borrow$relative_weight), ]
y <- seq_len(nrow(borrow))
point_color <- ifelse(borrow$endpoint == "toxicity", orange, green)
point_shape <- ifelse(borrow$endpoint == "toxicity", 21, 22)

open_png("borrowing-diagnostics.png", height = 650)
graphics::par(
  mar = c(4.3, 7.4, 3.2, 1.5), family = "sans", las = 1, bty = "n",
  col.axis = ink, col.lab = ink, cex.axis = 0.86, cex.lab = 0.92
)
graphics::plot(
  borrow$relative_weight, y, type = "n", xlim = c(0, 0.28), ylim = c(0.5, 4.7),
  yaxt = "n", xlab = "Borrowed fraction of available donor information", ylab = "",
  main = "Selective borrowing by donor and endpoint"
)
graphics::axis(2, at = y, labels = borrow$label, tick = FALSE, las = 1, col.axis = ink)
graphics::segments(0, y, borrow$relative_weight, y, col = grDevices::adjustcolor(point_color, 0.45), lwd = 2)
for (i in seq_len(nrow(borrow))) {
  graphics::points(
    borrow$relative_weight[i], y[i], pch = point_shape[i],
    bg = point_color[i], col = point_color[i], cex = 1.25
  )
  graphics::text(
    borrow$relative_weight[i] + 0.012, y[i],
    labels = sprintf("BESS %.1f", borrow$borrowed_n[i]),
    pos = 4, cex = 0.78, col = ink
  )
}
graphics::legend(
  "bottomright", c("Toxicity", "Efficacy"), pch = c(21, 22),
  pt.bg = c(orange, green), col = c(orange, green), bty = "n", cex = 0.82
)
graphics::mtext(
  "More compatible endpoint-specific donors receive more weight; BESS caps remain enforced.",
  side = 3, line = 0.35, cex = 0.75, col = muted
)
grDevices::dev.off()

decision <- recommend_dose(fit)
tab <- decision$dose_table

open_png("dose-recommendation.png", height = 700)
graphics::par(
  mar = c(4.4, 4.4, 3.2, 1.5), family = "sans", las = 1, bty = "l",
  col.axis = ink, col.lab = ink, cex.axis = 0.86, cex.lab = 0.92
)
graphics::plot(
  tab$dose, tab$toxicity, type = "n", ylim = c(0, 0.72),
  xlab = "Candidate dose", ylab = "Posterior event probability",
  xaxt = "n", main = "Safety-constrained dose recommendation"
)
graphics::axis(1, at = tab$dose, labels = tab$dose)
graphics::abline(h = 0.30, col = muted, lty = 3, lwd = 1.5)
graphics::text(6.7, 0.315, "Toxicity limit", pos = 2, cex = 0.72, col = muted)
plot_interval(tab$dose - 0.12, tab$toxicity, tab$toxicity_lower, tab$toxicity_upper, orange, 21)
plot_interval(tab$dose + 0.12, tab$efficacy, tab$efficacy_lower, tab$efficacy_upper, blue, 22)
graphics::points(
  tab$dose, rep(0.025, nrow(tab)),
  pch = ifelse(tab$admissible, 16, 4),
  col = ifelse(tab$admissible, green, pink), cex = 1.05, lwd = 1.7
)
selected <- decision$selected_row
graphics::arrows(
  tab$dose[selected] + 1.05, 0.09,
  tab$dose[selected] + 0.08, 0.032,
  length = 0.10, angle = 24, code = 2,
  col = red, lwd = 1.8
)
graphics::text(
  tab$dose[selected] + 1.15, 0.09, "Recommended", pos = 4,
  col = red, font = 2, cex = 0.79
)
graphics::legend(
  "topleft",
  c("Toxicity", "Efficacy", "Admissible", "Excluded"),
  pch = c(21, 22, 16, 4), pt.bg = c("white", "white", NA, NA),
  col = c(orange, blue, green, pink), bty = "n", cex = 0.78, ncol = 2
)
graphics::mtext(
  "Intervals show 95% posterior credible intervals; the local safety veto excludes dose 8.",
  side = 3, line = 0.35, cex = 0.75, col = muted
)
grDevices::dev.off()
