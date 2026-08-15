test_that("SAFE-AB fits, predicts, and respects BESS caps", {
  skip_if_not_installed("aghq")
  dat <- safeab_data(example_safeab_data(), target = "target")
  control <- safeab_control(
    donor_cap = 0.20, total_cap = 0.30, quadrature_points = 3
  )
  fit <- fit_safeab(dat, control = control)
  expect_s3_class(fit, "safeab_fit")
  predictions <- predict(fit, endpoint = "toxicity", model = "both")
  expect_true(all(predictions$estimate > 0 & predictions$estimate < 1))
  expect_true(all(predictions$lower <= predictions$estimate))
  expect_true(all(predictions$estimate <= predictions$upper))

  bess <- effective_sample_size(fit)
  target <- dat$data[dat$data$study == dat$target, ]
  expect_lte(unname(bess["toxicity"]), 0.30 * sum(target$n_toxicity) + 1e-8)
  expect_lte(unname(bess["efficacy"]), 0.30 * sum(target$n_efficacy) + 1e-8)

  by_donor <- borrowing_summary(fit)
  expect_true(all(by_donor$borrowed_n <= 0.20 * 80 + 1e-8))
})

test_that("the local veto has the non-override property", {
  skip_if_not_installed("aghq")
  dat <- safeab_data(example_safeab_data(), target = "target")
  fit <- fit_safeab(dat, control = safeab_control(quadrature_points = 3))
  safety <- admissible_doses(
    fit, toxicity_limit = 0.15, borrowed_cutoff = 0.50,
    local_cutoff = 0.10, local_veto = TRUE
  )
  expect_false(any(safety$admissible[safety$local_exceedance > 0.10]))
  recommendation <- recommend_dose(fit)
  expect_s3_class(recommendation, "safeab_recommendation")
  expect_equal(nrow(recommendation$dose_table), 4)

  image_file <- tempfile(fileext = ".png")
  grDevices::png(image_file)
  expect_no_error(plot(fit))
  expect_no_error(plot_borrowing(fit))
  expect_no_error(plot(recommendation))
  grDevices::dev.off()
  expect_true(file.exists(image_file))
})

test_that("target-only fitting gives zero borrowed sample size", {
  skip_if_not_installed("aghq")
  dat <- safeab_data(example_safeab_data(), target = "target")
  fit <- fit_safeab(
    dat, donors = character(), endpoints = "toxicity",
    control = safeab_control(quadrature_points = 3)
  )
  expect_equal(unname(effective_sample_size(fit)), 0)
})
