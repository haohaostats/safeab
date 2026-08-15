test_that("safeab_data validates and standardizes aggregate data", {
  raw <- example_safeab_data()
  dat <- safeab_data(raw, target = "target")
  expect_s3_class(dat, "safeab_data")
  expect_equal(dat$target, "target")
  expect_setequal(dat$donors, c("donor_A", "donor_B"))
  expect_true(all(dat$data$x >= 0 & dat$data$x <= 1))
  expect_error(
    safeab_data(transform(raw, toxicity = n_toxicity + 1), target = "target"),
    "0 <= events <= total"
  )
})

test_that("endpoint columns can be omitted", {
  raw <- example_safeab_data()
  dat <- safeab_data(
    raw, target = "target",
    efficacy_events = NULL, efficacy_total = NULL
  )
  expect_equal(dat$endpoints, "toxicity")
})
