test_that("safeab_data validates and standardizes aggregate data", {
  raw <- example_safeab_data()
  dat <- safeab_data(raw, target = "target")
  expect_s3_class(dat, "safeab_data")
  expect_equal(dat$target, "target")
  expect_setequal(dat$donors, c("donor_A", "donor_B"))
  expect_true(all(dat$data$x >= 0 & dat$data$x <= 1))
  expect_equal(dat$dose_transform$type, "rank")
  for (study_id in unique(dat$data$study)) {
    observed <- dat$data$x[dat$data$study == study_id]
    expect_equal(observed, seq(0, 1, length.out = length(observed)))
  }
  expect_error(
    safeab_data(transform(raw, toxicity = n_toxicity + 1), target = "target"),
    "0 <= events <= total"
  )
})

test_that("rank standardization is performed separately within each study", {
  raw <- example_safeab_data()
  raw$dose[raw$study == "donor_A"] <- c(10, 20, 40, 80)
  raw$dose[raw$study == "donor_B"] <- c(100, 200, 400, 800)
  dat <- safeab_data(raw, target = "target")
  expected <- seq(0, 1, length.out = 4)
  for (study_id in unique(dat$data$study)) {
    expect_equal(dat$data$x[dat$data$study == study_id], expected)
  }

  ranged <- safeab_data(raw, target = "target", dose_transform = "range")
  expect_false(isTRUE(all.equal(
    ranged$data$x[ranged$data$study == "target"], expected
  )))
})

test_that("rank standardization requires at least two doses per study", {
  raw <- example_safeab_data()
  raw <- raw[!(raw$study == "donor_B" & raw$dose != 1), ]
  expect_error(
    safeab_data(raw, target = "target"),
    "donor_B has only one"
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
