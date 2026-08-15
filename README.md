# safeab

`safeab` implements safety-aware selective Bayesian borrowing for a target
dose-response study with one or more external donor studies. It is an analysis
package for new datasets, not a reproduction package for the SAFE-AB paper.

## Installation

Install the development release from GitHub:

```r
# install.packages("remotes")
remotes::install_github("haohaostats/safeab")
```

## Basic workflow

The package provides a compact workflow:

```r
dat <- safeab_data(
  trial_data,
  target = "target",
  study = "study",
  dose = "dose",
  toxicity_events = "toxicity",
  toxicity_total = "n_toxicity",
  efficacy_events = "efficacy",
  efficacy_total = "n_efficacy"
)

fit <- fit_safeab(dat)
summary(fit)
borrowing_summary(fit)
predict(fit, newdata = seq(0, 1, length.out = 21))
recommend_dose(fit, toxicity_limit = 0.30)
```

The initial release supports aggregated binary toxicity and efficacy outcomes,
endpoint-specific borrowing, donor and total BESS caps, posterior prediction,
and a target-only local safety veto.
