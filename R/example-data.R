#' Generate a small synthetic SAFE-AB example
#'
#' This deterministic dataset is intended for examples and interface checks.
#'
#' @return A data frame with one target and two donor studies.
#' @export
example_safeab_data <- function() {
  dose <- rep(c(1, 2, 4, 8), 3)
  study <- rep(c("target", "donor_A", "donor_B"), each = 4)
  data.frame(
    study = study,
    dose = dose,
    toxicity = c(0, 1, 3, 6, 0, 1, 2, 5, 1, 3, 6, 8),
    n_toxicity = rep(20, 12),
    efficacy = c(2, 5, 9, 10, 3, 6, 10, 11, 1, 3, 6, 8),
    n_efficacy = rep(20, 12),
    stringsAsFactors = FALSE
  )
}
