# Table and logistic helpers based on existing clinical analyses.
suppressPackageStartupMessages({
  library(tableone)
  library(survival)
})

# Build TableOne by strata.
build_table_one <- function(data, strata, cat_vars = NULL) {
  CreateTableOne(data = data, strata = strata, factorVars = cat_vars)
}

# Logistic regression with OR/CI output.
logistic_or_table <- function(formula, data, family = "binomial") {
  fit <- glm(formula = formula, data = data, family = family)
  ors <- exp(coef(fit))
  cis <- exp(confint(fit))
  list(model = fit, OR = ors, CI = cis)
}
