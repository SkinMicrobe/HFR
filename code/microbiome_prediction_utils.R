# Microbiome prediction helpers (Adonis, logistic, binning) from the prediction script.
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(tableone)
  library(stats)
})

# Adonis with FDR adjustment, returning the filtered table.
adonis_with_fdr <- function(feature_mat, metadata, fdr_cut = 1) {
  res <- adonis(feature_mat ~ ., data = metadata)
  tab <- as.data.frame(res$aov.tab)
  tab <- tab[seq_len(ncol(metadata)), , drop = FALSE]
  tab$fdr <- p.adjust(tab$`Pr(>F)`, method = "fdr")
  if (!is.null(fdr_cut)) {
    tab <- tab[tab$fdr < fdr_cut, , drop = FALSE]
  }
  tab
}

# Build TableOne for prediction variables.
prediction_table_one <- function(data, strata, cat_vars) {
  CreateTableOne(data = data, strata = strata, factorVars = cat_vars)
}

# Logistic regression with stepwise selection (both directions).
stepwise_logistic <- function(formula, data, family = "binomial", maxit = 200) {
  fit <- glm(formula = formula, data = data, family = family, maxit = maxit)
  step(fit, direction = "both", data = data)
}

# Forestplot-style data frame (mean, lower, upper) from a logistic model.
logistic_or_forest_df <- function(fit) {
  ors <- exp(coef(fit))
  cis <- exp(confint(fit))
  df <- data.frame(term = names(ors),
                   mean = ors,
                   min = cis[, 1],
                   max = cis[, 2],
                   row.names = NULL)
  df
}
