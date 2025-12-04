# Limma differential expression helper.
suppressPackageStartupMessages({
  library(limma)
  library(tibble)
})

# Run limma with a contrast formula and return a tibble of results.
run_limma <- function(expr_mat, design, contrast_formula, coef = 1) {
  stopifnot(all(rownames(design) == colnames(expr_mat)))
  fit <- lmFit(expr_mat, design)
  cm <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  topTable(fit2, coef = coef, number = Inf) |>
    tibble::rownames_to_column("gene")
}
