# IOBR scoring helpers (wrapping calculate_sig_score and count2tpm patterns).
suppressPackageStartupMessages({
  library(IOBR)
})

# Calculate signature scores with a chosen method (pca, ssgsea, etc.).
iobr_signature_score <- function(eset, signature, method = "pca", mini_gene_count = 3) {
  calculate_sig_score(eset = eset,
                      signature = signature,
                      method = method,
                      mini_gene_count = mini_gene_count)
}

# Convert counts to TPM with optional organism/id type settings.
iobr_count_to_tpm <- function(count_mat,
                              id_type = "SYMBOL",
                              org = "hsa",
                              source = "web",
                              eff_length = NULL,
                              id_col = "id",
                              symbol_col = "symbol",
                              length_col = "eff_length",
                              remove_redundancy = "mean") {
  count2tpm(countMat = count_mat,
            idType = id_type,
            org = org,
            source = source,
            effLength = eff_length,
            id = id_col,
            gene_symbol = symbol_col,
            length = length_col,
            remove_redundancy = remove_redundancy)
}
