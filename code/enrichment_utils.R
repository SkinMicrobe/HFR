# Enrichment and GSEA helpers based on existing workflows.
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(stringr)
  library(enrichplot)
})

# KEGG enrichment from ranked list (GSEA style).
gsea_kegg_from_ranked <- function(gene_list,
                                  organism = "hsa",
                                  p_cut = 0.05,
                                  keywords = NULL,
                                  top_n = 6) {
  gsea_res <- gseKEGG(geneList = gene_list,
                      organism = organism,
                      pvalueCutoff = p_cut)
  if (is.null(keywords)) {
    return(list(result = gsea_res, plot_ids = head(gsea_res@result$ID, top_n)))
  }
  matched <- gsea_res@result %>%
    filter(str_detect(Description, paste(keywords, collapse = "|"))) %>%
    arrange(p.adjust) %>%
    slice(seq_len(min(top_n, n())))
  list(result = gsea_res, plot_ids = matched$ID)
}

# Extract top pathways table from enrichResult.
extract_top_pathways <- function(enrich_obj, n = 20) {
  as_tibble(enrich_obj@result) %>%
    arrange(p.adjust) %>%
    slice(seq_len(min(n, nrow(.))))
}
