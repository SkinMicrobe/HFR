# KEGG enrichment for human SYMBOL lists.
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

kegg_enrich_symbols <- function(symbols, p_cut = 0.05, q_cut = 0.2, show = 20) {
  genes <- clusterProfiler::bitr(symbols,
                                 fromType = "SYMBOL",
                                 toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db)
  ekk <- enrichKEGG(gene = genes$ENTREZID,
                    organism = "hsa",
                    pvalueCutoff = p_cut,
                    qvalueCutoff = q_cut)
  list(
    result = ekk@result,
    plot = dotplot(ekk, showCategory = show) +
      scale_y_discrete(labels = function(x) gsub("\n", " ", x))
  )
}
