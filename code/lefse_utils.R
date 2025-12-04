# LEfSe (microeco) helpers reflecting the LDA workflows.
suppressPackageStartupMessages({
  library(microeco)
  library(ggsci)
})

# Run LEfSe via microeco::trans_diff and return results plus a bar plot.
run_lefse <- function(feature_table, sample_table, tax_table,
                      group_col = "group", alpha = 0.05,
                      lda_cutoff = 2.0,
                      bar_use_number = 1:50,
                      group_order = NULL) {
  dataset <- microtable$new(sample_table = sample_table,
                            otu_table = feature_table,
                            tax_table = tax_table)
  lefse <- trans_diff$new(dataset = dataset,
                          method = "lefse",
                          group = group_col,
                          alpha = alpha,
                          lda_cutoff = lda_cutoff,
                          lefse_subgroup = NULL)
  plot_bar <- lefse$plot_diff_bar(use_number = bar_use_number,
                                  width = 0.3,
                                  group_order = group_order) +
    ggsci::scale_color_npg() +
    ggsci::scale_fill_npg()
  list(result = lefse$res_diff, plot = plot_bar)
}

# Optional cladogram helper.
plot_lefse_cladogram <- function(lefse_obj,
                                 taxa_num = 200,
                                 feature_num = 50,
                                 clade_label_level = 5,
                                 group_order = NULL,
                                 select_labels = NULL) {
  lefse_obj$plot_diff_cladogram(use_taxa_num = taxa_num,
                                use_feature_num = feature_num,
                                clade_label_level = clade_label_level,
                                group_order = group_order,
                                select_show_labels = select_labels)
}
