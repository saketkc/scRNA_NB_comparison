set.seed(42)
suppressMessages({
  library(tidyverse)
  library(iCOBRA)
})


`%notin%` <- Negate(`%in%`)





perf_datasets <- list()
all_fdr_tpr <- list()
all_fdrtprcurve <- list()
markers_list_df <- readRDS(here::here("output/muscat_simulation/results/muscat__DE_markers_list_df.rds"))

seeds <- c(20141015, 20160823, 20210827)
for (seed in seeds) {
  for (test.use in c("wilcox", "MAST")) {
    for (de_percent in c(0.01, 0.05, 0.1, 0.2)) {
      print(paste(seed, test.use, de_percent))
      df_shortlist <- markers_list_df[markers_list_df$seed == seed, ]
      df_shortlist <- df_shortlist[df_shortlist$de_percent == de_percent, ]
      df_shortlist <- df_shortlist[df_shortlist$test.use == test.use, ]
      p_vals <- list()
      p_adjs <- list()
      truths <- list()
      for (method in unique(df_shortlist$method)) {
        df_method <- df_shortlist[df_shortlist$method == method, ]
        df_method$gene_clusterid <- paste0(df_method$gene, "_", df_method$cluster_id)

        p_val <- df_method[, c("gene_clusterid", "p_val")]
        rownames(p_val) <- p_val$gene_clusterid
        p_val$gene_clusterid <- NULL
        p_vals[[method]] <- p_val

        p_val_adj <- df_method[, c("gene_clusterid", "p_val_adj")]
        rownames(p_val_adj) <- p_val_adj$gene_clusterid
        p_val_adj$gene_clusterid <- NULL
        p_adjs[[method]] <- p_val_adj

        truths <- df_method[, c("gene_clusterid", "is_DE", "category", "logFC", "sim_gene", "sim_disp", "sim_mean.A", "sim_mean.B", "seed", "test.use", "de_percent")]
        rownames(truths) <- truths$gene_clusterid
      }
      p_adjs_df <- bind_cols(p_adjs)
      colnames(p_adjs_df) <- unique(df_shortlist$method)
      p_vals_df <- bind_cols(p_vals)
      colnames(p_vals_df) <- unique(df_shortlist$method)

      cobra_data <- COBRAData(pval = p_vals_df, padj = p_adjs_df, truth = truths)
      performance <- calculate_performance(cobra_data,
                                           binary_truth = "is_DE",
                                           cont_truth = "none",
                                           aspects = c(
                                             "fdrtpr", "fdrtprcurve",
                                             "tpr", "roc"
                                           ),
                                           thrs = c(0.01, 0.05, 0.1), splv = "none"
      )

      perf_datasets[[paste(seed, test.use, as.character(de_percent), sep = "_")]] <- performance
      fdrtpr_df <- performance@fdrtpr
      fdrtpr_df$seed <- seed
      fdrtpr_df$test.use <- test.use
      fdrtpr_df$de_percent <- de_percent

      fdrtprcurve_df <- performance@fdrtprcurve
      fdrtprcurve_df$seed <- seed
      fdrtprcurve_df$test.use <- test.use
      fdrtprcurve_df$de_percent <- de_percent

      all_fdr_tpr[[paste(seed, test.use, as.character(de_percent), sep = "_")]] <- fdrtpr_df
      all_fdrtprcurve[[paste(seed, test.use, as.character(de_percent), sep = "_")]] <- fdrtprcurve_df
    }
  }
}
all_fdr_tpr_df <- bind_rows(all_fdr_tpr)
all_fdrtprcurve_df <- bind_rows(all_fdrtprcurve)

## For a seed and a test average over the TPR/FDR
all_fdr_tpr_df_melt <- all_fdr_tpr_df %>%
  select(thr, method, TPR, FDR, seed, test.use, de_percent) %>%
  mutate_at("thr", function(u) as.numeric(gsub("thr", "", u))) %>%
  group_by(test.use, de_percent, thr, method) %>%
  summarise_at(c("FDR", "TPR"), mean) %>%
  ungroup()
all_fdr_tpr_df_melt$FDR[all_fdr_tpr_df_melt$FDR==0] <- 0.001
all_fdr_tpr_df_melt <- all_fdr_tpr_df_melt[all_fdr_tpr_df_melt$method %in% c("lognorm", "scran", "SCT pr", "SCT2"),]
all_fdr_tpr_df_melt[all_fdr_tpr_df_melt$method == "SCT2", "method"] <- "SCT v2"
all_fdr_tpr_df_melt[all_fdr_tpr_df_melt$method == "SCT pr", "method"] <- "SCT v1"
saveRDS(all_fdr_tpr_df_melt, here::here("output/muscat_simulation/results/all_fdr_tpr_df_melt.rds"))



