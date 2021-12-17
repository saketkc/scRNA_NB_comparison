set.seed(42)
suppressMessages({
  library(muscat)
  library(SingleCellExperiment)
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(jsonlite)
  library(muscat)
  library(scater)
  library(sctransform)
  library(ggplot2)
  library(ggpubr)
  library(ROCR)
  library(Seurat)
  library(scran)
  library(PRROC)
  library(iCOBRA)
  library(future)
})
plan("multicore", workers = 30)


theme_set(theme_pubr())
`%notin%` <- Negate(`%in%`)

NormalizeDataSCT <- function(object, scale_factor = 'median'){
  corrected_counts <- GetAssayData(object, assay="SCT", slot="counts")
  if (scale_factor == 'median'){
    scale_factor <- median(object$nCount_RNA)
  }
  if (!is.numeric(scale_factor)){ stop("scale_factor should be 'median' for median UMI or a numeric vaule.")}
  normalized_counts <- log1p(scale_factor*corrected_counts/colSums(corrected_counts))
  object <- SetAssayData(object, slot="data", new.data = normalized_counts)
  return (object)
}


GetMarkers <- function(obj, test.use, default.assay = "SCT", slot = "data") {
  Idents(obj) <- "cluster_id"
  clusterwise_markers <- list()
  for (cluster_id in unique(obj$cluster_id)) {
    seu_subset <- subset(obj, idents = cluster_id)
    # do testing on A vs B
    print(cluster_id)
    Idents(seu_subset) <- "group_id"
    DefaultAssay(seu_subset) <- default.assay
    print(dim(seu_subset))
    markers <- FindMarkers(seu_subset,
                           ident.1 = "B",
                           slot = slot,
                           logfc.threshold = 0.,
                           test.use = test.use,
                           min.pct = 0.
    ) # , densify = TRUE)
    colnames(markers)[2] <- "avg_log2FC"
    markers$gene <- rownames(markers)
    clusterwise_markers[[cluster_id]] <- markers
  }
  clusterwise_markers_df <- bind_rows(clusterwise_markers, .id = "cluster_id")
  return(clusterwise_markers_df)
}


DoSCT <- function(object, sct_method, split_by = "sample_id") {
  object_split <- SplitObject(object, split.by = split_by)
  print(lapply(object_split, FUN = function(x) median(x$nCount_RNA)))
  min_median_umi <- min(sapply(object_split, FUN = function(x) median(x$nCount_RNA)))
  message(paste("Min median umi: ", min_median_umi))
  for (name in names(object_split)) {
    obj <- object_split[[name]]
    if (sct_method == "SCT") {
      obj <- SCTransform(obj, method = "glmGamPoi", verbose = FALSE, do.center = FALSE, do.scale = FALSE)
    } else if (sct_method == "SCT2") {
      obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE, scale_factor = min_median_umi, do.center = FALSE, do.scale = FALSE)
    } else if (sct_method == "SCT2medianscale") {
      obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE, do.center = FALSE, do.scale = FALSE)
    }
    object_split[[name]] <- obj
  }
  merged_obj <- merge(object_split[[1]], object_split[2:length(object_split)])

  return(merged_obj)
}




set.seed(42)

prep.sce <- readRDS(here::here("output/ss3_dropseq_prepsim.rds"))

roc_list <- list()
pr_list <- list()
markers_list <- list()
aupr_roc_list <- list()
seeds <- c(20141015, 20160823, 20210827)

for (seed in seeds) {

  for (test.use in c("wilcox", "MAST")) {
    for (de_percent in c(0.01, 0.05, 0.1, 0.2)) {
      t1 <- Sys.time()
      message(paste(seed, test.use, de_percent))
      set.seed(seed) # if the de percentage is same (across the two tests) the seed is the same to ensure same genes are DE across the two tests
      sim <- simData(prep.sce,
                     ng = nrow(prep.sce),
                     paired = FALSE,
                     nk = 3,
                     ns = 3,
                     nc = 200 * (3+3) * 3,
                     p_dd = c(1 - de_percent, 0, de_percent, 0., 0., 0.),
                     force = TRUE
      )
      genes_to_keep <- sample(rownames(sim), size=4000, replace=FALSE)
      sim <- sim[genes_to_keep,]
      gene_info <- sim@metadata$gene_info
      gene_info <- gene_info[gene_info$gene %in% genes_to_keep,]

      cm <- as(object = counts(sim)[, rownames(sim@colData)], Class = "dgCMatrix")
      cm <- cm[rowSums(cm)>0,]
      metadata <- sim@colData
      metadata$cluster_id <- as.character(metadata$cluster_id)
      metadata$sample_id <- as.character(metadata$sample_id)
      seu <- CreateSeuratObject(counts = cm, meta.data = as.data.frame(metadata),
                                min.cells = -1, min.features = -1)
      seu$original_sample <- stringr::str_split_fixed(seu$sample_id, pattern = "\\.", n = 2)[, 1]

      write_dir <- here::here("output/muscat_simulated/seurat_objects", seed, de_percent)
      dir.create(write_dir, showWarnings = FALSE, recursive = T)
      #saveRDS(seu, paste0(write_dir, "/", "seurat_object.rds"))
      #saveRDS(gene_info, paste0(write_dir, "/", "gene_info.rds"))

      clusters <- quickCluster(sim, min.size = 100)
      scran.sce <- computeSumFactors(sim, cluster = clusters)
      scran.sce <- logNormCounts(scran.sce)
      scran.seu <- CreateSeuratObject(counts = cm, meta.data = as.data.frame(metadata),
                                      min.cells = -1, min.features = -1) %>% NormalizeData()
      scran.seu$original_sample <- stringr::str_split_fixed(scran.seu$sample_id, pattern = "\\.", n = 2)[, 1]

      #scran.seu@assays$RNA@data <- logcounts(scran.sce)[, colnames(scran.seu@assays$RNA@counts)]
      scran.seu <- SetAssayData(scran.seu, slot = "data", new.data = logcounts(scran.sce)[rownames(scran.seu), colnames(scran.seu@assays$RNA@counts)])

      sct_obj <- DoSCT(seu, "SCT", "original_sample")
      sct2_obj <- DoSCT(seu, "SCT2", "original_sample")
      sct2medscale_obj <- DoSCT(seu, "SCT2medianscale", "original_sample")
      sct2medscale_obj <- NormalizeDataSCT(sct2medscale_obj)

      sct1_markers <- GetMarkers(sct_obj, test.use, "SCT", "data")
      sct1pr_markers <- GetMarkers(sct_obj, test.use, "SCT", "scale.data")
      sct2_markers <- GetMarkers(sct2_obj, test.use, "SCT", "data")

      sct2medscale_obj <- NormalizeDataSCT(sct2medscale_obj)
      sct2norm_markers <- GetMarkers(sct2medscale_obj, test.use, "SCT", "data")


      DefaultAssay(seu) <- "RNA"
      seu <- NormalizeData(seu)
      lognorm_markers <- GetMarkers(seu, test.use, "RNA")
      scran_markers <- GetMarkers(scran.seu, test.use, "RNA")


      sct_combined_markers <- left_join(gene_info, sct1_markers)
      sct_combined_markers$p_val_adj[is.na(sct_combined_markers$p_val_adj)] <- 1
      sct_combined_markers$is_DE <- 0
      sct_combined_markers$is_DE[sct_combined_markers$category %notin% c("ee", "ep")] <- 1
      sct_combined_markers$method <- "SCT"


      sct2_combined_markers <- left_join(gene_info, sct2_markers)
      sct2_combined_markers$p_val_adj[is.na(sct2_combined_markers$p_val_adj)] <- 1
      sct2_combined_markers$is_DE <- 0
      sct2_combined_markers$is_DE[sct2_combined_markers$category %notin% c("ee", "ep")] <- 1
      sct2_combined_markers$method <- "SCT2"

      sctpr_combined_markers <- left_join(gene_info, sct1pr_markers)
      sctpr_combined_markers$p_val_adj[is.na(sctpr_combined_markers$p_val_adj)] <- 1
      sctpr_combined_markers$is_DE <- 0
      sctpr_combined_markers$is_DE[sctpr_combined_markers$category %notin% c("ee", "ep")] <- 1
      sctpr_combined_markers$method <- "SCT pr"

      sct2norm_combined_markers <- left_join(gene_info, sct2norm_markers)
      sct2norm_combined_markers$p_val_adj[is.na(sct2norm_combined_markers$p_val_adj)] <- 1
      sct2norm_combined_markers$is_DE <- 0
      sct2norm_combined_markers$is_DE[sct2norm_combined_markers$category %notin% c("ee", "ep")] <- 1
      sct2norm_combined_markers$method <- "SCT2 normalized counts"

      scran_combined_markers <- left_join(gene_info, scran_markers)
      scran_combined_markers$p_val_adj[is.na(scran_combined_markers$p_val_adj)] <- 1
      scran_combined_markers$is_DE <- 0
      scran_combined_markers$is_DE[scran_combined_markers$category %notin% c("ee", "ep")] <- 1
      scran_combined_markers$method <- "scran"

      lognorm_combined_markers <- left_join(gene_info, lognorm_markers)
      lognorm_combined_markers$p_val_adj[is.na(lognorm_combined_markers$p_val_adj)] <- 1
      lognorm_combined_markers$is_DE <- 0
      lognorm_combined_markers$is_DE[lognorm_combined_markers$category %notin% c("ee", "ep")] <- 1
      lognorm_combined_markers$method <- "lognorm"

      markers_df <- rbind(rbind(sct_combined_markers, sct2_combined_markers),
                          rbind(sctpr_combined_markers, sct2norm_combined_markers))
      markers_df <- rbind(markers_df, rbind(scran_combined_markers, lognorm_combined_markers))
      markers_df$seed <- seed
      markers_df$test.use <- test.use
      markers_df$de_percent <- de_percent
      markers_list[[paste0(seed, test.use, as.character(de_percent))]] <- markers_df

      markers_list_df <- bind_rows(markers_list)

      t2 <- Sys.time()
      diff <- t2 - t1
      message(diff)
    }
  }
}
markers_list_df <- bind_rows(markers_list)
dir.create(here::here("output", "muscat_simulation", "results"), showWarnings = F, recursive = T)

saveRDS(markers_list_df, here::here("output/muscat_simulation/results/muscat__DE_markers_list_df.rds"))

