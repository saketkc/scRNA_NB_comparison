library(Seurat)
library(metacell)

writeMetaCells <- function(metacell_name, outputfile) {
  mc <- scdb_mc(metacell_name)
  lfp <- log2(mc@mc_fp)
  df <- data.frame(
    barcode = names(mc@mc),
    metacell_id = mc@mc
  )
  write.table(df,
    file = outputfile,
    quote = FALSE, sep = "\t", row.names = FALSE
  )
  return (df)
}

ProcessMC <- function(metacell_name, dataset_name) {
  # selecting markers
  mcell_gset_from_mc_markers(gset_id = metacell_name, mc_id = metacell_name)
  # plot heatmap
  mcell_mc_plot_marks(mc_id = metacell_name, gset_id = metacell_name, mat_id = dataset_name)
  # 2dproj
  mcell_mc2d_force_knn(metacell_name, metacell_name, dataset_name)
  #tgconfig::set_param("mcell_mc2d_height", 1000, "metacell")
  #tgconfig::set_param("mcell_mc2d_width", 1000, "metacell")
  #mcell_mc2d_plot(mc2d_id = metacell_name)

}

#' @import metacell
#' @export

RunMetaCell <- function(obj.seurat, dataset_name, output.prefix,
                        color_key = NULL,
                        min_UMI = 800, dsamp = T,
                        T_tot = 100, T_top3 = 2, T_lfc = 3,
                        Knn = 100, K_coclust = 30,
                        min_mc_size = 20, T_vm = 0.1,
                        p_resamp = 0.75, n_resamp = 500,
                        alpha = 2,
                        filter_bad_genes = TRUE,
                        extra_bad_genes = NULL) {
  outpath <- file.path(output.prefix, dataset_name)
  mc_name <- sprintf("%s_bs_%s", dataset_name, n_resamp)
  mcf_name <- sprintf("%s_bs_%s_filtered", dataset_name, n_resamp)
  coc_name <- sprintf("%s_coc_%s", dataset_name, n_resamp)
  if (!dir.exists(outpath)) dir.create(outpath, showWarnings = F, recursive = T)

  scdb_init(outpath, force_reinit = T)
  scfigs_init(outpath)

  if (dim(obj.seurat)[2]<1000){
  Knn <- 20
  K_coclust <- 20
  }
  sce <- as.SingleCellExperiment(obj.seurat)
  
  mat <- scm_import_sce_to_mat(sce)
  scdb_add_mat(dataset_name, mat)

  mcell_plot_umis_per_cell(dataset_name)
  # ignore cells
  mcell_mat_ignore_small_cells(dataset_name, dataset_name, min_UMI)

  # compute stats for each gene
  mcell_add_gene_stat(
    gstat_id = dataset_name,
    mat_id = dataset_name,
    force = T
  )

  # create a new object of type gset adding all genes whose scaled variance
  # exceeds a given threshold
  mcell_gset_filter_varmean(
    gstat_id = dataset_name,
    gset_id = dataset_name,
    T_vm = T_vm,
    force_new = T
  )

  # restruct gene set to genes with at least 100 UMIs and have atleast 3 cells for more than 2UMIs
  mcell_gset_filter_cov(gset_id = dataset_name, gstat_id = dataset_name, T_tot = T_tot, T_top3 = T_top3)

  if (isTRUE(filter_bad_genes)) {
    # all_genes <- c(rownames(mat@mat), rownames(mat@ignore_gmat))
    gset <- scdb_gset(dataset_name)
    all_genes <- names(gset@gene_set)

    ig_genes <- c(
      grep("^IGJ", all_genes, v = T),
      grep("^IGH", all_genes, v = T),
      grep("^IGK", all_genes, v = T),
      grep("^IGL", all_genes, v = T)
    )
    bad_genes <- unique(c(
      grep("^MT-", all_genes, ignore.case = T, value = T),
      # grep("^MTMR", all_genes, v = T),
      # grep("^MTND", all_genes, v = T),
      # "NEAT1", "TMSB4X", "TMSB10",
      grep("^RPL", all_genes, ignore.case = T, value = T),
      grep("^GM", all_genes, ignore.case = T, value = T), # mouse
      grep("RPS", all_genes, ignore.case = T, value = T)
    ))
    # ignore genes
    if (!is.null(extra_bad_genes)) {
      bad_genes <- c(bad_genes, extra_bad_genes)
    }
    #gset_final <- gset_new_restrict_nms(gset = gset, bad_genes, inverse = T, "filtered genelist")
    #scdb_add_gset(dataset_name, gset_final)
    mcell_mat_ignore_genes(new_mat_id = dataset_name, mat_id = dataset_name, bad_genes, reverse = F)
  }



  # plot all genes vs selected gene
  mcell_plot_gstats(gstat_id = dataset_name, gset_id = dataset_name)


  # create similarity graph
  mcell_add_cgraph_from_mat_bknn(
    mat_id = dataset_name,
    gset_id = dataset_name,
    graph_id = dataset_name,
    K = Knn,
    dsamp = dsamp
  )

  # resample graph to get consensus
  mcell_coclust_from_graph_resamp(
    coc_id = coc_name,
    graph_id = dataset_name,
    min_mc_size = min_mc_size,
    n_resamp = n_resamp,
    p_resamp = p_resamp
  )


  # create a new similarity graph using the co-clustering stats
  # before calling the final metacells
  mcell_mc_from_coclust_balanced(
    mc_id = mc_name,
    coc_id = coc_name,
    mat_id = dataset_name,
    K = K_coclust,
    min_mc_size = min_mc_size, alpha = alpha
  )

  # make sure metacells are homogeneous by removing outliers
  # plot outliers
  mcell_plot_outlier_heatmap(mc_id = mc_name, mat_id = dataset_name, T_lfc = T_lfc)
  # filter
  mcell_mc_split_filt(
    new_mc_id = mcf_name,
    mc_id = mc_name,
    mat_id = dataset_name,
    T_lfc = T_lfc, plot_mats = F
  )


  if (!is.null(color_key)) {
    marks_colors <- read.table(color_key, sep = "\t", header = T, stringsAsFactors = F)
    mc_colorize(mc_name, marker_colors = marks_colors)
    mc_colorize(mcf_name, marker_colors = marks_colors)
  } else {
    mc_colorize_default(mc_name)
    mc_colorize_default(mcf_name)
  }


  ProcessMC(mc_name, dataset_name)
  metadata <- writeMetaCells(mc_name, file.path(output.prefix, "metacell_clusters_raw.tsv"))
  rownames(metadata) <- metadata$barcode
  colnames(metadata)[2] <- "metacell_clusters_raw"
  obj.seurat <- AddMetaData(obj.seurat, metadata=metadata)

  ProcessMC(mcf_name, dataset_name)
  metadata <- writeMetaCells(mcf_name, file.path(output.prefix, "metacell_clusters_filtered.tsv"))
  rownames(metadata) <- metadata$barcode
  colnames(metadata)[2] <- "metacell_clusters_filtered"
  obj.seurat <- AddMetaData(obj.seurat, metadata=metadata)
  saveRDS(obj.seurat, file.path(output.prefix, "seurat_metacell_object.rds"))
  
  obj.seurat_subset <- subset(obj.seurat, cells=rownames(metadata))
  saveRDS(obj.seurat_subset, file.path(output.prefix, "seurat_metacell_filtered_subset.rds"))
  return (obj.seurat)

}


args <- commandArgs(trailingOnly = TRUE)

rds <- args[1]
Knn <- as.numeric(args[2])
output_prefix <- args[3]
object <- readRDS(rds)
seed <- 42

seu <- RunMetaCell(object, gsub(".rds", "", basename(rds)), output_prefix, Knn=Knn)

seu.raw_agg <- AggregateExpression(seu, assays=c("RNA"), slot="counts", group.by="metacell_clusters_raw", return.seurat=TRUE)
seu.raw_agg <- RenameCells(seu.raw_agg, add.cell.id="metacell_clusters_raw")
saveRDS(seu.raw_agg, file.path(output_prefix, "seurat_metacell_aggregatedcounts_raw.rds"))

seu.filtered_agg <- AggregateExpression(seu, assays=c("RNA"), slot="counts", group.by="metacell_clusters_filtered", return.seurat=TRUE)
seu.filtered_agg <- RenameCells(seu.filtered_agg, add.cell.id="metacell_clusters_filtered")
saveRDS(seu.filtered_agg, file.path(output_prefix, "seurat_metacell_aggregatedcounts_filtered.rds"))
