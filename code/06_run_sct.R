
library(future)
library(Seurat)
library(scattermore)
library(ggplot2)
library(patchwork)
library(magrittr)
options(future.globals.maxSize = 8991289600)
args <- commandArgs(trailingOnly = TRUE)

rds <- args[1]
sct_method <- args[2]
output_prefix <- args[3]

object <- readRDS(rds)
n_genes <- dim(object)[1]
ncells  <- dim(object)[2]
if (ncells> 50000){
  #only use 5k cells for large datasets
  ncells <- 5000
}
seed <- 42
if (sct_method == "vst2"){
  times <- list(start_time = Sys.time())
  object <- SCTransform(object, 
                        ncells = ncells , 
                        n_genes = n_genes,
                        method = "glmGamPoi_offset",
                        exclude_poisson = TRUE,
                        conserve.memory = TRUE,
                        seed.use = seed,
                        verbosity = 0)
  times$stop_time = Sys.time()
} else if (grepl("offset", sct_method)) {
  method <- strsplit(sct_method, "-")[[1]][[1]]
  theta_given <- as.numeric(strsplit(sct_method, "-")[[1]][[2]])
  times <- list(start_time = Sys.time())
  object <- SCTransform(object, 
                        ncells = ncells,
                        n_genes = n_genes,
                        method = method,
                        theta_given = theta_given,
                        exclude_poisson = TRUE,
                        use_geometric_mean = FALSE,
                        conserve.memory = TRUE,
                        seed.use = seed,
                        verbosity = 0)
  times$stop_time = Sys.time()
} else if (sct_method == "fixed-intercept"){
  times <- list(start_time = Sys.time())
  object <- SCTransform(object, 
                        ncells = ncells , 
                        n_genes = n_genes,
                        method = "glmGamPoi", 
                        exclude_poisson = TRUE, 
                        fix_intercept = TRUE, 
                        use_geometric_mean = FALSE,
                        conserve.memory = TRUE,
                        seed.use = seed,
                        verbosity = 0)
  times$stop_time = Sys.time()
} else if (sct_method == "fixed-slope"){
  times <- list(start_time = Sys.time())
  object <- SCTransform(object,
                        ncells = ncells, 
                        n_genes = n_genes,
                        method = "glmGamPoi", 
                        exclude_poisson = TRUE, 
                        fix_slope = TRUE, 
                        use_geometric_mean = FALSE,
                        conserve.memory = TRUE,
                        seed.use = seed,
                        verbosity = 0)
  times$stop_time = Sys.time()
} else if (sct_method == "fixed-interceptandslope"){
  times <- list(start_time = Sys.time())
  object <- SCTransform(object, 
                        ncells = ncells, 
                        n_genes = n_genes,
                        method = "glmGamPoi", 
                        exclude_poisson = TRUE, 
                        fix_slope = TRUE, 
                        fix_intercept = TRUE,
                        use_geometric_mean = FALSE,
                        conserve.memory = TRUE,
                        seed.use = seed,
                        verbosity = 0)
  times$stop_time = Sys.time()
} else {
  times <- list(start_time = Sys.time())
  object <- SCTransform(object,
                        method = sct_method,
                        ncells = ncells,
                        n_genes = n_genes,
                        seed.use = seed,
                        conserve.memory = TRUE,
                        verbosity = 0
                        )
  times$stop_time = Sys.time()
}

#saveRDS(object,  file.path(output_prefix, "seurat_sct_object.rds"))




residualVarPlot <- function(gene_var, xaxis="gmean", max_resvar=100, ntop = 30, annotate = F) {
  gene_var$gene <- rownames(gene_var)
  topn <- subset(gene_var, rank(-gene_var[, "residual_variance"]) <= ntop)$gene
  gene_var[gene_var$residual_variance>max_resvar, "residual_variance"] <- max_resvar
  p <- ggplot(gene_var, aes_string(xaxis, "residual_variance")) +
    geom_scattermore(pointsize = 1.1, shape = 16, alpha = 0.5, color = "#43a2ca") +
    geom_scattermore(data = subset(gene_var, gene %in% topn), pointsize = 1.1, shape = 16, alpha = 1.0, color = "deeppink") +
    geom_hline(yintercept = 1, color = "#4daf4a", size = 0.9, linetype = "dashed") +
    geom_smooth(method = "loess", span = 0.1, size = 0.9, formula = "y ~ x", color = "#e41a1c") +
    scale_y_continuous(trans = "sqrt", breaks = c(0, 1, 10, 25, 50, 100, 150), limits = c(0, max_resvar+1)) +
    scale_x_continuous(trans = "log10", breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = MASS::rational) +
    # facet_wrap(~ model, ncol=3, scales = 'free_y') +
    xlab("Gene mean") +
    ylab("Residual variance")
  if (annotate) {
    p <- p + geom_text_repel(
      data = subset(gene_var, gene %in% topn), aes(label = gene), color = "gray25",
      size = 1.8,
      nudge_y = 230 - subset(gene_var, gene %in% topn)[, col],
      direction = "x",
      angle = 90,
      vjust = 0.5,
      hjust = 0.5,
      segment.size = 0.2,
      segment.alpha = 0.2
    )
  }

  return(p)
}


maxdims <- min(30, dim(object)[2]-1)
object <- RunPCA(object = object, npcs=min(maxdims, 50), verbose = FALSE)
dir.create(output_prefix, showWarnings = FALSE)
saveRDS(object,  file.path(output_prefix, "seurat_sct_object.rds"))

gene_attr <- SCTResults(object, slot = "feature.attributes", assay="SCT")
output_geneattr <- file.path(output_prefix, "gene_attr.csv")
write.csv(gene_attr, 
          output_geneattr)
