# edit this line
library(Seurat)
library(sctransform) 
library(tidyverse)
library(future)
options(future.globals.maxSize = 9991289600)

args <- commandArgs(trailingOnly = TRUE)
# argument 1 => location to RDS
# argument 2 => UMIs to keep
# argument 3 => Mode of SCT (glmGamPoi/glmGamPoi2/offset-10/offset-100/offset-Inf/offset_shared_theta-100/offset_shared_theta-500/offset_shared_theta-1000)
# argument 4 => prefix to output

rds <- args[1]
seed <- as.numeric(args[2])
downsample_cells <- as.numeric(args[3])
output_dir <- args[4]

set.seed(seed)
obj <- readRDS(rds)
obj <- subset(obj, cells=sample(Cells(obj), downsample_cells))
cm <- GetAssayData(object = obj, assay = "RNA", slot = "counts")
## sct 

start_time <- Sys.time()
vst.out.sct <- sctransform::vst(umi = cm, method="glmGamPoi", n_genes=2000, n_cells = NULL, verbosity =0)
end_time <- Sys.time()

sct_difftime <-  as.numeric( as.numeric(end_time-start_time, units="secs"))

start_time <- Sys.time()
vst.out.sct2 <- sctransform::vst(
  umi = cm, 
  method = "glmGamPoi_offset",
  exclude_poisson = TRUE,
  n_genes = 2000, n_cells = 5000,
  verbosity = 0
)
end_time <- Sys.time()
sct2_difftime <-  as.numeric( as.numeric(end_time-start_time, units="secs"))

dir.create(file.path(output_dir), showWarnings = FALSE)
output_prefix <- file.path(output_dir)

#saveRDS(cm, file.path(output_prefix, "count_matrix.rds"))
#saveRDS(vst.out.sct, file.path(output_prefix, "vst_out_sct.rds"))
#saveRDS(vst.out.sct2, file.path(output_prefix, "vst_out_sct2.rds"))

times <- data.frame(ncells=downsample_cells, time_sct = sct_difftime, time_sct2=sct2_difftime)
write.csv(times, file.path(output_prefix, "times.csv"))


model_pars <- as.data.frame(vst.out.sct$model_pars)
model_pars$gene <- rownames(model_pars)
write.csv(model_pars, file.path(output_prefix, "model_pars_sct.csv"))

model_pars_fit <- as.data.frame(vst.out.sct$model_pars_fit)
model_pars_fit$gene <- rownames(model_pars_fit)
write.csv(model_pars_fit, file.path(output_prefix, "model_fit_sct.csv"))

gene_attr <- as.data.frame(vst.out.sct$gene_attr)
gene_attr$gene <- rownames(gene_attr)
write.csv(gene_attr, file.path(output_prefix, "gene_attr_sct.csv"))


cell_attr <- as.data.frame(vst.out.sct$cell_attr)
cell_attr$gene <- rownames(cell_attr)
write.csv(cell_attr, file.path(output_prefix, "cell_attr_sct.csv"))




model_pars <- as.data.frame(vst.out.sct2$model_pars)
model_pars$gene <- rownames(model_pars)
write.csv(model_pars, file.path(output_prefix, "model_pars_sct2.csv"))

model_pars_fit <- as.data.frame(vst.out.sct2$model_pars_fit)
model_pars_fit$gene <- rownames(model_pars_fit)
write.csv(model_pars_fit, file.path(output_prefix, "model_fit_sct2.csv"))

gene_attr <- as.data.frame(vst.out.sct2$gene_attr)
gene_attr$gene <- rownames(gene_attr)
write.csv(gene_attr, file.path(output_prefix, "gene_attr_sct2.csv"))


cell_attr <- as.data.frame(vst.out.sct2$cell_attr)
cell_attr$gene <- rownames(cell_attr)
write.csv(cell_attr, file.path(output_prefix, "cell_attr_sct2.csv"))
