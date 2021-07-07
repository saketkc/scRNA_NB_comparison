library(Seurat)
library(sctransform) 
library(future)
options(future.globals.maxSize = 8991289600)

args <- commandArgs(trailingOnly = TRUE)
# argument 1 => location to RDS
# argument 2 => number of cells
# argument 3 => prefix to output

rds <- args[1]
ncells <- as.numeric(args[2])
output_dir <- args[3]
dir.create(file.path(output_dir), showWarnings = FALSE)
output_prefix <- file.path(output_dir)

set.seed(42)
obj <- readRDS(rds)
cm <- GetAssayData(object = obj, assay = "RNA", slot = "counts")
# Use all cells
if (is.na(ncells)){
  ncells <- dim(cm)[2]
}
## sct 

start_time <- Sys.time()
vst.out.sct <- SCTransform(obj, method="glmGamPoi", ncells = ncells, do.correct.umi=FALSE, conserve.memory=TRUE, verbosity =2)
end_time <- Sys.time()
sct_difftime <-  as.numeric( as.numeric(end_time-start_time, units="secs"))
gene_attr <- SCTResults(vst.out.sct, slot = "feature.attributes", assay = "SCT")
gene_attr$gene <- rownames(gene_attr)
write.csv(gene_attr, file.path(output_prefix, "gene_attr_sct.csv"))

vst.out.sct <- NULL 
gc()

## sct2
start_time <- Sys.time()
vst.out.sct2 <- SCTransform(obj, method = "glmGamPoi_offset", exclude_poisson=TRUE, 
                             ncells = ncells,  do.correct.umi=FALSE, conserve.memory = TRUE, verbosity = 2)
end_time <- Sys.time()
sct2_difftime <-  as.numeric( as.numeric(end_time-start_time, units="secs"))
gene_attr <- SCTResults(vst.out.sct2, slot = "feature.attributes", assay = "SCT")
gene_attr$gene <- rownames(gene_attr)
write.csv(gene_attr, file.path(output_prefix, "gene_attr_sct2.csv"))
times <- data.frame(ncells=ncells, time_sct = sct_difftime, time_sct2=sct2_difftime)
write.csv(times, file.path(output_prefix, "times.csv"))

vst.out.sct2 <- NULL 
gc()
