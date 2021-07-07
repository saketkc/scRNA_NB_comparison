library(dplyr)
library(MASS)
library(moments)
library(parallel)
library(qvalue)
library(reshape2)
library(Seurat)
library(sparseMatrixStats)
library(statmod)
set.seed(42)
mc.cores <- detectCores()
# 10

GetQvalues <- function(pvalues) {
  qvalues <- tryCatch(
    {
      qvalue(pvalues)$qvalues
    },
    error = function(x) {
      qvalue(pvalues, pi0 = 1)$qvalues
    }
  )
  return(qvalues)
}

DoShapiroTest <- function(residuals) {
  pvalues <- tryCatch(
    {
      shapiro.test(residuals)$p.value
    },
    error = function(x) {
      rep(1, length(residuals))
    }
  )
  return(pvalues)
}

GeneSummary <- function(cm) {

  # remove genes and cells with zero counts
  cm <- cm[rowSums(cm) > 0, colSums(cm) > 0]

  total_umi_per_gene <- rowSums(cm)
  expressed_cells_per_gene <- rowSums(cm > 0)
  n_cells <- dim(cm)[2]
  nonexpressed_cells_per_gene <- n_cells - expressed_cells_per_gene

  median_umi_per_gene <- median(total_umi_per_gene)

  avg_umi_per_gene <- total_umi_per_gene / n_cells
  avg_umi_per_gene_expressedcells <- total_umi_per_gene / expressed_cells_per_gene

  gene_amean <- rowMeans(cm)
  gene_var <- rowVars(cm)
  gene_gmean <- sctransform:::row_gmean(cm)

  gene_attr <- data.frame(
    total_umi = total_umi_per_gene,
    n_expressed_cells = expressed_cells_per_gene,
    n_nonexpressed_cells = nonexpressed_cells_per_gene,
    prop_expressed_cells = expressed_cells_per_gene / n_cells,
    prop_nonexpressed_cells = nonexpressed_cells_per_gene / n_cells,
    avg_umi = avg_umi_per_gene,
    avg_umi_expressedcells = avg_umi_per_gene_expressedcells,
    gene_amean = gene_amean,
    gene_gmean = gene_gmean,
    gene_variance = gene_var
  )

  return(gene_attr)
}

CellSummary <- function(cm) {
  total_umi_per_cell <- colSums(cm)
  expressed_features_per_cell <- colSums(x = cm > 0)
  n_features <- dim(cm)[1]

  nonexpressed_features_per_cell <- n_features - expressed_features_per_cell
  median_umi_per_cell <- median(total_umi_per_cell)
  avg_umi_per_cell <- total_umi_per_cell / n_features
  avg_umi_per_cell_expressedgenes <- total_umi_per_cell / expressed_features_per_cell
  cell_amean <- colMeans2(cm)
  cell_variance <- colVars(cm)
  cell_attr <- data.frame(
    total_umi = total_umi_per_cell, n_expressed_genes = expressed_features_per_cell,
    n_nonexpressed_cells = nonexpressed_features_per_cell,
    prop_expressed_genes = expressed_features_per_cell / n_features,
    prop_nonexpressed_genes = nonexpressed_features_per_cell / n_features,
    avg_umi = avg_umi_per_cell, avg_umi_expressedgenes = avg_umi_per_cell_expressedgenes, cell_amean = cell_amean,
    cell_variance = cell_variance
  )

  return(cell_attr)
}

MeanVarFit <- function(counts) {
  means <- sparseMatrixStats::rowMeans2(counts)
  variance <- sparseMatrixStats::rowVars(counts)
  df <- data.frame(gene = rownames(counts), mean = means, variance = variance)
  return(df)
}

dds <- function(libsizes, cells, ncells) {
  libsize_density <- density(x = libsizes, bw = "nrd", adjust = 1)
  sampling_prob <- 1 / (approx(x = libsize_density$x, y = libsize_density$y, xout = libsizes)$y + .Machine$double.eps)
  sampled_cells <- sample(x = cells, size = ncells, prob = sampling_prob)
  return(sampled_cells)
}

GetChiSqPvalue <- function(ssr_residuals, df) {
  pchisq(ssr_residuals, df = df, lower.tail = FALSE)
}

ssr <- function(r) {
  sum(r^2)
}

RunGLM <- function(gene_umi, cell_umi, type = "poisson") {
  df <- data.frame(gene_umi = gene_umi, cell_umi = cell_umi)

  if (type == "poisson") {
    fit <- glm(gene_umi ~ offset(log(cell_umi)), data = df, family = poisson(link = "log"))
    fit_type <- "poisson"
  } else if (type == "negbin") {
    # this might fail sometimes
    fit_type <- "negbin"
    fit <- tryCatch(
      {
        glm.nb(gene_umi ~ offset(log(cell_umi)), data = df)
      },
      error = function(x) {
        fit_type <- "negbinfailed-poisson"
        glm(gene_umi ~ offset(log(cell_umi)), data = df, family = poisson(link = "log"))
      }
    )
  } else if (type == "quasipoisson") {
    fit_type <- "quasipoisson"
    fit <- glm(gene_umi ~ offset(log(cell_umi)), data = df, family = quasipoisson(link = "log"))
  }

  residuals_deviance <- residuals(fit, type = "deviance")
  residuals_pearson <- residuals(fit, type = "pearson")
  residuals_quantile <- statmod::qresid(fit)

  resid_mean_deviance <- mean(residuals_deviance)
  resid_mean_pearson <- mean(residuals_pearson)
  resid_mean_quantile <- mean(residuals_quantile)

  resid_sd_deviance <- sd(residuals_deviance)
  resid_sd_pearson <- sd(residuals_pearson)
  resid_sd_quantile <- sd(residuals_quantile)

  resid_kurtosis_deviance <- kurtosis(residuals_deviance)
  resid_kurtosis_pearson <- kurtosis(residuals_pearson)
  resid_kurtosis_quantile <- kurtosis(residuals_quantile)

  resid_skew_deviance <- skewness(residuals_deviance)
  resid_skew_pearson <- skewness(residuals_pearson)
  resid_skew_quantile <- skewness(residuals_quantile)


  resid_swpval_deviance <- DoShapiroTest(residuals_deviance)
  resid_swpval_pearson <- DoShapiroTest(residuals_pearson)
  resid_swpval_quantile <- DoShapiroTest(residuals_quantile)




  phi.est <- summary(fit)$dispersion # 1
  h <- hatvalues(fit)

  residstd_deviance <- rstandard(fit)
  residstd_pearson <- residuals_pearson / sqrt(phi.est * (1 - h))
  residstd_quantile <- residuals_quantile / sqrt(1 - h)

  ssr_deviance <- ssr(residuals_deviance)
  ssr_pearson <- ssr(residuals_pearson)
  ssr_quantile <- ssr(residuals_quantile)

  ssrstd_deviance <- ssr(residstd_deviance)
  ssrstd_pearson <- ssr(residstd_pearson)
  ssrstd_quantile <- ssr(residstd_quantile)



  pval_ssr_dev <- GetChiSqPvalue(ssr_deviance, fit$df.residual)
  pval_ssr_pear <- GetChiSqPvalue(ssr_pearson, fit$df.residual)
  pval_ssr_quan <- GetChiSqPvalue(ssr_quantile, fit$df.residual)


  pval_ssrstd_dev <- GetChiSqPvalue(ssrstd_deviance, fit$df.residual)
  pval_ssrstd_pear <- GetChiSqPvalue(ssrstd_pearson, fit$df.residual)
  pval_ssrstd_quan <- GetChiSqPvalue(ssrstd_quantile, fit$df.residual)


  pvals <- data.frame(
    dof = fit$df.residual,
    fit_type = fit_type,
    ssr_deviance = ssr_deviance,
    ssr_deviance_normalized = ssr_deviance / fit$df.residual,
    pval_ssr_deviance = pval_ssr_dev,
    qval_ssr_deviance = pval_ssr_dev,
    ssr_pearson = ssr_pearson,
    ssr_pearson_normalized = ssr_pearson / fit$df.residual,
    pval_ssr_pearson = pval_ssr_pear,
    qval_ssr_pearson = pval_ssr_pear,
    ssr_quantile = ssr_quantile,
    ssr_quantile_normalized = ssr_quantile / fit$df.residual,
    pval_ssr_quantile = pval_ssr_quan,
    qval_ssr_quantile = pval_ssr_quan,
    ssrstd_deviance = ssrstd_deviance,
    ssrstd_deviance_normalized = ssrstd_deviance / fit$df.residual,
    pval_ssrstd_deviance = pval_ssrstd_dev,
    qval_ssrstd_deviance = pval_ssrstd_dev,
    ssrstd_pearson = ssrstd_pearson,
    ssrstd_pearson_normalized = ssrstd_pearson / fit$df.residual,
    pval_ssrstd_pearson = pval_ssrstd_pear,
    qval_ssrstd_pearson = pval_ssrstd_pear,
    ssrstd_quantile = ssrstd_quantile,
    ssrstd_quantile_normalized = ssrstd_quantile / fit$df.residual,
    pval_ssrstd_quantile = pval_ssrstd_quan,
    qval_ssrstd_quantile = pval_ssrstd_quan,
    resid_mean_deviance = resid_mean_deviance,
    resid_mean_pearson = resid_mean_pearson,
    resid_mean_quantile = resid_mean_quantile,
    resid_sd_deviance = resid_sd_deviance,
    resid_sd_pearson = resid_sd_pearson,
    resid_sd_quantile = resid_sd_quantile,
    resid_kurtosis_deviance = resid_kurtosis_deviance,
    resid_kurtosis_pearson = resid_kurtosis_pearson,
    resid_kurtosis_quantile = resid_kurtosis_quantile,
    resid_skew_deviance = resid_skew_deviance,
    resid_skew_pearson = resid_skew_pearson,
    resid_skew_quantile = resid_skew_quantile,
    resid_swpval_deviance = resid_swpval_deviance,
    resid_swpval_pearson = resid_swpval_deviance,
    resid_swpval_quantile = resid_swpval_quantile
  )

  return(pvals)
}

GetFit <- function(cm, ncells = dim(cm)[2], type = "poisson") {
  cm <- cm[rowSums(cm) > 0, colSums(cm) > 0]
  ncells <- min(ncells, dim(cm)[2])
  total_umi <- colSums(cm)
  if (ncells == dim(cm)[2]) {
    # No sampling
    sampled_cells <- colnames(cm)
  } else {
    sampled_cells <- dds(total_umi, colnames(cm), ncells)
  }

  cm_sampled <- cm[, sampled_cells]
  total_umi_sampled <- colSums(cm_sampled)

  genes <- rownames(cm_sampled)
  #pvals_list <- mclapply(genes, FUN = function(gene) {
  #  gene_umi <- as.vector(cm_sampled[gene, ])

   # pvals <- RunGLM(gene_umi = gene_umi, cell_umi = total_umi_sampled, type = type)
   # pvals
  #}, mc.cores = mc.cores)
  pvals_list <- lapply(genes, FUN = function(gene) {
    gene_umi <- as.vector(cm_sampled[gene, ])

    pvals <- RunGLM(gene_umi = gene_umi, cell_umi = total_umi_sampled, type = type)
    pvals
  })
  names(pvals_list) <- genes

  pvalsfit <- bind_rows(pvals_list, .id = "gene")
  meanvar <- MeanVarFit(cm)
  meanvar$total_gene_umi <- rowSums(cm)

  meanvar$total_cell_umi <- sum(total_umi_sampled)
  meanvar$median_cell_umi <- median(total_umi_sampled)
  meanvar$mean_cell_umi <- mean(total_umi_sampled)

  meanvarfit <- inner_join(meanvar, pvalsfit, by = "gene")

  qval_cols <- c("qval_ssr_deviance", "qval_ssr_pearson", "qval_ssr_quantile", "qval_ssrstd_deviance", "qval_ssrstd_pearson", "qval_ssrstd_quantile")
  for (qval_col in qval_cols) {
    pval_col <- gsub("qval", "pval", qval_col)
    meanvarfit[, qval_col] <- GetQvalues(meanvarfit[, pval_col])
  }


  return(list(meanvarfit = meanvarfit, cell_attr = CellSummary(cm_sampled), gene_attr = GeneSummary(cm_sampled)))
}



args <- commandArgs(trailingOnly = TRUE)
# argument 1 => location to RDS
# argument 2 => Number of cells to use
# argument 3 => prefix to output

rds <- args[1]
output_path <- args[4]
glmtype <- args[3]
ncells <- as.numeric(args[2])


object <- readRDS(rds)
seed <- 42

cm <- GetAssayData(object, assay = "RNA", slot = "counts")
cm <- cm[rowSums(cm) > 0, colSums(cm) > 0]
if (is.na(ncells)) {
  ncells <- dim(cm)[2]
}
obj <- GetFit(cm, ncells, glmtype)
fits <- obj$meanvarfit
cell_attr <- obj$cell_attr
gene_attr <- obj$gene_attr
write.csv(fits, output_path)
write.csv(cell_attr, gsub(".csv", "_cell_attr.csv", output_path))
write.csv(gene_attr, gsub(".csv", "_gene_attr.csv", output_path))
