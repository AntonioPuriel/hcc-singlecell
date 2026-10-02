# Progenitor-like (stem-like candidate) tumour states --------------------------------
#
# Revised after the second run. The first definition intersected the top quartile
# of a CSC marker score with the top quartile of a within-sample gene-count rank.
# The two readouts were uncorrelated in every patient (Spearman -0.17 to 0.28), so
# the intersection held ~6% of cells, which is what chance gives (0.25 x 0.25);
# the gene count also depends on sequencing depth and cell size.
#
# Current definition, within CNV-confirmed malignant cells:
#   prog_score:  module score of the hepatic progenitor / biliary programme
#                (EPCAM, KRT19, KRT7, SOX9, CD24, PROM1); in the malignant-cell
#                clustering these genes are co-expressed in the same clusters
#   progenitor-like: top quartile of prog_score within each patient, so every
#                patient contributes both groups to the paired pseudobulk DE
# Independent check, not used in the definition:
#   CytoTRACE2 (Kang et al. 2024), a trained model of developmental potential
#   whose input is the whole transcriptome. Per patient we report the Spearman
#   correlation with prog_score and the difference in median CytoTRACE2 score
#   between progenitor-like and other malignant cells.
# The state is called "progenitor-like" because it is defined by a lineage
# programme; it is described as stem-like only where CytoTRACE2 agrees.

progenitor_markers <- c("EPCAM", "KRT19", "KRT7", "SOX9", "CD24", "PROM1")

malignant_cells <- function(seu) {
  mal <- subset(seu, cells = colnames(seu)[seu$cnv_call %in% "malignant"])
  mal <- Seurat::DietSeurat(mal, layers = "counts")
  Seurat::NormalizeData(mal, verbose = FALSE)
}

#' CytoTRACE2 on the raw counts of the malignant cells. Few cores: the package
#' forks per model, and each fork holds a dense copy of its batch.
run_cytotrace2 <- function(mal, seed = 14) {
  counts <- SeuratObject::LayerData(mal, assay = "RNA", layer = "counts")
  counts <- counts[Matrix::rowSums(counts) > 0, ]
  suppressPackageStartupMessages(library(CytoTRACE2))   # internals rely on its Depends being attached
  res <- CytoTRACE2::cytotrace2(as.matrix(counts), species = "human", is_seurat = FALSE,
                                batch_size = 5000, smooth_batch_size = 1000,
                                ncores = min(4L, n_workers()), seed = seed)
  data.frame(cell = rownames(res), cytotrace2 = res$CytoTRACE2_Score,
             cytotrace2_potency = as.character(res$CytoTRACE2_Potency))
}

score_progenitor <- function(mal, cyto, markers = progenitor_markers, q = 0.75, seed = 1234) {
  genes <- intersect(markers, rownames(mal))
  mal <- Seurat::AddModuleScore(mal, features = list(genes), name = "prog", seed = seed,
                                ctrl = min(100L, floor(nrow(mal) / 50)))
  mal$prog_score <- mal$prog1; mal$prog1 <- NULL
  i <- match(colnames(mal), cyto$cell)
  mal$cytotrace2 <- cyto$cytotrace2[i]
  mal$cytotrace2_potency <- cyto$cytotrace2_potency[i]
  md <- mal[[]]
  mal$gene_count_rank <- stats::ave(md$nFeature_RNA, md$sample, FUN = function(x) rank(x) / length(x))
  mal$progenitor_like <- as.logical(stats::ave(md$prog_score, md$patient,
                                               FUN = function(v) v >= stats::quantile(v, q)))
  mal$prog_group <- ifelse(mal$progenitor_like, "progenitor", "other")
  mal
}

cluster_malignant <- function(mal, dims = 1:20, resolution = 0.5, seed = 1234, batch_correct = TRUE) {
  mal <- normalise_and_pca(mal, n_hvg = 2000, n_pcs = max(dims), seed = seed)
  red <- "pca"
  if (batch_correct) {
    set.seed(seed)
    mal <- harmony::RunHarmony(mal, group.by.vars = "patient", reduction.use = "pca",
                               dims.use = dims, project.dim = FALSE, verbose = FALSE)
    red <- "harmony"
  }
  mal <- Seurat::RunUMAP(mal, reduction = red, dims = dims, seed.use = seed, verbose = FALSE)
  mal <- Seurat::FindNeighbors(mal, reduction = red, dims = dims, verbose = FALSE)
  mal$malignant_state <- leiden_clusters(mal, resolution, seed)
  mal
}

progenitor_summary <- function(mal) {
  md <- mal[[]]
  do.call(rbind, lapply(split(md, md$patient), function(d) {
    p <- d$progenitor_like
    w <- if (sum(p) >= 5 && sum(!p) >= 5)
      stats::wilcox.test(d$cytotrace2[p], d$cytotrace2[!p])$p.value else NA
    data.frame(
      patient = d$patient[1], sites = paste(sort(unique(d$site)), collapse = ","),
      n_malignant = nrow(d), n_progenitor_like = sum(p),
      spearman_prog_vs_cytotrace2 = round(stats::cor(d$prog_score, d$cytotrace2, method = "spearman"), 2),
      median_cytotrace2_progenitor = round(stats::median(d$cytotrace2[p]), 3),
      median_cytotrace2_other = round(stats::median(d$cytotrace2[!p]), 3),
      wilcox_p = signif(w, 2),
      spearman_prog_vs_gene_count = round(stats::cor(d$prog_score, d$gene_count_rank, method = "spearman"), 2))
  }))
}

plot_progenitor <- function(mal, markers = progenitor_markers) {
  p1 <- Seurat::DimPlot(mal, group.by = "patient", raster = TRUE, pt.size = 2) + ggplot2::ggtitle("Patient")
  p2 <- Seurat::FeaturePlot(mal, "prog_score", raster = TRUE, pt.size = 2, order = TRUE) +
    ggplot2::scale_colour_viridis_c() + ggplot2::ggtitle("Progenitor programme score")
  p3 <- Seurat::FeaturePlot(mal, "cytotrace2", raster = TRUE, pt.size = 2, order = TRUE) +
    ggplot2::scale_colour_viridis_c() + ggplot2::ggtitle("CytoTRACE2 potency")
  p4 <- Seurat::DimPlot(mal, group.by = "progenitor_like", raster = TRUE, pt.size = 2,
                        cols = c(`FALSE` = "grey85", `TRUE` = "#D55E00")) + ggplot2::ggtitle("Progenitor-like call")
  p5 <- ggplot2::ggplot(mal[[]], ggplot2::aes(prog_group, cytotrace2, fill = prog_group)) +
    ggplot2::geom_violin(scale = "width", linewidth = 0.2) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = c(other = "grey80", progenitor = "#D55E00"), guide = "none") +
    ggplot2::facet_wrap(~ patient, nrow = 2) + ggplot2::theme_bw(base_size = 8) +
    ggplot2::labs(x = NULL, y = "CytoTRACE2 score")
  p6 <- Seurat::DotPlot(mal, features = intersect(markers, rownames(mal)), group.by = "malignant_state")
  patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 2)
}
