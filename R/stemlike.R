# Stem-like tumour states -----------------------------------------------------------
#
# Cancer stem cells in HCC are a debated concept, and a score built on a handful of
# surface markers is weak evidence on its own. Two independent readouts are
# therefore combined, within malignant cells (inferCNV) only:
#   csc_score:     module score of reported HCC CSC markers
#   potency_rank:  percentile, within sample, of the number of detected genes;
#                  gene counts are the core signal of CytoTRACE (Gulati et al. 2020)
#                  and do not depend on any marker list
# A cell is "stem-like" when both are in the top quartile of its patient's
# malignant cells. The per-patient Spearman correlation between the two readouts
# is reported: if they do not agree, the state is not supported.

csc_markers <- c("EPCAM", "PROM1", "CD24", "KRT19", "CD44", "THY1", "ANPEP", "SOX9")

score_stemlike <- function(seu, markers = csc_markers, q = 0.75, seed = 1234) {
  mal <- subset(seu, cells = colnames(seu)[seu$cnv_call %in% "malignant"])
  mal <- Seurat::DietSeurat(mal, layers = "counts")
  mal <- Seurat::NormalizeData(mal, verbose = FALSE)
  genes <- intersect(markers, rownames(mal))
  mal <- Seurat::AddModuleScore(mal, features = list(genes), name = "csc", seed = seed,
                                ctrl = min(100L, floor(nrow(mal) / 50)))
  mal$csc_score <- mal$csc1; mal$csc1 <- NULL
  md <- mal[[]]
  mal$potency_rank <- stats::ave(md$nFeature_RNA, md$sample, FUN = function(x) rank(x) / length(x))
  md <- mal[[]]
  hi <- function(x, g) stats::ave(x, g, FUN = function(v) v >= stats::quantile(v, q))
  mal$stem_like <- as.logical(hi(md$csc_score, md$patient)) & as.logical(hi(md$potency_rank, md$patient))
  mal$stem_group <- ifelse(mal$stem_like, "stem", "other")
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

stemlike_summary <- function(mal) {
  md <- mal[[]]
  do.call(rbind, lapply(split(md, md$patient), function(d) data.frame(
    patient = d$patient[1], n_malignant = nrow(d),
    n_stem_like = sum(d$stem_like), pct_stem_like = round(100 * mean(d$stem_like), 1),
    spearman_csc_vs_potency = round(stats::cor(d$csc_score, d$potency_rank, method = "spearman"), 2),
    sites = paste(sort(unique(d$site)), collapse = ","))))
}

plot_stemlike <- function(mal, markers = csc_markers) {
  p1 <- Seurat::DimPlot(mal, group.by = "patient", raster = TRUE, pt.size = 2) + ggplot2::ggtitle("Patient")
  p2 <- Seurat::FeaturePlot(mal, "csc_score", raster = TRUE, pt.size = 2, order = TRUE) +
    ggplot2::scale_colour_viridis_c() + ggplot2::ggtitle("CSC marker score")
  p3 <- Seurat::FeaturePlot(mal, "potency_rank", raster = TRUE, pt.size = 2, order = TRUE) +
    ggplot2::scale_colour_viridis_c() + ggplot2::ggtitle("Potency (gene-count rank)")
  p4 <- Seurat::DimPlot(mal, group.by = "stem_like", raster = TRUE, pt.size = 2,
                        cols = c(`FALSE` = "grey85", `TRUE` = "#D55E00")) + ggplot2::ggtitle("Stem-like call")
  p5 <- ggplot2::ggplot(mal[[]], ggplot2::aes(potency_rank, csc_score)) +
    ggplot2::geom_point(size = 0.2, alpha = 0.3) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "#D55E00") +
    ggplot2::facet_wrap(~ patient) + ggplot2::theme_bw(base_size = 8) +
    ggplot2::labs(x = "potency rank (within sample)", y = "CSC marker score")
  p6 <- Seurat::DotPlot(mal, features = intersect(markers, rownames(mal)), group.by = "malignant_state")
  patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 2)
}
