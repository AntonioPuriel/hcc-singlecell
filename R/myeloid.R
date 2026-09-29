# Myeloid states ------------------------------------------------------------------
#
# Myeloid cells are re-analysed on their own (new HVGs, PCA, Harmony by patient,
# Leiden), because states inside a compartment are invisible at atlas resolution.
# "Immunosuppressive" is defined by an explicit gene programme, not by assumed
# identity: MDSC-like cells cannot be separated reliably from monocytes or
# neutrophils with scRNA-seq alone.

immunosuppression_genes <- c("CD274", "PDCD1LG2", "IL10", "TGFB1", "VEGFA", "SPP1", "TREM2",
                             "APOE", "CD163", "MRC1", "LGALS9", "IDO1", "ARG1", "CCL18",
                             "VSIG4", "SIGLEC10", "HAVCR2", "NT5E")

myeloid_panel <- c("FCN1", "S100A8", "VCAN", "LYZ", "CD14", "FCGR3A", "CD68", "C1QA", "APOE",
                   "SPP1", "TREM2", "CD163", "MARCO", "CD5L", "VSIG4", "CD1C", "CLEC9A",
                   "LAMP3", "CCR7", "CD274", "IL10", "VEGFA", "MKI67")

subcluster_myeloid <- function(seu, dims = 1:20, resolution = 0.6, seed = 1234,
                               batch_correct = TRUE, suppression = immunosuppression_genes) {
  sub <- subset(seu, cells = colnames(seu)[seu$compartment == "Myeloid"])
  sub <- Seurat::DietSeurat(sub, layers = "counts")
  sub <- normalise_and_pca(sub, n_hvg = 2000, n_pcs = max(dims), seed = seed)
  red <- "pca"
  if (batch_correct) {
    set.seed(seed)
    sub <- harmony::RunHarmony(sub, group.by.vars = "patient", reduction.use = "pca",
                               dims.use = dims, project.dim = FALSE, verbose = FALSE)
    red <- "harmony"
  }
  sub <- Seurat::RunUMAP(sub, reduction = red, dims = dims, seed.use = seed, verbose = FALSE)
  sub <- Seurat::FindNeighbors(sub, reduction = red, dims = dims, verbose = FALSE)
  sub$myeloid_state <- leiden_clusters(sub, resolution, seed)
  Seurat::Idents(sub) <- "myeloid_state"
  genes <- intersect(suppression, rownames(sub))
  sub <- Seurat::AddModuleScore(sub, features = list(genes), name = "suppression_score", seed = seed,
                                ctrl = min(100L, floor(nrow(sub) / 50)))
  sub$suppression_score <- sub$suppression_score1; sub$suppression_score1 <- NULL
  sub
}

myeloid_markers <- function(sub, n_top = 25) {
  res <- presto::wilcoxauc(sub, group_by = "myeloid_state", seurat_assay = "RNA")
  res <- res[res$padj < 0.05 & res$logFC > 0.25, ]
  res <- res[order(res$group, -res$auc), ]
  do.call(rbind, lapply(split(res, res$group), utils::head, n_top))
}

#' Per state: size, patients, site mix, mean immunosuppression score, top markers,
#' and tumour enrichment = median over paired patients of log2(fraction of the
#' patient's tumour myeloid cells in the state / same fraction in non-tumour).
myeloid_summary <- function(sub, markers, pseudo = 0.005) {
  md <- sub[[]]
  paired <- names(which(tapply(md$site, md$patient, function(s) all(c("Tumor", "Normal") %in% s))))
  frac <- function(p, s) { d <- md[md$patient == p & md$site == s, ]
                           table(factor(d$myeloid_state, levels = levels(md$myeloid_state))) / nrow(d) }
  lfc <- sapply(paired, function(p) log2((frac(p, "Tumor") + pseudo) / (frac(p, "Normal") + pseudo)))
  top <- tapply(markers$feature, markers$group, function(g) paste(utils::head(g, 8), collapse = ", "))
  out <- do.call(rbind, lapply(levels(md$myeloid_state), function(s) {
    d <- md[md$myeloid_state == s, ]
    data.frame(state = s, n_cells = nrow(d), n_patients = length(unique(d$patient)),
               pct_tumour_sites = round(100 * mean(d$site != "Normal"), 1),
               median_log2_tumour_vs_normal = round(stats::median(lfc[s, ]), 2),
               n_paired_patients = length(paired),
               mean_suppression_score = round(mean(d$suppression_score), 3),
               top_markers = unname(top[s]))
  }))
  out[order(-out$mean_suppression_score), ]
}

plot_myeloid <- function(sub) {
  p1 <- Seurat::DimPlot(sub, group.by = "myeloid_state", label = TRUE, raster = TRUE, pt.size = 2) +
    Seurat::NoLegend() + ggplot2::ggtitle("Myeloid states")
  p2 <- Seurat::DimPlot(sub, group.by = "site", raster = TRUE, pt.size = 2) + ggplot2::ggtitle("Tissue site")
  p3 <- Seurat::FeaturePlot(sub, "suppression_score", raster = TRUE, pt.size = 2, order = TRUE) +
    ggplot2::scale_colour_viridis_c() + ggplot2::ggtitle("Immunosuppression programme")
  p4 <- Seurat::DotPlot(sub, features = intersect(myeloid_panel, rownames(sub)),
                        group.by = "myeloid_state") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))
  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2, heights = c(1, 1))
}
