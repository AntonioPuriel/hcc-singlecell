# Normalisation, integration across patients, clustering ------------------------
#
# Log-normalisation + PCA + Harmony (batch = patient). SCTransform was not used:
# differential expression is done on raw counts at the pseudobulk level, so the
# gain from SCT would be limited to the embedding, at a large memory cost on a
# cluster limited to 2 GB per CPU. Tissue site is NOT corrected: it is the
# biological signal of interest.

normalise_and_pca <- function(seu, n_hvg = 3000, n_pcs = 50, seed = 1234) {
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = n_hvg, verbose = FALSE)
  # Mitochondrial and ribosomal genes would drive PCs by technical/metabolic state
  hvg <- Seurat::VariableFeatures(seu)
  Seurat::VariableFeatures(seu) <- hvg[!grepl("^MT-|^RP[SL]", hvg)]
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  set.seed(seed)
  seu <- Seurat::RunPCA(seu, npcs = n_pcs, verbose = FALSE)
  seu[["RNA"]]$scale.data <- NULL   # ~2 GB, not needed downstream; saves disk in the target store
  seu
}

integrate_harmony <- function(seu, batch = "patient", dims = 1:30, seed = 1234) {
  set.seed(seed)
  # project.dim = FALSE: gene loadings are not projected (they would need the
  # scale.data layer, which is dropped after PCA to save memory and disk)
  seu <- harmony::RunHarmony(seu, group.by.vars = batch, reduction.use = "pca",
                             dims.use = dims, project.dim = FALSE, verbose = FALSE)
  # Unintegrated and integrated UMAPs, both kept for comparison
  seu <- Seurat::RunUMAP(seu, reduction = "pca", dims = dims, reduction.name = "umap_pca",
                         seed.use = seed, verbose = FALSE)
  seu <- Seurat::RunUMAP(seu, reduction = "harmony", dims = dims, reduction.name = "umap",
                         seed.use = seed, verbose = FALSE)
  seu
}

cluster_cells <- function(seu, dims = 1:30, resolutions = c(0.3, 0.5, 1), main = 0.5,
                          seed = 1234) {
  seu <- Seurat::FindNeighbors(seu, reduction = "harmony", dims = dims, verbose = FALSE)
  for (r in resolutions) seu[[paste0("RNA_snn_res.", r)]] <- leiden_clusters(seu, r, seed)
  seu$cluster <- seu[[paste0("RNA_snn_res.", main)]][, 1]
  Seurat::Idents(seu) <- "cluster"
  seu
}

#' Leiden clustering on the shared-nearest-neighbour graph. Uses Seurat's
#' algorithm 4 (leidenbase in Seurat >= 5.1); falls back to igraph's
#' implementation if that path needs Python's leidenalg. Clusters are relabelled
#' by decreasing size so that labels are stable and interpretable.
leiden_clusters <- function(seu, resolution, seed = 1234) {
  cl <- tryCatch({
    s <- Seurat::FindClusters(seu, resolution = resolution, algorithm = 4,
                              random.seed = seed, verbose = FALSE)
    as.character(s$seurat_clusters)
  }, error = function(e) {
    message("Seurat Leiden unavailable (", conditionMessage(e), "); using igraph")
    g <- igraph::graph_from_adjacency_matrix(seu@graphs$RNA_snn, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
    set.seed(seed)
    as.character(igraph::cluster_leiden(g, objective_function = "modularity",
                                        resolution = resolution, n_iterations = 10)$membership)
  })
  ord <- names(sort(table(cl), decreasing = TRUE))
  factor(match(cl, ord) - 1L, levels = seq_along(ord) - 1L)
}

# Integration assessment: local inverse Simpson index (LISI) -------------------
#
# For each cell, the effective number of distinct labels among its k nearest
# neighbours. Patient LISI (iLISI) should increase after integration (better
# mixing of patients); cell-type LISI (cLISI) should stay close to 1 (cell types
# not merged). Unweighted kNN version of Korsunsky et al. (2019).

lisi_knn <- function(embedding, labels, k = 30) {
  nn <- BiocNeighbors::findKNN(embedding, k = k, warn.ties = FALSE)$index
  lab <- as.integer(factor(labels))
  apply(nn, 1, function(i) {
    p <- tabulate(lab[i], nbins = max(lab)) / k
    1 / sum(p^2)
  })
}

integration_metrics <- function(seu, dims = 1:30, label_col = NULL, k = 30) {
  emb <- list(unintegrated = Seurat::Embeddings(seu, "pca")[, dims],
              harmony      = Seurat::Embeddings(seu, "harmony")[, dims])
  vars <- c(patient = "patient", site = "site")
  if (!is.null(label_col) && label_col %in% colnames(seu[[]])) vars <- c(vars, cell_type = label_col)
  out <- list()
  for (e in names(emb)) for (v in names(vars)) {
    l <- lisi_knn(emb[[e]], seu[[vars[[v]]]][, 1], k = k)
    out[[length(out) + 1]] <- data.frame(
      embedding = e, label = v, n_labels = length(unique(seu[[vars[[v]]]][, 1])),
      median_lisi = round(stats::median(l), 3), mean_lisi = round(mean(l), 3))
  }
  do.call(rbind, out)
}

# Annotation ---------------------------------------------------------------------

canonical_markers <- list(
  Hepatocyte  = c("ALB", "APOA1", "TTR", "SERPINA1"),
  Tumour      = c("GPC3", "AFP"),
  Cholangio   = c("EPCAM", "KRT19", "KRT7"),
  Endothelial = c("PECAM1", "VWF", "CLEC4G"),
  Fibro_HSC   = c("COL1A1", "DCN", "ACTA2", "RGS5"),
  T_NK        = c("CD3D", "CD8A", "CD4", "NKG7", "GNLY"),
  B_plasma    = c("MS4A1", "CD79A", "JCHAIN", "MZB1"),
  Myeloid     = c("LYZ", "CD14", "FCGR3A", "CD68", "C1QA", "CD163", "SPP1"),
  DC          = c("FCER1A", "CLEC9A", "LAMP3"),
  Mast        = c("TPSAB1", "KIT"),
  Cycling     = c("MKI67", "TOP2A")
)

#' Reference-based labels with SingleR (Human Primary Cell Atlas, main labels),
#' per cell, then summarised per cluster.
run_singler <- function(seu, workers = n_workers()) {
  ref <- celldex::HumanPrimaryCellAtlasData()
  bp <- if (workers > 1L) BiocParallel::MulticoreParam(workers) else BiocParallel::SerialParam()
  pred <- SingleR::SingleR(test = SeuratObject::LayerData(seu, layer = "data"),
                           ref = ref, labels = ref$label.main, BPPARAM = bp)
  data.frame(cell = rownames(pred), singler_label = pred$pruned.labels)
}

add_singler <- function(seu, singler) {
  seu$singler_label <- singler$singler_label[match(colnames(seu), singler$cell)]
  seu
}

#' Top markers per cluster (Wilcoxon AUC via presto; fast on ~70k cells).
cluster_markers <- function(seu, n_top = 20) {
  res <- presto::wilcoxauc(seu, group_by = "cluster", seurat_assay = "RNA")
  res <- res[res$padj < 0.05 & res$logFC > 0.25, ]
  res <- res[order(res$group, -res$auc), ]
  do.call(rbind, lapply(split(res, res$group), utils::head, n_top))
}

#' Per cluster: size, patients/sites represented, majority SingleR label and,
#' if available, majority author label with its agreement fraction.
annotation_table <- function(seu, author_col = NULL) {
  md <- seu[[]]
  maj <- function(x) { t <- sort(table(x, useNA = "no"), decreasing = TRUE)
                       if (length(t)) c(names(t)[1], round(t[[1]] / sum(t), 2)) else c(NA, NA) }
  do.call(rbind, lapply(split(md, md$cluster), function(d) {
    s <- maj(d$singler_label)
    row <- data.frame(cluster = d$cluster[1], n_cells = nrow(d),
                      n_patients = length(unique(d$patient)),
                      top_patient_frac = round(max(table(d$patient)) / nrow(d), 2),
                      sites = paste(names(sort(table(d$site), decreasing = TRUE)), collapse = ","),
                      singler_majority = s[1], singler_frac = as.numeric(s[2]))
    if (!is.null(author_col) && author_col %in% colnames(d)) {
      a <- maj(d[[author_col]])
      row$author_majority <- a[1]; row$author_frac <- as.numeric(a[2])
    }
    row
  }))
}

# Figures --------------------------------------------------------------------------

umap_panel <- function(seu, reduction, group, title, label = FALSE) {
  Seurat::DimPlot(seu, reduction = reduction, group.by = group, label = label,
                  raster = TRUE, pt.size = 1.5, label.size = 3) +
    ggplot2::ggtitle(title) + ggplot2::theme(legend.text = ggplot2::element_text(size = 7))
}

plot_integration <- function(seu) {
  patchwork::wrap_plots(
    umap_panel(seu, "umap_pca", "patient", "Unintegrated: patient"),
    umap_panel(seu, "umap",     "patient", "Harmony: patient"),
    umap_panel(seu, "umap_pca", "site",    "Unintegrated: tissue site"),
    umap_panel(seu, "umap",     "site",    "Harmony: tissue site"),
    ncol = 2)
}

plot_annotation <- function(seu, author_col = NULL) {
  p <- list(umap_panel(seu, "umap", "cluster", "Leiden clusters", label = TRUE) + Seurat::NoLegend(),
            umap_panel(seu, "umap", "singler_label", "SingleR (HPCA)"))
  if (!is.null(author_col) && author_col %in% colnames(seu[[]]))
    p <- c(p, list(umap_panel(seu, "umap", author_col, "Author annotation")))
  patchwork::wrap_plots(p, ncol = 2)
}

plot_marker_dotplot <- function(seu, markers = canonical_markers) {
  feats <- lapply(markers, intersect, rownames(seu))
  feats <- feats[lengths(feats) > 0]
  Seurat::DotPlot(seu, features = feats, group.by = "cluster") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
                   strip.text.x = ggplot2::element_text(angle = 90, size = 7))
}
