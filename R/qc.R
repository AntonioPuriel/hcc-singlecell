# Quality control --------------------------------------------------------------
#
# The GEO matrix was already filtered by the authors (71,915 cells), so QC here
# is a re-assessment, not a first pass.
#
# A first run with per-sample MAD thresholds on both tails (UMIs and % mito)
# removed 15% of cells, and the removal was composition-dependent: high-UMI
# outliers were almost exclusive to non-tumour samples (8.8% vs 0.8% in tumour).
# Per-sample MADs assume one unimodal population, but liver samples mix small
# immune cells with large, mitochondria-rich hepatocytes, so the upper tail is a
# cell type, not a technical artefact. The final rules are therefore:
#   - low-quality cells: lower-tail MAD outliers on log UMIs / log genes, per sample
#   - % mitochondrial: a single absolute cap (no per-sample MAD)
#   - no upper UMI filter: multiplets are handled by scDblFinder
# The MAD-based upper-tail flags are still computed (prefixed `diag_`) and
# tabulated by author cell type to document why they were not used.

add_qc_metrics <- function(seu) {
  seu[["percent_mt"]]   <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent_ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu[["log10_umi"]]    <- log10(seu$nCount_RNA)
  seu[["log10_genes"]]  <- log10(seu$nFeature_RNA)
  seu
}

# TRUE where x lies more than `nmads` MADs from the median of its group.
mad_outlier <- function(x, group, nmads = 3, type = c("both", "lower", "higher")) {
  type <- match.arg(type)
  out <- logical(length(x))
  for (g in unique(group)) {
    i <- which(group == g)
    med <- stats::median(x[i]); dev <- stats::mad(x[i])
    lo <- x[i] < med - nmads * dev
    hi <- x[i] > med + nmads * dev
    out[i] <- switch(type, both = lo | hi, lower = lo, higher = hi)
  }
  out
}

flag_qc_outliers <- function(seu, nmads = 3, max_mt = 30) {
  md <- seu[[]]
  seu$qc_low_umi   <- mad_outlier(md$log10_umi,   md$sample, nmads, "lower")
  seu$qc_low_genes <- mad_outlier(md$log10_genes, md$sample, nmads, "lower")
  seu$qc_high_mt   <- md$percent_mt > max_mt
  seu$qc_fail <- seu$qc_low_umi | seu$qc_low_genes | seu$qc_high_mt
  # Diagnostics only (not used for filtering, see header)
  seu$diag_mad_high_umi <- mad_outlier(md$log10_umi,  md$sample, nmads, "higher")
  seu$diag_mad_high_mt  <- mad_outlier(md$percent_mt, md$sample, nmads, "higher")
  seu
}

#' Doublet detection with scDblFinder, run per sample (doublets form within a
#' capture, never across samples). Returns a data frame keyed by cell.
#' Workers are read from the SLURM allocation at run time (not tracked by targets)
#' and capped: CPUs are requested mainly to obtain memory (MaxMemPerCPU = 2 GB)
#' and each worker holds its own copy of a sample.
n_workers <- function(max_workers = 6L) {
  min(max_workers, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
}

run_doublets <- function(seu, seed = 1234, workers = n_workers()) {
  set.seed(seed)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"))
  )
  bp <- if (workers > 1L) BiocParallel::MulticoreParam(workers, RNGseed = seed)
        else BiocParallel::SerialParam(RNGseed = seed)
  sce <- scDblFinder::scDblFinder(sce, samples = seu$sample, BPPARAM = bp)
  data.frame(cell = colnames(sce),
             doublet_score = sce$scDblFinder.score,
             doublet_class = as.character(sce$scDblFinder.class))
}

add_doublets <- function(seu, doublets) {
  i <- match(colnames(seu), doublets$cell)
  seu$doublet_score <- doublets$doublet_score[i]
  seu$doublet_class <- doublets$doublet_class[i]
  seu
}

apply_qc <- function(seu) {
  keep <- !seu$qc_fail & seu$doublet_class == "singlet"
  subset(seu, cells = colnames(seu)[keep])
}

#' Cells before/after QC per sample, with the reason for removal.
qc_summary <- function(seu_flagged) {
  md <- seu_flagged[[]]
  do.call(rbind, lapply(split(md, md$sample), function(d) data.frame(
    sample = d$sample[1], patient = d$patient[1], site = d$site[1],
    cells_in = nrow(d),
    low_umi = sum(d$qc_low_umi), low_genes = sum(d$qc_low_genes), high_mt = sum(d$qc_high_mt),
    doublets = sum(d$doublet_class == "doublet"),
    cells_out = sum(!d$qc_fail & d$doublet_class == "singlet"),
    median_umi = stats::median(d$nCount_RNA),
    median_genes = stats::median(d$nFeature_RNA),
    median_pct_mt = round(stats::median(d$percent_mt), 2)
  )))
}

#' How the rejected upper-tail MAD rules would have hit each author cell type,
#' compared with the rules actually applied. Documents the QC decision.
qc_flags_by_celltype <- function(seu_flagged, label_col) {
  md <- seu_flagged[[]]
  if (!label_col %in% colnames(md)) return(data.frame(note = paste("column", label_col, "not found")))
  pct <- function(x) round(100 * mean(x), 1)
  out <- do.call(rbind, lapply(split(md, md[[label_col]]), function(d) data.frame(
    cell_type = d[[label_col]][1], n_cells = nrow(d),
    median_umi = stats::median(d$nCount_RNA), median_pct_mt = round(stats::median(d$percent_mt), 2),
    pct_diag_mad_high_umi = pct(d$diag_mad_high_umi),
    pct_diag_mad_high_mt  = pct(d$diag_mad_high_mt),
    pct_removed_qc        = pct(d$qc_fail),
    pct_doublet           = pct(d$doublet_class == "doublet"))))
  out[order(-out$n_cells), ]
}

plot_qc <- function(seu_flagged) {
  md <- seu_flagged[[]]
  md$status <- ifelse(md$doublet_class == "doublet", "doublet",
                      ifelse(md$qc_fail, "QC outlier", "kept"))
  md$sample <- factor(md$sample, levels = sort(unique(md$sample)))
  base <- function(y, lab) {
    ggplot2::ggplot(md, ggplot2::aes(sample, .data[[y]])) +
      ggplot2::geom_violin(fill = "grey85", colour = NA, scale = "width") +
      ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA) +
      ggplot2::facet_grid(~ site, scales = "free_x", space = "free_x") +
      ggplot2::labs(x = NULL, y = lab) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  p1 <- base("log10_umi", "log10 UMIs")
  p2 <- base("log10_genes", "log10 genes")
  p3 <- base("percent_mt", "% mitochondrial")
  p4 <- ggplot2::ggplot(md, ggplot2::aes(sample, fill = status)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::scale_fill_manual(values = c(kept = "grey70", `QC outlier` = "#D55E00",
                                          doublet = "#0072B2")) +
    ggplot2::facet_grid(~ site, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = NULL, y = "fraction of cells", fill = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 1)
}
