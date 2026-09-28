# Quality control --------------------------------------------------------------
#
# The GEO matrix was already filtered by the authors (71,915 cells), so QC here
# is a re-assessment, not a first pass: outliers are flagged per sample with
# robust (MAD-based) thresholds instead of fixed global cut-offs.
#
# Mitochondrial content: hepatocytes and many HCC tumour cells are
# mitochondria-rich. A fixed 10-20% cut-off would preferentially remove them, so
# only a per-sample upper MAD threshold plus a lenient absolute cap are used.

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

flag_qc_outliers <- function(seu, nmads = 3, max_mt = 50) {
  md <- seu[[]]
  seu$qc_low_umi   <- mad_outlier(md$log10_umi,   md$sample, nmads, "lower")
  seu$qc_high_umi  <- mad_outlier(md$log10_umi,   md$sample, nmads, "higher")
  seu$qc_low_genes <- mad_outlier(md$log10_genes, md$sample, nmads, "lower")
  seu$qc_high_mt   <- mad_outlier(md$percent_mt,  md$sample, nmads, "higher") |
                      md$percent_mt > max_mt
  seu$qc_fail <- seu$qc_low_umi | seu$qc_high_umi | seu$qc_low_genes | seu$qc_high_mt
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
    low_umi = sum(d$qc_low_umi), high_umi = sum(d$qc_high_umi),
    low_genes = sum(d$qc_low_genes), high_mt = sum(d$qc_high_mt),
    doublets = sum(d$doublet_class == "doublet"),
    cells_out = sum(!d$qc_fail & d$doublet_class == "singlet"),
    median_umi = stats::median(d$nCount_RNA),
    median_genes = stats::median(d$nFeature_RNA),
    median_pct_mt = round(stats::median(d$percent_mt), 2)
  )))
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
