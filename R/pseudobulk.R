# Pseudobulk differential expression -----------------------------------------------
#
# Cells from one sample are not independent replicates. Counts are summed per
# sample x group, and DESeq2 is fitted with the patient as blocking factor
# (~ patient + group), so every contrast is made within patients.

pseudobulk <- function(seu, group_col, cells = colnames(seu), min_cells = 20) {
  md <- seu[[]][cells, ]
  key <- paste(md$sample, md[[group_col]], sep = "__")
  keep <- names(which(table(key) >= min_cells))
  cells <- cells[key %in% keep]; key <- key[key %in% keep]
  f <- factor(key)
  ind <- Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                              dimnames = list(cells, levels(f)))
  counts <- SeuratObject::LayerData(seu, assay = "RNA", layer = "counts")[, cells] %*% ind
  idx <- match(levels(f), key)
  col <- data.frame(pb = levels(f), sample = md[cells[idx], "sample"],
                    patient = md[cells[idx], "patient"], group = sub(".*__", "", levels(f)),
                    n_cells = as.integer(table(f)[levels(f)]), row.names = levels(f))
  list(counts = as.matrix(counts), coldata = col)
}

#' Paired DE: keeps patients that have both levels, ~ patient + group.
de_paired <- function(pb, test, reference, min_count = 10) {
  cd <- pb$coldata[pb$coldata$group %in% c(test, reference), ]
  both <- names(which(tapply(cd$group, cd$patient, function(g) all(c(test, reference) %in% g))))
  cd <- cd[cd$patient %in% both, ]
  if (length(both) < 3) stop("Fewer than 3 patients with both ", test, " and ", reference)
  counts <- round(pb$counts[, rownames(cd)])
  counts <- counts[rowSums(counts >= min_count) >= 3, ]
  cd$group <- factor(ifelse(cd$group == test, "test", "reference"), levels = c("reference", "test"))
  cd$patient <- factor(cd$patient)
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, design = ~ patient + group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- as.data.frame(DESeq2::results(dds, name = "group_test_vs_reference"))
  res$gene <- rownames(res)
  res$n_patients <- length(both)
  res$contrast <- paste(test, "vs", reference)
  res[order(res$pvalue), c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                           "n_patients", "contrast")]
}

hallmark_gsea <- function(de, seed = 1234) {
  hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  sets <- split(hs$gene_symbol, hs$gs_name)
  stats <- stats::setNames(de$stat, de$gene)
  stats <- sort(stats[!is.na(stats)], decreasing = TRUE)
  set.seed(seed)
  res <- fgsea::fgsea(pathways = sets, stats = stats, minSize = 15, maxSize = 500)
  res$leadingEdge <- vapply(res$leadingEdge, function(g) paste(utils::head(g, 10), collapse = ","), "")
  res$contrast <- de$contrast[1]
  as.data.frame(res[order(res$padj), ])
}

plot_gsea <- function(gsea_list, top = 12) {
  d <- do.call(rbind, lapply(gsea_list, function(g) utils::head(g[g$padj < 0.25, ], top)))
  if (!nrow(d)) d <- do.call(rbind, lapply(gsea_list, utils::head, 5))
  d$pathway <- sub("HALLMARK_", "", d$pathway)
  ggplot2::ggplot(d, ggplot2::aes(NES, stats::reorder(pathway, NES), fill = -log10(padj))) +
    ggplot2::geom_col() + ggplot2::facet_wrap(~ contrast, scales = "free_y") +
    ggplot2::scale_fill_viridis_c() + ggplot2::labs(y = NULL, x = "normalised enrichment score") +
    ggplot2::theme_bw(base_size = 9)
}
