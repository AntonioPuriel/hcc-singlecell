# Cell-cell communication and candidate prioritisation ----------------------------------
#
# LIANA consensus (Dimitrov et al. 2022) on tumour-site cells, with stem-like and
# other malignant cells, the myeloid states and the remaining compartments as
# groups. Methods without permutations (NATMI, Connectome, log2FC, SingleCellSignalR)
# keep runtime and memory manageable on ~35k cells; their ranks are aggregated.
#
# Prioritisation of stem-like <-> tumour-associated macrophage (TAM) interactions
# combines three evidence lines, each converted to a rank:
#   1. LIANA aggregate rank of the interaction
#   2. specificity: the same ligand-receptor pair ranks better with stem-like
#      cells than with other malignant cells on the stem-like side
#   3. pseudobulk DE of the stem-like-side gene (stem-like vs other malignant,
#      paired within patients)
# The priority score is the mean of the three ranks (lower = stronger). It is a
# way to order hypotheses for validation, not a statistical test.

myeloid_state_labels <- c(`0` = "APOE+LGMN+ TAM", `1` = "Kupffer", `2` = "C1Q+APOE+ TAM",
                          `3` = "SPP1+ TAM", `4` = "Classical monocyte", `5` = "Stress-high monocyte",
                          `6` = "cDC2/mature DC", `7` = "MHCII-high mac/DC", `8` = "cDC1",
                          `9` = "Cycling macrophage", `10` = "artefact", `11` = "artefact")
tam_states <- c("SPP1+ TAM", "C1Q+APOE+ TAM", "APOE+LGMN+ TAM")

#' Guard: myeloid state numbers are only valid for this run; stop if the marker
#' that defines each named state is no longer highest in it.
check_myeloid_labels <- function(markers, expected = c(`3` = "SPP1", `1` = "CD5L", `4` = "FCN1",
                                                       `8` = "CLEC9A", `9` = "STMN1")) {
  ok <- vapply(names(expected), function(s)
    expected[[s]] %in% utils::head(markers$feature[markers$group == s], 25), logical(1))
  if (!all(ok)) stop("Myeloid state labels no longer match markers for state(s): ",
                     paste(names(expected)[!ok], collapse = ", "))
  data.frame(state = names(expected), marker = unname(expected), found_in_top25 = ok)
}

label_cell_groups <- function(seu, mal, sub, label_check, labels = myeloid_state_labels) {
  stopifnot(all(label_check$found_in_top25))
  g <- seu$compartment
  g[g == "Epithelial"] <- "Other epithelium"          # non-malignant or unresolved by inferCNV
  in_mal <- colnames(seu) %in% colnames(mal)
  g[in_mal] <- ifelse(mal$stem_like[match(colnames(seu)[in_mal], colnames(mal))],
                      "Stem-like tumour", "Other tumour")
  i <- match(colnames(seu), colnames(sub))
  my <- !is.na(i)
  g[my] <- unname(labels[as.character(sub$myeloid_state[i[my]])])
  seu$cell_group <- g
  seu
}

run_liana_tumour <- function(seu, sites = c("Tumor", "PVTT"), drop = c("Flagged", "Unassigned", "artefact",
                                                                        "Other epithelium"),
                             min_cells = 50) {
  keep <- seu$site %in% sites & !(seu$cell_group %in% drop)
  sub <- subset(seu, cells = colnames(seu)[keep])
  sizes <- table(sub$cell_group)
  sub <- subset(sub, cells = colnames(sub)[sub$cell_group %in% names(sizes)[sizes >= min_cells]])
  sub <- Seurat::NormalizeData(sub, verbose = FALSE)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = SeuratObject::LayerData(sub, layer = "counts"),
                  logcounts = SeuratObject::LayerData(sub, layer = "data")),
    colData = S4Vectors::DataFrame(cell_group = sub$cell_group))
  res <- liana::liana_wrap(sce, idents_col = "cell_group", assay.type = "logcounts",
                           method = c("natmi", "connectome", "logfc", "sca"),
                           resource = "Consensus", expr_prop = 0.1)
  as.data.frame(liana::liana_aggregate(res))
}

prioritise_mediators <- function(liana_res, de_stem, targets = tam_states, top = 40) {
  lr <- liana_res
  lr$pair <- paste(lr$ligand.complex, lr$receptor.complex, sep = " -> ")
  other <- lr[lr$source == "Other tumour" | lr$target == "Other tumour", ]
  out <- list()
  for (dir in c("stem_to_tam", "tam_to_stem")) {
    if (dir == "stem_to_tam") {
      d <- lr[lr$source == "Stem-like tumour" & lr$target %in% targets, ]
      ref <- other[other$source == "Other tumour", c("target", "pair", "aggregate_rank")]
      d <- merge(d, ref, by = c("target", "pair"), all.x = TRUE, suffixes = c("", "_other"))
      d$stem_gene <- d$ligand.complex
    } else {
      d <- lr[lr$target == "Stem-like tumour" & lr$source %in% targets, ]
      ref <- other[other$target == "Other tumour", c("source", "pair", "aggregate_rank")]
      d <- merge(d, ref, by = c("source", "pair"), all.x = TRUE, suffixes = c("", "_other"))
      d$stem_gene <- d$receptor.complex
    }
    if (!nrow(d)) next
    d$direction <- dir
    out[[dir]] <- d
  }
  d <- do.call(rbind, out)
  d$aggregate_rank_other[is.na(d$aggregate_rank_other)] <- 1
  d$specificity <- log10(d$aggregate_rank_other) - log10(d$aggregate_rank)   # > 0: stronger with stem-like
  first_gene <- sub("_.*", "", d$stem_gene)                                   # complexes: first subunit
  i <- match(first_gene, de_stem$gene)
  d$stem_gene_log2FC <- de_stem$log2FoldChange[i]
  d$stem_gene_padj <- de_stem$padj[i]
  d$stem_gene_stat <- de_stem$stat[i]
  d$rank_liana <- rank(d$aggregate_rank)
  d$rank_specificity <- rank(-d$specificity)
  d$rank_de <- rank(-ifelse(is.na(d$stem_gene_stat), -Inf, d$stem_gene_stat))
  d$priority_score <- (d$rank_liana + d$rank_specificity + d$rank_de) / 3
  d <- d[order(d$priority_score), c("direction", "source", "target", "ligand.complex", "receptor.complex",
                                    "aggregate_rank", "aggregate_rank_other", "specificity",
                                    "stem_gene", "stem_gene_log2FC", "stem_gene_padj", "priority_score")]
  utils::head(d, top)
}

plot_mediators <- function(prio, n = 25) {
  d <- utils::head(prio, n)
  d$pair <- paste(d$ligand.complex, "->", d$receptor.complex)
  d$partner <- ifelse(d$direction == "stem_to_tam", d$target, d$source)
  d$pair <- factor(d$pair, levels = rev(unique(d$pair)))
  ggplot2::ggplot(d, ggplot2::aes(partner, pair, size = -log10(aggregate_rank), colour = specificity)) +
    ggplot2::geom_point() + ggplot2::facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_gradient2(low = "grey60", mid = "grey80", high = "#D55E00", midpoint = 0) +
    ggplot2::labs(x = "TAM state", y = NULL, size = "-log10 LIANA rank",
                  colour = "specificity\n(vs other tumour)") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
