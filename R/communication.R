# Cell-cell communication and candidate prioritisation ----------------------------------
#
# LIANA consensus (Dimitrov et al. 2022) on tumour-site cells, with progenitor-like
# and other malignant cells, the myeloid states and the remaining compartments as
# groups. Methods without permutations (NATMI, Connectome, log2FC, SingleCellSignalR)
# keep runtime and memory manageable on ~35k cells; their ranks are aggregated.
#
# Prioritisation of progenitor-like <-> tumour-associated macrophage (TAM)
# interactions. Revised after the second run, where the third evidence line
# (pseudobulk DE of the tumour-side gene) was noise: almost every gene had
# padj ~ 1, and the top-ranked gene (CD24) was one of the genes that define the
# group. The score now uses two ranks:
#   1. LIANA aggregate rank of the interaction
#   2. specificity: how much better the same ligand-receptor pair ranks with
#      progenitor-like cells than with other malignant cells on the tumour side
# Only interactions with positive specificity are kept. DE of the tumour-side
# gene is reported as annotation, and genes of the defining programme are
# flagged (circular by construction). The score orders hypotheses for
# validation; it is not a statistical test.

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
  g[in_mal] <- ifelse(mal$progenitor_like[match(colnames(seu)[in_mal], colnames(mal))],
                      "Progenitor-like tumour", "Other tumour")
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

prioritise_mediators <- function(liana_res, de_prog, targets = tam_states,
                                 defining = progenitor_markers, top = 40) {
  lr <- liana_res
  lr$pair <- paste(lr$ligand.complex, lr$receptor.complex, sep = " -> ")
  prog <- "Progenitor-like tumour"
  other <- lr[lr$source == "Other tumour" | lr$target == "Other tumour", ]
  out <- list()
  for (dir in c("prog_to_tam", "tam_to_prog")) {
    if (dir == "prog_to_tam") {
      d <- lr[lr$source == prog & lr$target %in% targets, ]
      ref <- other[other$source == "Other tumour", c("target", "pair", "aggregate_rank")]
      d <- merge(d, ref, by = c("target", "pair"), all.x = TRUE, suffixes = c("", "_other"))
      d$tumour_gene <- d$ligand.complex
    } else {
      d <- lr[lr$target == prog & lr$source %in% targets, ]
      ref <- other[other$target == "Other tumour", c("source", "pair", "aggregate_rank")]
      d <- merge(d, ref, by = c("source", "pair"), all.x = TRUE, suffixes = c("", "_other"))
      d$tumour_gene <- d$receptor.complex
    }
    if (!nrow(d)) next
    d$direction <- dir
    out[[dir]] <- d
  }
  d <- do.call(rbind, out)
  d$aggregate_rank_other[is.na(d$aggregate_rank_other)] <- 1
  d$specificity <- log10(d$aggregate_rank_other) - log10(d$aggregate_rank)   # > 0: stronger with progenitor-like
  d <- d[d$specificity > 0, ]
  subunits <- strsplit(d$tumour_gene, "_", fixed = TRUE)
  d$defining_gene <- vapply(subunits, function(g) any(g %in% defining), logical(1))
  i <- match(vapply(subunits, `[`, "", 1), de_prog$gene)                     # complexes: first subunit
  d$tumour_gene_log2FC <- de_prog$log2FoldChange[i]
  d$tumour_gene_padj <- de_prog$padj[i]
  d$rank_liana <- rank(d$aggregate_rank)
  d$rank_specificity <- rank(-d$specificity)
  d$priority_score <- (d$rank_liana + d$rank_specificity) / 2
  d <- d[order(d$priority_score), c("direction", "source", "target", "ligand.complex", "receptor.complex",
                                    "aggregate_rank", "aggregate_rank_other", "specificity",
                                    "tumour_gene", "defining_gene", "tumour_gene_log2FC",
                                    "tumour_gene_padj", "priority_score")]
  utils::head(d, top)
}

plot_mediators <- function(prio, n = 25) {
  d <- utils::head(prio, n)
  d$pair <- paste(d$ligand.complex, "->", d$receptor.complex)
  d$partner <- ifelse(d$direction == "prog_to_tam", d$target, d$source)
  d$pair <- ifelse(d$defining_gene, paste0(d$pair, " *"), d$pair)
  d$pair <- factor(d$pair, levels = rev(unique(d$pair)))
  ggplot2::ggplot(d, ggplot2::aes(partner, pair, size = -log10(aggregate_rank), colour = specificity)) +
    ggplot2::geom_point() + ggplot2::facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_gradient(low = "grey75", high = "#D55E00") +
    ggplot2::labs(x = "TAM state", y = NULL, size = "-log10 LIANA rank",
                  colour = "specificity\n(vs other tumour)",
                  caption = "* tumour-side gene belongs to the programme that defines the progenitor-like group") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
