# Myeloid trajectories ---------------------------------------------------------------
#
# Slingshot (Street et al. 2018) on the Harmony embedding of the myeloid subclustering,
# rooted in classical monocytes. Only states that can plausibly derive from circulating
# monocytes are kept:
#   - Kupffer cells are excluded: they are yolk-sac-derived resident macrophages, so a
#     monocyte -> Kupffer curve would be an artefact of transcriptional similarity;
#   - dendritic cells (cDC1, cDC2, MHCII-high mac/DC) belong to a separate lineage;
#   - the stress-high monocyte state (dissociation signature) and artefact clusters are
#     excluded.
# Pseudotime orders cells by transcriptional similarity; it is not a time course and
# does not prove differentiation. Genes are ranked by their Spearman correlation with
# pseudotime, and each correlation is recomputed within each patient: a gene is only
# reported as consistent when most patients agree on the sign, so the ranking is not
# driven by one patient (cells are not independent observations).

trajectory_states <- c("Classical monocyte", "SPP1+ TAM", "C1Q+APOE+ TAM", "APOE+LGMN+ TAM",
                       "Cycling macrophage")

#' Slingshot lineages on the myeloid object. Returns the cell table (state, site,
#' patient, pseudotime per lineage, UMAP), the lineage paths and curves embedded in
#' the UMAP for plotting.
run_trajectory <- function(sub, label_check, labels = myeloid_state_labels,
                           keep = trajectory_states, start = "Classical monocyte",
                           reduction = "harmony", dims = 1:20) {
  stopifnot(all(label_check$found_in_top25))
  state <- unname(labels[as.character(sub$myeloid_state)])
  cells <- colnames(sub)[state %in% keep]
  emb <- SeuratObject::Embeddings(sub, reduction)[cells, dims]
  umap <- SeuratObject::Embeddings(sub, "umap")[cells, 1:2]
  cl <- state[match(cells, colnames(sub))]
  sds <- slingshot::slingshot(emb, clusterLabels = cl, start.clus = start, approx_points = 150)
  pt <- slingshot::slingPseudotime(sds)
  w <- slingshot::slingCurveWeights(sds)
  lin <- slingshot::slingLineages(sds)
  colnames(pt) <- colnames(w) <- paste0("lineage", seq_along(lin))
  # Each cell is assigned to the lineage with the highest curve weight; cells shared
  # by several lineages (the monocyte root) appear in all of them in the gene tests.
  md <- sub[[]][cells, ]
  cells_df <- data.frame(cell = cells, state = cl, site = md$site, patient = md$patient,
                         suppression_score = md$suppression_score,
                         umap_1 = umap[, 1], umap_2 = umap[, 2], pt, check.names = FALSE)
  curves <- slingshot::slingCurves(slingshot::embedCurves(sds, umap))
  curves_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
    s <- curves[[i]]$s[curves[[i]]$ord, , drop = FALSE]
    data.frame(lineage = paste0("lineage", i), umap_1 = s[, 1], umap_2 = s[, 2])
  }))
  list(cells = cells_df,
       lineages = data.frame(lineage = paste0("lineage", seq_along(lin)),
                             path = vapply(lin, paste, "", collapse = " -> ")),
       curves = curves_df)
}

#' Per lineage: endpoint, number of cells, fraction of tumour-site cells in the first
#' and last pseudotime decile, and the per-patient Spearman correlation of the
#' immunosuppression score with pseudotime (median and patients with rho > 0).
trajectory_summary <- function(traj, min_cells = 50) {
  d <- traj$cells
  do.call(rbind, lapply(traj$lineages$lineage, function(l) {
    x <- d[!is.na(d[[l]]), ]
    dec <- cut(rank(x[[l]], ties.method = "first"), 10, labels = FALSE)
    by_pat <- split(x, x$patient)
    by_pat <- by_pat[vapply(by_pat, nrow, 1L) >= min_cells]
    rho <- vapply(by_pat, function(p) stats::cor(p[[l]], p$suppression_score, method = "spearman"), 1)
    data.frame(lineage = l, path = traj$lineages$path[traj$lineages$lineage == l],
               n_cells = nrow(x),
               pct_tumour_sites_first_decile = round(100 * mean(x$site[dec == 1] != "Normal"), 1),
               pct_tumour_sites_last_decile = round(100 * mean(x$site[dec == 10] != "Normal"), 1),
               n_patients_tested = length(rho),
               median_rho_suppression = round(stats::median(rho), 2),
               n_patients_rho_positive = sum(rho > 0))
  }))
}

#' Genes associated with pseudotime, per lineage. Candidate genes are the myeloid HVGs.
#' Overall Spearman rho over all cells of the lineage plus per-patient rho; reports the
#' top genes in each direction with the number of patients that agree on the sign.
trajectory_genes <- function(sub, traj, n_top = 30, min_cells = 50, min_pct = 0.05) {
  genes <- SeuratObject::VariableFeatures(sub)
  expr <- SeuratObject::LayerData(sub, assay = "RNA", layer = "data")
  d <- traj$cells
  out <- lapply(traj$lineages$lineage, function(l) {
    x <- d[!is.na(d[[l]]), ]
    m <- expr[genes, x$cell, drop = FALSE]
    m <- m[Matrix::rowMeans(m > 0) >= min_pct, , drop = FALSE]
    rk <- function(cells) {               # Spearman = Pearson on ranks; dense per block of genes
      pt <- rank(x[[l]][match(cells, x$cell)])
      mm <- as.matrix(m[, cells, drop = FALSE])
      r <- t(apply(mm, 1, rank))
      as.vector(stats::cor(t(r), pt))
    }
    rho <- rk(x$cell)
    pats <- split(x$cell, x$patient)
    pats <- pats[lengths(pats) >= min_cells]
    rho_p <- vapply(pats, rk, numeric(nrow(m)))
    rho_p[is.na(rho_p)] <- 0
    res <- data.frame(lineage = l, gene = rownames(m), rho = round(rho, 3),
                      n_patients = length(pats),
                      n_patients_same_sign = rowSums(sign(rho_p) == sign(rho)))
    res <- res[!is.na(res$rho), ]
    up <- utils::head(res[order(-res$rho), ], n_top)
    down <- utils::head(res[order(res$rho), ], n_top)
    rbind(transform(up, direction = "increases"), transform(down, direction = "decreases"))
  })
  do.call(rbind, out)
}

plot_trajectory <- function(traj, genes, n_genes = 12) {
  d <- traj$cells
  d$pseudotime <- apply(d[, traj$lineages$lineage, drop = FALSE], 1, function(v) suppressWarnings(min(v, na.rm = TRUE)))
  p1 <- ggplot2::ggplot(d, ggplot2::aes(umap_1, umap_2, colour = state)) +
    ggplot2::geom_point(size = 0.3) +
    ggplot2::geom_path(data = traj$curves, ggplot2::aes(umap_1, umap_2, group = lineage),
                       colour = "black", linewidth = 0.8, inherit.aes = FALSE) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::theme_bw(base_size = 9) + ggplot2::ggtitle("Slingshot lineages (root: classical monocyte)")
  p2 <- ggplot2::ggplot(d, ggplot2::aes(umap_1, umap_2, colour = pseudotime)) +
    ggplot2::geom_point(size = 0.3) + ggplot2::scale_colour_viridis_c() +
    ggplot2::theme_bw(base_size = 9) + ggplot2::ggtitle("Pseudotime")
  long <- do.call(rbind, lapply(traj$lineages$lineage, function(l) {
    x <- d[!is.na(d[[l]]), ]
    data.frame(lineage = l, pseudotime = x[[l]], site = x$site, suppression_score = x$suppression_score)
  }))
  p3 <- ggplot2::ggplot(long, ggplot2::aes(pseudotime, suppression_score)) +
    ggplot2::geom_point(ggplot2::aes(colour = site), size = 0.2, alpha = 0.3) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "black", formula = y ~ x) +
    ggplot2::facet_wrap(~ lineage, scales = "free_x") +
    ggplot2::theme_bw(base_size = 9) + ggplot2::ggtitle("Immunosuppression programme along pseudotime")
  g <- genes[genes$n_patients_same_sign >= ceiling(0.75 * genes$n_patients), ]
  g <- do.call(rbind, lapply(split(g, list(g$lineage, g$direction)), utils::head, n_genes / 2))
  g$gene <- factor(paste(g$lineage, g$gene), levels = rev(paste(g$lineage, g$gene)))
  p4 <- ggplot2::ggplot(g, ggplot2::aes(rho, gene, fill = direction)) +
    ggplot2::geom_col() + ggplot2::facet_wrap(~ lineage, scales = "free_y") +
    ggplot2::scale_y_discrete(labels = function(x) sub("^\\S+ ", "", x)) +
    ggplot2::scale_fill_manual(values = c(increases = "#D55E00", decreases = "#0072B2")) +
    ggplot2::theme_bw(base_size = 8) + ggplot2::labs(x = "Spearman rho with pseudotime", y = NULL) +
    ggplot2::ggtitle("Top genes (sign agreed in >= 75% of patients)")
  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)
}
