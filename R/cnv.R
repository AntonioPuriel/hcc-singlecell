# Malignant-cell identification with inferCNV ------------------------------------
#
# Run per patient: copy-number profiles are patient-specific, and per-patient runs
# keep memory within the 2 GB-per-CPU limit. Observations = that patient's
# epithelial cells; reference = up to 1,000 of that patient's immune/stromal cells.
# Each cell gets two numbers from the denoised inferCNV matrix (Tirosh et al. 2016;
# Puram et al. 2017):
#   cnv_score: mean squared deviation from the reference (how much CNV signal)
#   cnv_cor:   correlation with the patient's consensus tumour profile (the mean of
#              the 10% highest-scoring epithelial cells)
# A cell is called malignant when both exceed the 99th percentile of the patient's
# reference cells. Non-tumour liver epithelium serves as a check: it should be
# called malignant rarely.

epithelial_patients <- function(seu, min_cells = 50) {
  tab <- table(seu$patient[seu$compartment == "Epithelial"])
  sort(names(tab)[tab >= min_cells])
}

cnv_scores_from_expr <- function(expr, obs, ref, top_frac = 0.1) {
  dev <- expr - 1
  score <- colMeans(dev^2)
  top <- obs[order(score[obs], decreasing = TRUE)][seq_len(max(10, ceiling(top_frac * length(obs))))]
  profile <- rowMeans(dev[, top, drop = FALSE])
  cor_v <- as.numeric(stats::cor(dev, profile))
  data.frame(cell = colnames(expr), cnv_score = score, cnv_cor = cor_v,
             is_ref = colnames(expr) %in% ref)
}

run_infercnv_patient <- function(seu, patient, gene_order_file, n_ref_max = 1000,
                                 seed = 1234, out_root = "tmp/infercnv") {
  md <- seu[[]]
  in_p <- md$patient == patient
  obs <- rownames(md)[in_p & md$compartment == "Epithelial"]
  pool <- rownames(md)[in_p & md$compartment %in% c("T_NK", "Myeloid", "Endothelial", "Fibroblast")]
  set.seed(seed)
  ref <- if (length(pool) > n_ref_max) sample(pool, n_ref_max) else pool

  counts <- SeuratObject::LayerData(seu, assay = "RNA", layer = "counts")[, c(obs, ref)]
  ann <- data.frame(group = c(paste0("obs_", md[obs, "site"]), rep("ref", length(ref))),
                    row.names = c(obs, ref))
  go <- utils::read.table(gene_order_file, sep = "\t", row.names = 1, header = FALSE,
                          stringsAsFactors = FALSE)
  go <- go[!duplicated(rownames(go)), ]

  obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = counts, annotations_file = ann,
                                        gene_order_file = go, ref_group_names = "ref")
  out_dir <- file.path(out_root, patient)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)   # intermediates are large
  args <- list(infercnv_obj = obj, cutoff = 0.1, out_dir = out_dir, cluster_by_groups = TRUE,
               denoise = TRUE, HMM = FALSE, analysis_mode = "samples",
               num_threads = n_workers(), no_plot = TRUE, no_prelim_plot = TRUE,
               plot_steps = FALSE, write_expr_matrix = FALSE, write_phylo = FALSE,
               save_rds = FALSE, save_final_rds = FALSE, resume_mode = FALSE)
  args <- args[names(args) %in% names(formals(infercnv::run))]  # robust to version differences
  obj <- do.call(infercnv::run, args)

  res <- cnv_scores_from_expr(obj@expr.data, obs = intersect(obs, colnames(obj@expr.data)),
                              ref = ref)
  res$patient <- patient
  res$site <- md[res$cell, "site"]
  res$n_genes <- nrow(obj@expr.data)
  res
}

call_malignant <- function(cnv, q = 0.99) {
  do.call(rbind, lapply(split(cnv, cnv$patient), function(d) {
    thr_s <- stats::quantile(d$cnv_score[d$is_ref], q)
    thr_c <- stats::quantile(d$cnv_cor[d$is_ref], q)
    hi_s <- d$cnv_score > thr_s; hi_c <- d$cnv_cor > thr_c
    d$cnv_call <- ifelse(d$is_ref, "reference",
                  ifelse(hi_s & hi_c, "malignant",
                  ifelse(!hi_s & !hi_c, "non-malignant", "unresolved")))
    d$thr_score <- thr_s; d$thr_cor <- thr_c
    d
  }))
}

add_cnv <- function(seu, calls) {
  obs <- calls[!calls$is_ref, ]
  i <- match(colnames(seu), obs$cell)
  seu$cnv_score <- obs$cnv_score[i]
  seu$cnv_cor   <- obs$cnv_cor[i]
  seu$cnv_call  <- obs$cnv_call[i]
  seu
}

cnv_summary <- function(calls) {
  obs <- calls[!calls$is_ref, ]
  tab <- as.data.frame.matrix(table(paste(obs$patient, obs$site, sep = "|"), obs$cnv_call))
  out <- cbind(do.call(rbind, strsplit(rownames(tab), "|", fixed = TRUE)), tab)
  colnames(out)[1:2] <- c("patient", "site")
  out$n_epithelial <- rowSums(tab)
  out$pct_malignant <- round(100 * out$malignant / out$n_epithelial, 1)
  rownames(out) <- NULL
  out[order(out$patient, out$site), ]
}

plot_cnv <- function(calls) {
  calls$group <- ifelse(calls$is_ref, "reference (immune/stromal)", paste("epithelial,", calls$site))
  thr <- unique(calls[, c("patient", "thr_score", "thr_cor")])
  ggplot2::ggplot(calls, ggplot2::aes(cnv_score, cnv_cor, colour = group)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.4) +
    ggplot2::geom_vline(data = thr, ggplot2::aes(xintercept = thr_score), linetype = 2) +
    ggplot2::geom_hline(data = thr, ggplot2::aes(yintercept = thr_cor), linetype = 2) +
    ggplot2::facet_wrap(~ patient, scales = "free_x") +
    ggplot2::scale_colour_manual(values = c("reference (immune/stromal)" = "grey60",
      "epithelial, Normal" = "#009E73", "epithelial, Tumor" = "#D55E00",
      "epithelial, PVTT" = "#CC79A7", "epithelial, Lymph" = "#0072B2")) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    ggplot2::labs(x = "CNV score", y = "correlation with tumour CNV profile", colour = NULL) +
    ggplot2::theme_bw(base_size = 9)
}
