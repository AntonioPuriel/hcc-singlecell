# hcc-singlecell pipeline ------------------------------------------------------
# Run from the project root:  targets::tar_make()   (on the cluster: slurm/10_run_pipeline.sh)
#
# Parameters are stored as individual targets so that changing one of them only
# invalidates the steps that use it. Computing resources (CPUs, workers) are read
# at run time inside functions and are deliberately NOT tracked: changing the
# job size must never trigger a recomputation.

library(targets)
library(tarchetypes)

tar_source("R")

tar_option_set(
  packages = c("Matrix", "data.table", "ggplot2"),
  format   = format_qs2,
  seed     = 1234,
  memory   = "transient",      # drop large objects from RAM once used
  garbage_collection = TRUE
)

list(
  # ---- Parameters -------------------------------------------------------------
  # Column names in the GEO metadata (checked by read_metadata())
  tar_target(meta_cols, list(cell = "Cell", sample = "sample", patient = "patient", site = "site")),
  tar_target(qc_params, list(nmads = 3, max_mt = 30)),
  tar_target(author_label, "celltype"),   # author annotation column, used only for comparison

  # ---- 1. Ingestion --------------------------------------------------------
  tar_target(counts_file,   "data/raw/GSE149614_HCC.scRNAseq.S71915.count.txt.gz", format = "file"),
  tar_target(metadata_file, "data/raw/GSE149614_HCC.metadata.updated.txt.gz",     format = "file"),
  tar_target(counts,   read_dense_counts(counts_file)),
  tar_target(metadata, read_metadata(metadata_file, meta_cols)),
  tar_target(seu_raw,  build_seurat(counts, metadata, meta_cols)),
  tar_target(ingest_tsv,
             write_tsv(ingest_summary(counts, metadata, seu_raw, meta_cols),
                       "results/tables/01_ingest_summary.tsv"),
             format = "file"),
  tar_target(design_tsv,
             write_tsv(design_table(seu_raw), "results/tables/01_cells_per_patient_site.tsv"),
             format = "file"),

  # ---- 2. Quality control ----------------------------------------------------
  tar_target(seu_metrics, flag_qc_outliers(add_qc_metrics(seu_raw),
                                           nmads = qc_params$nmads, max_mt = qc_params$max_mt)),
  tar_target(doublets, run_doublets(seu_metrics, seed = 1234)),
  tar_target(seu_flagged, add_doublets(seu_metrics, doublets)),
  tar_target(seu_qc, apply_qc(seu_flagged)),
  tar_target(qc_tsv,
             write_tsv(qc_summary(seu_flagged), "results/tables/02_qc_summary.tsv"),
             format = "file"),
  tar_target(qc_celltype_tsv,
             write_tsv(qc_flags_by_celltype(seu_flagged, author_label),
                       "results/tables/02_qc_flags_by_celltype.tsv"),
             format = "file"),
  tar_target(qc_png,
             save_plot(plot_qc(seu_flagged), "results/figures/02_qc_per_sample.png",
                       width = 12, height = 13),
             format = "file"),

  # ---- 3. Normalisation, integration, clustering -------------------------------
  tar_target(seu_pca, normalise_and_pca(seu_qc)),
  tar_target(seu_int, cluster_cells(integrate_harmony(seu_pca))),
  tar_target(lisi_tsv,
             write_tsv(integration_metrics(seu_int, label_col = author_label),
                       "results/tables/03_integration_lisi.tsv"),
             format = "file"),
  tar_target(integration_png,
             save_plot(plot_integration(seu_int), "results/figures/03_integration_umap.png",
                       width = 14, height = 12),
             format = "file"),

  # ---- 4. Annotation -------------------------------------------------------------
  tar_target(singler, run_singler(seu_int)),
  tar_target(seu_annot, add_singler(seu_int, singler)),
  tar_target(markers_tsv,
             write_tsv(cluster_markers(seu_annot), "results/tables/04_cluster_markers.tsv"),
             format = "file"),
  tar_target(annotation_tsv,
             write_tsv(annotation_table(seu_annot, author_label),
                       "results/tables/04_cluster_annotation.tsv"),
             format = "file"),
  tar_target(annotation_png,
             save_plot(plot_annotation(seu_annot, author_label),
                       "results/figures/04_annotation_umap.png", width = 16, height = 12),
             format = "file"),
  tar_target(dotplot_png,
             save_plot(plot_marker_dotplot(seu_annot), "results/figures/04_marker_dotplot.png",
                       width = 16, height = 9),
             format = "file"),

  # ---- 5. Compartments ------------------------------------------------------------
  # Cluster numbers refer to the first-pass Leiden clustering (res 0.5, seed 1234);
  # check_compartments() stops the pipeline if they no longer match.
  tar_target(compartment_map, list(
    T_NK        = c(0, 2, 9, 12),
    Myeloid     = c(1, 5, 7, 18),
    pDC         = 22,
    B_plasma    = c(14, 15),
    Endothelial = c(8, 16),
    Fibroblast  = 10,
    Epithelial  = c(3, 4, 6, 11, 13, 17, 21),
    Flagged     = c(19, 20)    # possible residual doublets / ambient RNA, excluded downstream
  )),
  tar_target(seu_comp, add_compartments(seu_annot, compartment_map)),
  tar_target(compartment_tsv,
             write_tsv(check_compartments(seu_comp, author_label, list(
               T_NK = "T/NK", Myeloid = "Myeloid", B_plasma = "B", Endothelial = "Endothelial",
               Fibroblast = "Fibroblast", Epithelial = "Hepatocyte")),
               "results/tables/05_compartment_check.tsv"),
             format = "file"),

  # ---- 6. Malignant cells (inferCNV, per patient) ---------------------------------
  tar_target(gene_order_file, "data/ref/gene_order_hg38.txt", format = "file"),
  tar_target(cnv_patients, epithelial_patients(seu_comp)),
  tar_target(cnv, run_infercnv_patient(seu_comp, cnv_patients, gene_order_file),
             pattern = map(cnv_patients)),
  tar_target(cnv_calls, call_malignant(cnv)),
  tar_target(seu_cnv, add_cnv(seu_comp, cnv_calls)),
  tar_target(cnv_tsv, write_tsv(cnv_summary(cnv_calls), "results/tables/06_cnv_calls_by_sample.tsv"),
             format = "file"),
  tar_target(cnv_png, save_plot(plot_cnv(cnv_calls), "results/figures/06_cnv_scores.png",
                                width = 14, height = 10),
             format = "file"),

  # ---- 7. Myeloid states -------------------------------------------------------------
  tar_target(seu_myeloid, subcluster_myeloid(seu_cnv)),
  tar_target(myeloid_marker_table, myeloid_markers(seu_myeloid)),
  tar_target(myeloid_markers_tsv,
             write_tsv(myeloid_marker_table, "results/tables/07_myeloid_markers.tsv"),
             format = "file"),
  tar_target(myeloid_summary_tsv,
             write_tsv(myeloid_summary(seu_myeloid, myeloid_marker_table),
                       "results/tables/07_myeloid_states.tsv"),
             format = "file"),
  tar_target(myeloid_png, save_plot(plot_myeloid(seu_myeloid), "results/figures/07_myeloid_states.png",
                                    width = 16, height = 13),
             format = "file")
)
