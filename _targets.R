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
  tar_target(qc_params, list(nmads = 3, max_mt = 50)),

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
  tar_target(qc_png,
             save_plot(plot_qc(seu_flagged), "results/figures/02_qc_per_sample.png",
                       width = 12, height = 13),
             format = "file")
)
