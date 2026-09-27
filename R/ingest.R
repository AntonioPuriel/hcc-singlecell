# Ingestion: dense GEO count matrix -> sparse matrix -> Seurat object ---------

#' Read a gzipped dense genes x cells text matrix in row chunks.
#'
#' The GEO file is ~25k genes x ~72k cells stored as dense text; reading it in
#' one go needs >14 GB of RAM. Streaming blocks of genes and converting each
#' block to a sparse matrix keeps peak memory to a few GB.
#'
#' Handles both header layouts written by R/Python: header with one field less
#' than data rows (gene column unnamed) or header with a leading gene-column name.
read_dense_counts <- function(path, chunk_genes = 1000L, sep = "\t") {
  con <- gzfile(path, open = "r")
  on.exit(close(con), add = TRUE)

  header <- strsplit(readLines(con, n = 1L), sep, fixed = TRUE)[[1]]
  header <- gsub('^"|"$', "", header)

  blocks <- list()
  genes <- character(0)
  n_cells <- NA_integer_
  repeat {
    lines <- readLines(con, n = chunk_genes)
    if (length(lines) == 0L) break
    dt <- data.table::fread(text = lines, sep = sep, header = FALSE,
                            quote = "\"", showProgress = FALSE)
    if (is.na(n_cells)) {
      n_cells <- ncol(dt) - 1L
      if (length(header) == n_cells + 1L) header <- header[-1L]
      if (length(header) != n_cells)
        stop(sprintf("Header has %d fields but rows have %d values", length(header), n_cells))
    }
    genes <- c(genes, as.character(dt[[1L]]))
    m <- as.matrix(dt[, -1L])
    storage.mode(m) <- "double"
    blocks[[length(blocks) + 1L]] <- methods::as(m, "CsparseMatrix")
    rm(dt, m); invisible(gc(verbose = FALSE))
  }

  counts <- do.call(rbind, blocks)
  dimnames(counts) <- list(genes, header)

  if (anyDuplicated(genes)) {
    warning(sum(duplicated(genes)), " duplicated gene names made unique")
    rownames(counts) <- make.unique(genes)
  }
  if (any(counts@x != round(counts@x)))
    warning("Non-integer values found: is this really a raw count matrix?")
  counts
}

#' Read cell metadata and check the columns the pipeline relies on.
read_metadata <- function(path, cols) {
  md <- data.table::fread(path, data.table = FALSE)
  needed <- unlist(cols[c("cell", "sample", "patient", "site")])
  missing <- setdiff(needed, colnames(md))
  if (length(missing))
    stop("Metadata columns not found: ", paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(md), collapse = ", "),
         "\nEdit `cfg$meta_cols` in _targets.R.")
  md
}

#' Align counts and metadata and build the Seurat object.
#' Standardised columns sample / patient / site are added next to the originals.
build_seurat <- function(counts, metadata, cols) {
  cells <- as.character(metadata[[cols$cell]])
  shared <- intersect(colnames(counts), cells)
  if (length(shared) == 0L)
    stop("No cell identifiers shared between count matrix and metadata")

  md <- metadata[match(shared, cells), , drop = FALSE]
  rownames(md) <- shared
  md$sample  <- as.character(md[[cols$sample]])
  md$patient <- as.character(md[[cols$patient]])
  md$site    <- as.character(md[[cols$site]])

  SeuratObject::CreateSeuratObject(counts = counts[, shared], meta.data = md,
                                   project = "hcc", min.cells = 0, min.features = 0)
}

#' One-row-per-check summary of what was read and how well it matched.
ingest_summary <- function(counts, metadata, seu, cols) {
  cells_md <- as.character(metadata[[cols$cell]])
  data.frame(
    check = c("genes_in_matrix", "cells_in_matrix", "cells_in_metadata",
              "cells_matched", "cells_only_in_matrix", "cells_only_in_metadata",
              "nonzero_fraction", "patients", "samples", "sites"),
    value = c(nrow(counts), ncol(counts), length(cells_md), ncol(seu),
              length(setdiff(colnames(counts), cells_md)),
              length(setdiff(cells_md, colnames(counts))),
              signif(length(counts@x) / (as.numeric(nrow(counts)) * ncol(counts)), 3),
              length(unique(seu$patient)), length(unique(seu$sample)),
              paste(sort(unique(seu$site)), collapse = ", "))
  )
}

#' Cells per patient x site: shows which comparisons can be paired.
design_table <- function(seu) {
  tab <- as.data.frame.matrix(table(seu$patient, seu$site))
  cbind(patient = rownames(tab), tab, row.names = NULL)
}
