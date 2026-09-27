# Shared helpers -------------------------------------------------------------

# Storage format for large targets: qs2 (fast, compact), defined explicitly so
# it does not depend on which serialiser a given `targets` version maps "qs" to.
format_qs2 <- targets::tar_format(
  read  = function(path) qs2::qs_read(path),
  write = function(object, path) qs2::qs_save(object, path, nthreads = 2L)
)

# Write a data frame as TSV and return the path (for format = "file" targets).
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t")
  path
}

# Save a ggplot / patchwork object as PNG and return the path.
save_plot <- function(plot, path, width = 10, height = 6, dpi = 150) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi, bg = "white")
  path
}
