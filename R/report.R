# HTML report --------------------------------------------------------------------------
#
# Rendered as the last step of the pipeline into docs/index.html (served by GitHub
# Pages from the docs/ folder). The report only reads files written by earlier
# targets; their paths are passed as `inputs`, so it is rebuilt whenever any
# result table or figure changes.

render_report <- function(rmd, inputs, output_dir = "docs") {
  root <- normalizePath(".")
  out_dir <- file.path(root, output_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rmarkdown::render(rmd, output_file = "index.html", output_dir = out_dir,
                    knit_root_dir = root, intermediates_dir = tempdir(),
                    params = list(root = root), envir = new.env(), quiet = TRUE)
  file.path(output_dir, "index.html")
}
