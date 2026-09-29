# Compartments -------------------------------------------------------------------
#
# Clusters from the first-pass annotation (Leiden, resolution 0.5, seed 1234) are
# grouped into compartments. Cluster numbers are only meaningful for this exact
# run, so check_compartments() verifies each compartment against the authors'
# labels and stops the pipeline if the mapping no longer fits the clusters.

add_compartments <- function(seu, compartment_map) {
  lut <- stats::setNames(rep(names(compartment_map), lengths(compartment_map)),
                         unlist(compartment_map))
  seu$compartment <- unname(lut[as.character(seu$cluster)])
  seu$compartment[is.na(seu$compartment)] <- "Unassigned"
  seu
}

check_compartments <- function(seu, author_col, expected, min_agreement = 0.8) {
  md <- seu[[]]
  out <- do.call(rbind, lapply(names(expected), function(comp) {
    d <- md[md$compartment == comp, ]
    agree <- mean(d[[author_col]] %in% expected[[comp]])
    data.frame(compartment = comp, n_cells = nrow(d),
               expected_author_labels = paste(expected[[comp]], collapse = "/"),
               agreement = round(agree, 3))
  }))
  bad <- out$compartment[out$n_cells > 0 & out$agreement < min_agreement]
  if (length(bad))
    stop("Compartment map no longer matches the clusters for: ", paste(bad, collapse = ", "),
         ". Re-check cluster numbers in the compartment_map target.")
  out
}
