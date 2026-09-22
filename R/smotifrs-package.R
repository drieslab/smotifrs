#' smotifrs: Rust backend for smotif motif enumeration
#'
#' Drop-in Rust acceleration for the [smotif] R package: parallel
#' enumeration of size-2 / size-3 / size-4 motifs in spatial omics
#' graphs, plus direct parquet I/O for million-cell datasets.
#'
#' @keywords internal
#' @importFrom data.table := .N .SD as.data.table data.table is.data.table setDT
"_PACKAGE"

# Suppress R CMD check NOTEs for data.table column references used via NSE.
utils::globalVariables(c(
  "cell_id", "instance_id", "motif_id", "orbit_id",
  "observed", "expected", "sd_null", "p_enrich", "p_deplete", "p_adj", "z",
  ".N", ".SD", "."
))
