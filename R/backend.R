#' @keywords internal
#' @noRd
.rs_list_to_motifs <- function(rs_out, size, colored, sg_hash = NA_character_,
                               anchored = FALSE, source = "rust") {
  if (!requireNamespace("smotif", quietly = TRUE)) {
    stop("smotif must be installed to construct a Motifs object",
         call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table must be installed", call. = FALSE)
  }
  catalog <- data.table::as.data.table(rs_out$catalog)
  incidence <- data.table::as.data.table(rs_out$incidence)

  instance_meta <- lapply(rs_out$instance_meta, data.table::as.data.table)

  smotif::Motifs(
    catalog       = catalog,
    incidence     = incidence,
    instance_meta = instance_meta,
    meta = list(
      size         = as.integer(size),
      backend      = source,
      colored      = isTRUE(colored),
      anchored     = isTRUE(anchored),
      sg_hash      = sg_hash,
      generated_at = Sys.time()
    )
  )
}

#' Find motifs in a SpatialGraph using the Rust backend.
#'
#' Drop-in replacement for `smotif::find_motifs()` reached via
#' `smotif::find_motifs(..., backend = "rust")`. Uses a CSR-backed Rust
#' kernel with sorted-merge triangle enumeration and rayon-friendly
#' loops over vertices. Output is bit-for-bit compatible with the
#' igraph backend so downstream `Motifs`-consumers (`null_models`,
#' `genes`, `compare`) work unchanged.
#'
#' @param sg a `smotif::SpatialGraph` object.
#' @param size integer 2, 3, or 4.
#' @param colored logical; if `TRUE`, distinguish motifs by color tuple.
#' @param anchored_on optional character vector of cell ids to restrict
#'   enumeration to motifs touching at least one anchor.
#' @param max_instances optional integer cap on enumerated instances.
#' @return a `smotif::Motifs` S3 object.
#' @export
find_motifs_rs <- function(sg,
                           size = 3L,
                           colored = TRUE,
                           anchored_on = NULL,
                           max_instances = NULL) {
  if (!inherits(sg, "SpatialGraph")) {
    stop("`sg` must be a smotif::SpatialGraph", call. = FALSE)
  }
  size <- as.integer(size)
  if (length(size) != 1L || !size %in% c(2L, 3L, 4L)) {
    stop("`size` must be one of 2, 3, 4", call. = FALSE)
  }

  cell_ids   <- as.character(sg$nodes$cell_id)
  cell_types <- as.character(sg$nodes$cell_type)
  src <- match(as.character(sg$edges$source), cell_ids)
  tgt <- match(as.character(sg$edges$target), cell_ids)
  if (anyNA(src) || anyNA(tgt)) {
    stop("some edge endpoints are not in sg$nodes$cell_id", call. = FALSE)
  }

  anchors_ix <- if (is.null(anchored_on)) NULL else {
    miss <- setdiff(anchored_on, cell_ids)
    if (length(miss)) {
      stop(length(miss),
           " cell_id(s) in anchored_on are not in sg$nodes$cell_id",
           call. = FALSE)
    }
    as.integer(match(anchored_on, cell_ids))
  }

  rs_out <- rs_find_motifs(
    cell_ids = cell_ids, cell_types = cell_types,
    edge_source = as.integer(src), edge_target = as.integer(tgt),
    size = size, colored = isTRUE(colored),
    anchored_on = anchors_ix
  )

  m <- .rs_list_to_motifs(rs_out, size = size, colored = colored,
                          anchored = !is.null(anchored_on))

  if (!is.null(max_instances)) {
    m <- .cap_instances_smotif(m, max_instances)
  }
  m
}

#' Find motifs by reading nodes/edges parquet directly in Rust.
#'
#' Skips the R-side `SpatialGraph` round-trip: opens both parquet files
#' inside Rust via `arrow-rs`, builds a CSR adjacency, and returns the
#' same `Motifs` object as [find_motifs_rs()]. Intended for million-cell
#' datasets where holding the graph in R memory would be wasteful.
#'
#' @param nodes_path,edges_path file paths to the parquet outputs of
#'   `smotif::write_spatial_graph()`.
#' @inheritParams find_motifs_rs
#' @return a `smotif::Motifs` S3 object.
#' @export
find_motifs_from_parquet <- function(nodes_path,
                                     edges_path,
                                     size = 3L,
                                     colored = TRUE,
                                     anchored_on = NULL,
                                     max_instances = NULL) {
  if (!file.exists(nodes_path)) {
    stop("nodes_path not found: ", nodes_path, call. = FALSE)
  }
  if (!file.exists(edges_path)) {
    stop("edges_path not found: ", edges_path, call. = FALSE)
  }
  size <- as.integer(size)
  if (length(size) != 1L || !size %in% c(2L, 3L, 4L)) {
    stop("`size` must be one of 2, 3, 4", call. = FALSE)
  }

  anchors_ix <- if (is.null(anchored_on)) NULL else {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("arrow is needed to resolve anchored_on against nodes.parquet",
           call. = FALSE)
    }
    nodes <- arrow::read_parquet(nodes_path)
    miss <- setdiff(anchored_on, nodes$cell_id)
    if (length(miss)) {
      stop(length(miss), " cell_id(s) in anchored_on not found in nodes.parquet",
           call. = FALSE)
    }
    as.integer(match(anchored_on, nodes$cell_id))
  }

  rs_out <- rs_find_motifs_from_parquet(
    nodes_path = normalizePath(nodes_path, mustWork = TRUE),
    edges_path = normalizePath(edges_path, mustWork = TRUE),
    size = size, colored = isTRUE(colored),
    anchored_on_idx = anchors_ix
  )
  m <- .rs_list_to_motifs(rs_out, size = size, colored = colored,
                          anchored = !is.null(anchored_on),
                          source = "rust_parquet")

  if (!is.null(max_instances)) {
    m <- .cap_instances_smotif(m, max_instances)
  }
  m
}

# Honor max_instances by random-sampling. Mirrors the smotif-side cap.
# We re-implement here rather than calling internal smotif helpers so
# smotifrs stays decoupled from smotif's private API.
.cap_instances_smotif <- function(m, max_instances) {
  total <- if (nrow(m$catalog)) sum(m$catalog$count) else 0L
  if (total <= max_instances) return(m)
  warning(sprintf(
    "max_instances = %d below total enumerated count (%d); subsampling",
    max_instances, total
  ), call. = FALSE)
  set.seed(1L)
  all_inst <- unique(m$incidence$instance_id)
  keep <- sort(sample(all_inst, max_instances, replace = FALSE))
  m$incidence <- m$incidence[m$incidence$instance_id %in% keep]
  m$catalog <- m$incidence[, .(count = data.table::uniqueN(instance_id)),
                           by = "motif_id"]
  for (cls in names(m$instance_meta)) {
    m$instance_meta[[cls]] <-
      m$instance_meta[[cls]][m$instance_meta[[cls]]$instance_id %in% keep]
  }
  m
}
