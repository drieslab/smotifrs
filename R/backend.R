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


#' Motif enrichment in one pass
#'
#' Enumerates every connected induced subgraph of the given size, then tests
#' each colored motif class against a label-permutation null. Enumeration
#' happens once and each permutation is a relabel-and-count over the fixed
#' instance set, so per-instance data never reaches R and the cost of the null
#' is integer arithmetic rather than string construction.
#'
#' The graph is supplied as plain vectors rather than as any package's graph
#' class, so callers outside \pkg{smotif} can use it directly. Cell types are
#' passed as a factor (or coerced to one), which fixes the level ordering on the
#' R side -- the backend never applies its own collation.
#'
#' @param from,to integer vectors of 1-based node indices, one pair per edge.
#'   Direction is ignored and duplicate or reversed edges are collapsed.
#' @param cell_type factor or character vector of cell type labels, one per
#'   node, in node-index order.
#' @param n_nodes number of nodes. Defaults to `length(cell_type)`.
#' @param size motif size: 2, 3 or 4.
#' @param n_perm number of label permutations.
#' @param seed integer seed. Results depend only on `(seed, draw index)`, never
#'   on thread scheduling.
#' @param null which null to test against.
#'   * `"label"` shuffles labels over all nodes. Answers "is this motif more
#'     common than the tissue's cell type composition implies".
#'   * `"stratified"` shuffles only within `strata`, preserving each stratum's
#'     composition exactly.
#'   * `"conditional"` holds the observed **pairwise** edge-type composition
#'     and randomizes everything else, so a size-3 or size-4 motif that is
#'     still enriched is enriched *beyond* what its constituent pairs already
#'     explain. Under `"label"` any motif built from an attracting pair looks
#'     enriched, and the higher-order signal cannot be separated from an echo
#'     of the pairwise one.
#'
#'   The conditional null conditions by tolerance, not exactly: a Metropolis
#'   chain over label swaps holds the pairwise table close rather than pinning
#'   it. Three attributes report what the chain actually did -- `cond_dev` (the
#'   table deviation, as a fraction of all edges), `cond_accept` (acceptance
#'   rate) and `cond_moved` (mean fraction of labels displaced from the observed
#'   assignment). Report them rather than assuming the constraint held.
#'
#'   Watch `cond_moved` in particular. A chain held too cold barely moves, so
#'   every draw is essentially the observed data and *everything* comes back
#'   non-significant -- a failure that looks like a clean result. Below 5%
#'   displacement this warns.
#' @param cond_temp Metropolis temperature for `null = "conditional"`. Lower
#'   holds the pairwise table tighter but mixes more slowly.
#' @param strata optional factor or integer vector, one per node. Required for
#'   `null = "stratified"`. When given, labels are only exchanged between nodes
#'   sharing a stratum, which preserves each stratum's cell type composition
#'   exactly.
#' @param anchored_on optional integer vector of 1-based node indices.
#'   Enumeration is restricted to subgraphs containing at least one of them.
#' @returns a `data.table`, one row per motif class, with columns `motif_id`,
#'   `topology`, `size`, `color_tuple` (a list column of labels in canonical
#'   orbit order), `observed`, `expected`, `sd_null`, `z`, `fold`, `p_enrich`,
#'   `p_deplete` and `p_adj`. `n_instances` is attached as an attribute.
#'
#'   Positions within `color_tuple` are structural roles, not an arbitrary
#'   sort, so `A-B-B` and `B-A-A` are different motifs. The first position is
#'   the distinguished vertex wherever a topology has one: the centre of a
#'   `open` wedge or a `claw`, and the degree-3 vertex of a `paw`. For `path`
#'   the tuple runs end to end, for `cycle` it runs around the ring, and for
#'   `closed` (triangle), `diamond` and `K4` the remaining positions are
#'   interchangeable within their degree class.
#' @examples
#' # a 6-cycle alternating between two cell types
#' from <- c(1L, 2L, 3L, 4L, 5L, 6L)
#' to <- c(2L, 3L, 4L, 5L, 6L, 1L)
#' ct <- rep(c("A", "B"), 3)
#' motif_enrichment_rs(from, to, ct, size = 3L, n_perm = 99L)
#' @export
motif_enrichment_rs <- function(from,
                                to,
                                cell_type,
                                n_nodes = NULL,
                                size = 3L,
                                n_perm = 1000L,
                                seed = 1L,
                                null = c("label", "stratified", "conditional"),
                                cond_temp = 1,
                                strata = NULL,
                                anchored_on = NULL) {
  null <- match.arg(null)
  size <- as.integer(size)
  n_perm <- as.integer(n_perm)
  if (length(size) != 1L || is.na(size) || !size %in% 2:4) {
    stop("size must be 2, 3 or 4", call. = FALSE)
  }
  if (length(n_perm) != 1L || is.na(n_perm) || n_perm < 1L) {
    stop("n_perm must be a positive integer", call. = FALSE)
  }
  if (length(from) != length(to)) {
    stop("from and to must have equal length", call. = FALSE)
  }

  ct <- if (is.factor(cell_type)) cell_type else factor(cell_type)
  if (anyNA(ct)) stop("cell_type may not contain NA", call. = FALSE)
  n <- if (is.null(n_nodes)) length(ct) else as.integer(n_nodes)
  if (length(ct) != n) {
    stop("cell_type must have one entry per node", call. = FALSE)
  }
  # With n_nodes left NULL the node count comes from cell_type, so an edge
  # naming a node beyond it is a length mismatch the caller cannot see. Catch
  # it here rather than letting it surface as an out-of-range error from Rust.
  if (length(from)) {
    hi <- max(max(as.integer(from)), max(as.integer(to)))
    if (is.finite(hi) && hi > n) {
      stop(sprintf(
        "edges reference node %d but cell_type covers only %d node%s",
        hi, n, if (n == 1L) "" else "s"
      ), call. = FALSE)
    }
    lo <- min(min(as.integer(from)), min(as.integer(to)))
    if (is.finite(lo) && lo < 1L) {
      stop("from and to must be 1-based node indices", call. = FALSE)
    }
  }

  if (identical(null, "stratified") && is.null(strata)) {
    stop('null = "stratified" needs a strata vector', call. = FALSE)
  }
  if (!is.null(strata) && identical(null, "label")) {
    null <- "stratified"
  }
  strata_i <- NULL
  if (!is.null(strata)) {
    strata_i <- as.integer(if (is.factor(strata)) strata else factor(strata))
    if (length(strata_i) != n) {
      stop("strata must have one entry per node", call. = FALSE)
    }
  }
  anchors_i <- if (is.null(anchored_on)) NULL else as.integer(anchored_on)

  rs <- rs_motif_enrichment(
    n_nodes = n,
    type_codes = as.integer(ct),
    type_levels = levels(ct),
    edge_source = as.integer(from),
    edge_target = as.integer(to),
    size = size,
    n_perm = n_perm,
    seed = as.integer(seed),
    strata = strata_i,
    anchored_on = anchors_i,
    null_kind = null,
    cond_temp = as.numeric(cond_temp)
  )

  nclass <- length(rs$motif_id)
  cols <- matrix(rs$color_flat, nrow = nclass, ncol = rs$k, byrow = TRUE)
  eps <- .Machine$double.eps

  out <- data.table::data.table(
    motif_id = rs$motif_id,
    topology = rs$topology,
    size = as.integer(rs$size),
    color_tuple = split(cols, row(cols)),
    observed = rs$observed,
    expected = rs$mean_null,
    sd_null = rs$sd_null,
    p_enrich = rs$p_enrich,
    p_deplete = rs$p_deplete
  )
  names(out$color_tuple) <- NULL
  out[, "z" := (observed - expected) / pmax(sd_null, sqrt(eps))]
  out[, "fold" := observed / pmax(expected, eps)]
  # two-sided p from the more extreme tail, then BH across all classes
  out[, "p_adj" := stats::p.adjust(
    pmin(1, 2 * pmin(p_enrich, p_deplete)),
    method = "BH"
  )]
  data.table::setcolorder(out, c(
    "motif_id", "topology", "size", "color_tuple", "observed", "expected",
    "sd_null", "z", "fold", "p_enrich", "p_deplete", "p_adj"
  ))
  data.table::setorder(out, p_adj, -z)
  attr(out, "n_instances") <- rs$n_instances
  attr(out, "n_perm") <- rs$n_perm
  attr(out, "null") <- null
  if (identical(null, "conditional")) {
    # how tightly the constraint held, as a fraction of all edges
    attr(out, "cond_dev") <- if (rs$n_edges > 0) rs$cond_dev / rs$n_edges else NA_real_
    attr(out, "cond_accept") <- rs$cond_accept
    attr(out, "cond_moved") <- rs$cond_moved
    # A chain held too cold barely moves, and then every draw is the observed
    # data. That reads as "nothing is significant" rather than as a failure, so
    # it has to be said out loud rather than left in an attribute.
    attr(out, "cond_moves_per_draw") <- rs$cond_moves_per_draw
    # Draws are only decorrelated once the chain makes on the order of one
    # accepted swap per node between samples. Below that, consecutive draws are
    # near-copies of each other and of the observed data -- which surfaces as
    # "nothing is significant" rather than as a failure, so say it out loud.
    # Two independent signals, and both must fire. Displacement alone is not
    # enough: when the pairwise table genuinely pins most labels, low
    # displacement is the correct answer rather than a mixing failure.
    if (is.finite(rs$cond_moves_per_draw) &&
        rs$cond_moves_per_draw < rs$n_nodes_out / 4 &&
        is.finite(rs$cond_moved) && rs$cond_moved < 0.25) {
      warning(sprintf(
        paste0(
          "the conditional null mixed poorly: %.0f accepted swaps between ",
          "draws for %d nodes (acceptance %.1f%%, %.1f%% of labels displaced). ",
          "Consecutive draws are near-copies, so p-values are unreliable and ",
          "will tend to look uniformly non-significant. Raise cond_temp above ",
          "%g, or raise n_perm."
        ),
        rs$cond_moves_per_draw, rs$n_nodes_out, 100 * rs$cond_accept,
        100 * rs$cond_moved, cond_temp
      ), call. = FALSE)
    }
  } else {
    attr(out, "cond_dev") <- NA_real_
    attr(out, "cond_accept") <- NA_real_
    attr(out, "cond_moved") <- NA_real_
  }
  out[]
}


#' Instances of selected motif classes
#'
#' Returns the cells making up each occurrence of the named motif classes.
#'
#' Deliberately narrow. [motif_enrichment_rs()] never returns per-instance data
#' because at size 4 on a real dataset that is tens of millions of rows; this
#' returns instances for *named classes only*, so the output is bounded by those
#' classes' observed counts. Pick the motifs worth looking at first, then ask
#' for their instances.
#'
#' @inheritParams motif_enrichment_rs
#' @param motif_ids character vector of `motif_id` values from
#'   [motif_enrichment_rs()].
#' @param max_per_class cap on instances returned per class. `Inf` for no cap.
#' @returns a `data.table` with `motif_id`, `instance`, `slot` (the structural
#'   position within the motif, 1-based) and `node` (the 1-based node index),
#'   in long form -- `k` rows per instance.
#' @examples
#' from <- c(1L, 2L, 3L, 4L, 5L, 6L)
#' to <- c(2L, 3L, 4L, 5L, 6L, 1L)
#' ct <- rep(c("A", "B"), 3)
#' e <- motif_enrichment_rs(from, to, ct, size = 3L, n_perm = 19L)
#' motif_instances_rs(from, to, ct, size = 3L, motif_ids = e$motif_id[1])
#' @export
motif_instances_rs <- function(from,
                               to,
                               cell_type,
                               motif_ids,
                               n_nodes = NULL,
                               size = 3L,
                               max_per_class = 5000) {
  size <- as.integer(size)
  if (length(size) != 1L || is.na(size) || !size %in% 2:4) {
    stop("size must be 2, 3 or 4", call. = FALSE)
  }
  motif_ids <- unique(as.character(motif_ids))
  if (!length(motif_ids)) {
    stop("motif_ids must name at least one motif class", call. = FALSE)
  }
  ct <- if (is.factor(cell_type)) cell_type else factor(cell_type)
  n <- if (is.null(n_nodes)) length(ct) else as.integer(n_nodes)
  if (length(ct) != n) {
    stop("cell_type must have one entry per node", call. = FALSE)
  }
  cap <- if (is.finite(max_per_class)) as.integer(max_per_class) else 0L

  rs <- rs_motif_instances(
    n_nodes = n,
    type_codes = as.integer(ct),
    type_levels = levels(ct),
    edge_source = as.integer(from),
    edge_target = as.integer(to),
    size = size,
    motif_ids = motif_ids,
    max_per_class = cap
  )
  k <- rs$k
  n_inst <- length(rs$which)
  if (n_inst == 0L) {
    return(data.table::data.table(
      motif_id = character(), instance = integer(),
      slot = integer(), node = integer()
    ))
  }
  data.table::data.table(
    motif_id = rep(motif_ids[rs$which], each = k),
    instance = rep(seq_len(n_inst), each = k),
    slot = rep(seq_len(k), times = n_inst),
    node = as.integer(rs$verts)
  )
}
