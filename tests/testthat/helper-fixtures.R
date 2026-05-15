# Reusable fixtures for cross-backend parity tests. We keep these
# light so the suite stays fast under R CMD check.

skip_if_no_smotif <- function() {
  testthat::skip_if_not_installed("smotif")
}

# Random scattered cells in a unit square; kNN graph; mixed cell types.
fixture_random <- function(N = 50L, k = 4L, seed = 1L,
                           types = c("T", "B", "M")) {
  set.seed(seed)
  coords <- cbind(x = stats::runif(N), y = stats::runif(N))
  cell_types <- sample(types, N, replace = TRUE)
  cell_id <- paste0("c", seq_len(N))
  smotif::build_spatial_graph(
    coords, cell_types, "s1",
    cell_id = cell_id, method = "knn", k = k
  )
}

# Two-region fixture with cell-type composition shift; useful when
# checking that anchored_on / region-aware behavior matches.
fixture_two_regions <- function(N_per = 30L, seed = 42L) {
  set.seed(seed)
  build <- function(prefix, xc, prob) {
    n <- N_per
    data.frame(
      cell_id   = paste0(prefix, "_", seq_len(n)),
      x         = stats::runif(n, xc, xc + 1),
      y         = stats::runif(n, 0, 1),
      cell_type = sample(c("T", "B"), n, replace = TRUE, prob = prob),
      sample_id = "s1",
      region_id = prefix,
      stringsAsFactors = FALSE
    )
  }
  fx <- rbind(build("r1", 0, c(0.85, 0.15)),
              build("r2", 5, c(0.15, 0.85)))
  smotif::build_spatial_graph(
    fx[, c("x", "y")], fx$cell_type, fx$sample_id,
    cell_id = fx$cell_id, region_id = fx$region_id,
    method = "knn", k = 4L
  )
}

# Sort catalog/incidence the same way before comparing across backends.
canon_catalog <- function(cat) {
  out <- cat[order(cat$motif_id), c("motif_id", "size",
                                    "canonical_iso_class",
                                    "color_tuple", "count")]
  rownames(out) <- NULL
  out
}
canon_incidence <- function(inc) {
  # incidence rows aren't ordered by instance position in the R backend
  # vs Rust backend deterministically, so compare per-(motif_id, cell_id)
  # frequency rather than row-for-row.
  data.table::setDT(inc)
  inc[, .N, by = .(motif_id, cell_id, orbit_id)][order(motif_id, cell_id, orbit_id)]
}
