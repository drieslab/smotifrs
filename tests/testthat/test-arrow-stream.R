# The path-based reader sees a store's files as they sit on disk. A stream
# moves that decision to the producer, so a caller whose view differs from the
# files -- a pending subset, a multi-file dataset -- gets the edge set it
# meant. These tests pin that the two paths agree when the producer hands over
# everything, and that the stream is what decides when it does not.

skip_if_no_stream <- function() {
    skip_if_not_installed("arrow")
    skip_if_not_installed("nanoarrow")
}

.edges_tbl <- function(from, to) {
    arrow::arrow_table(from_id = as.integer(from), to_id = as.integer(to))
}

test_that("a stream gives the same answer as the same edges by path", {
    skip_if_no_stream()
    set.seed(1)
    n <- 120L
    g <- igraph::sample_gnp(n, 0.06)
    el <- igraph::as_edgelist(g, names = FALSE)
    ct <- factor(sample(c("A", "B", "C"), n, TRUE))

    d <- withr::local_tempdir()
    dir.create(file.path(d, "nodes")); dir.create(file.path(d, "edges"))
    np <- file.path(d, "nodes", "nodes.parquet")
    ep <- file.path(d, "edges", "edges.parquet")
    arrow::write_parquet(
        data.frame(node_id = as.character(seq_len(n)), int_id = seq_len(n)), np)
    arrow::write_parquet(
        data.frame(from_id = as.integer(el[, 1]), to_id = as.integer(el[, 2])), ep)

    by_path <- motif_enrichment_edge_store(np, ep, ct, size = 3L,
                                           n_perm = 49L, seed = 7L)
    by_stream <- motif_enrichment_stream(.edges_tbl(el[, 1], el[, 2]),
                                         int_ids = seq_len(n), cell_type = ct,
                                         size = 3L, n_perm = 49L, seed = 7L)
    data.table::setorder(by_path, motif_id)
    data.table::setorder(by_stream, motif_id)
    expect_identical(by_path$motif_id, by_stream$motif_id)
    expect_equal(by_path$observed, by_stream$observed)
    expect_equal(by_path$p_enrich, by_stream$p_enrich)
})

test_that("the stream decides the edge set, not the file", {
    skip_if_no_stream()
    set.seed(2)
    n <- 100L
    g <- igraph::sample_gnp(n, 0.08)
    el <- igraph::as_edgelist(g, names = FALSE)
    ct <- factor(sample(c("A", "B"), n, TRUE))

    full <- motif_enrichment_stream(.edges_tbl(el[, 1], el[, 2]),
                                    int_ids = seq_len(n), cell_type = ct,
                                    size = 3L, n_perm = 19L, seed = 1L)
    # the narrowing a producer would have applied before yielding batches
    keep <- el[, 1] <= 50L & el[, 2] <= 50L
    part <- motif_enrichment_stream(.edges_tbl(el[keep, 1], el[keep, 2]),
                                    int_ids = seq_len(n), cell_type = ct,
                                    size = 3L, n_perm = 19L, seed = 1L)
    expect_lt(sum(part$observed), sum(full$observed))
})

test_that("a multi-batch stream is pulled to exhaustion", {
    skip_if_no_stream()
    set.seed(3)
    n <- 80L
    el <- igraph::as_edgelist(igraph::sample_gnp(n, 0.1), names = FALSE)
    ct <- factor(sample(c("A", "B"), n, TRUE))
    one <- .edges_tbl(el[, 1], el[, 2])

    half <- nrow(el) %/% 2L
    chunked <- arrow::RecordBatchReader$create(
        arrow::record_batch(from_id = as.integer(el[seq_len(half), 1]),
                            to_id = as.integer(el[seq_len(half), 2])),
        arrow::record_batch(from_id = as.integer(el[-seq_len(half), 1]),
                            to_id = as.integer(el[-seq_len(half), 2]))
    )
    a <- motif_enrichment_stream(one, seq_len(n), ct, size = 3L, n_perm = 19L, seed = 5L)
    b <- motif_enrichment_stream(chunked, seq_len(n), ct, size = 3L, n_perm = 19L, seed = 5L)
    data.table::setorder(a, motif_id); data.table::setorder(b, motif_id)
    expect_identical(a$motif_id, b$motif_id)
    expect_equal(a$observed, b$observed)
})

test_that("a gappy node universe is remapped, not rejected", {
    skip_if_no_stream()
    ints <- c(3L, 9L, 14L, 20L)          # deliberately non-contiguous
    ct <- factor(c("A", "B", "A", "B"))
    e <- .edges_tbl(c(3L, 9L, 3L), c(9L, 14L, 14L))
    r <- motif_enrichment_stream(e, ints, ct, size = 3L, n_perm = 9L, seed = 1L)
    expect_gt(nrow(r), 0L)
})

test_that("an endpoint outside the node universe is an error", {
    skip_if_no_stream()
    expect_error(
        motif_enrichment_stream(.edges_tbl(c(1L, 2L), c(2L, 99L)),
                                int_ids = 1:3, cell_type = factor(c("A","B","A")),
                                size = 3L, n_perm = 9L),
        "not in the node sidecar"
    )
})

test_that("cell_type must match the node universe", {
    skip_if_no_stream()
    expect_error(
        motif_enrichment_stream(.edges_tbl(1L, 2L), int_ids = 1:3,
                                cell_type = factor(c("A", "B"))),
        "cell_type has 2 entries but int_ids has 3"
    )
})
