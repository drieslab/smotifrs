# GiottoDisk's parquetEdgeStore interns node ids at write time, so its edge
# parquet already holds integer endpoints. The reader here takes those columns
# straight through rather than re-hashing strings, which is work GiottoDisk
# already did. These tests pin the schema it depends on -- documented as locked
# on the GiottoDisk side -- and that the disk path and the in-memory path give
# the same answer, which is the only reason to have two.

skip_if_no_arrow <- function() skip_if_not_installed("arrow")

# write a store by hand, so the test does not need GiottoDisk installed
.write_store <- function(dir, from, to, node_ids, int_ids = seq_along(node_ids)) {
    dir.create(file.path(dir, "nodes"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(dir, "edges"), recursive = TRUE, showWarnings = FALSE)
    arrow::write_parquet(
        data.frame(
            row_index = seq_along(node_ids),
            node_id = node_ids,
            int_id = as.integer(int_ids)
        ),
        file.path(dir, "nodes", "nodes.parquet")
    )
    arrow::write_parquet(
        data.frame(from_id = as.integer(from), to_id = as.integer(to)),
        file.path(dir, "edges", "edges.parquet")
    )
    list(
        nodes = file.path(dir, "nodes", "nodes.parquet"),
        edges = file.path(dir, "edges", "edges.parquet")
    )
}

test_that("the sidecar is read back in its own order", {
    skip_if_no_arrow()
    d <- withr::local_tempdir()
    p <- .write_store(d, c(1L, 2L), c(2L, 3L), c("a", "b", "c"))
    sc <- edge_store_nodes(p$nodes)
    expect_identical(sc$node_id, c("a", "b", "c"))
    expect_identical(sc$int_id, 1:3)
})

test_that("the disk path and the in-memory path agree exactly", {
    skip_if_no_arrow()
    set.seed(4)
    n <- 200L
    ig <- igraph::sample_gnm(n, 700, directed = FALSE)
    el <- igraph::as_edgelist(ig, names = FALSE)
    ct <- sample(c("T", "B", "M"), n, TRUE)

    d <- withr::local_tempdir()
    p <- .write_store(d, el[, 1], el[, 2], sprintf("c%03d", seq_len(n)))

    for (sz in 2:4) {
        a <- motif_enrichment_edge_store(p$nodes, p$edges,
            cell_type = ct, size = sz, n_perm = 50L, seed = 7L
        )
        b <- motif_enrichment_rs(el[, 1], el[, 2], ct,
            size = sz, n_perm = 50L, seed = 7L
        )
        expect_identical(attr(a, "n_instances"), attr(b, "n_instances"))
        expect_setequal(a$motif_id, b$motif_id)
        m <- merge(
            a[, list(motif_id, x = observed)],
            b[, list(motif_id, y = observed)],
            by = "motif_id"
        )
        expect_equal(m$x, m$y, info = paste("size", sz))
        # same seed must give the same null, not just the same observed counts
        ma <- a[order(a$motif_id), ]
        mb <- b[order(b$motif_id), ]
        expect_equal(ma$expected, mb$expected, info = paste("size", sz))
    }
})

test_that("a non-contiguous int_id universe is remapped, not assumed", {
    # a subset store's sidecar is neither 0-based nor gap-free
    skip_if_no_arrow()
    d <- withr::local_tempdir()
    p <- .write_store(d,
        from = c(10L, 20L, 30L), to = c(20L, 30L, 10L),
        node_ids = c("x", "y", "z"), int_ids = c(10L, 20L, 30L)
    )
    r <- motif_enrichment_edge_store(p$nodes, p$edges,
        cell_type = c("A", "B", "C"), size = 3L, n_perm = 19L
    )
    expect_equal(sum(r$observed), 1) # one triangle
    expect_identical(r$topology[1], "closed")
})

test_that("int64 id columns are accepted", {
    skip_if_no_arrow()
    d <- withr::local_tempdir()
    dir.create(file.path(d, "nodes"), recursive = TRUE)
    dir.create(file.path(d, "edges"), recursive = TRUE)
    arrow::write_parquet(
        arrow::arrow_table(
            row_index = 1:3,
            node_id = c("a", "b", "c"),
            int_id = arrow::Array$create(1:3, type = arrow::int64())
        ),
        file.path(d, "nodes", "nodes.parquet")
    )
    arrow::write_parquet(
        arrow::arrow_table(
            from_id = arrow::Array$create(c(1L, 2L, 3L), type = arrow::int64()),
            to_id = arrow::Array$create(c(2L, 3L, 1L), type = arrow::int64())
        ),
        file.path(d, "edges", "edges.parquet")
    )
    r <- motif_enrichment_edge_store(
        file.path(d, "nodes", "nodes.parquet"),
        file.path(d, "edges", "edges.parquet"),
        cell_type = c("A", "B", "C"), size = 3L, n_perm = 19L
    )
    expect_equal(sum(r$observed), 1)
})

test_that("mismatched label length and missing files are reported", {
    skip_if_no_arrow()
    d <- withr::local_tempdir()
    p <- .write_store(d, c(1L, 2L), c(2L, 3L), c("a", "b", "c"))
    expect_error(
        motif_enrichment_edge_store(p$nodes, p$edges, cell_type = c("A", "B")),
        "colors length"
    )
    expect_error(
        motif_enrichment_edge_store("nope.parquet", p$edges, cell_type = "A"),
        "no such file"
    )
})

test_that("a real GiottoDisk store round-trips", {
    skip_if_not_installed("GiottoDisk")
    skip_if_no_arrow()
    set.seed(11)
    n <- 150L
    ig <- igraph::sample_gnm(n, 500, directed = FALSE)
    igraph::V(ig)$name <- sprintf("cell_%03d", seq_len(n))
    ct <- sample(c("T", "B"), n, TRUE)

    root <- file.path(withr::local_tempdir(), "store")
    st <- GiottoDisk::storeWrite(
        GiottoDisk::storeCreate(type = "parquetEdgeStore", path = root),
        ig,
        type = "spatial"
    )
    np <- list.files(file.path(root, "nodes"), "parquet$", full.names = TRUE)[1]
    ep <- list.files(file.path(root, "edges"), "parquet$", full.names = TRUE)[1]

    sc <- edge_store_nodes(np)
    lab <- stats::setNames(ct, igraph::V(ig)$name)[sc$node_id]
    a <- motif_enrichment_edge_store(np, ep,
        cell_type = lab, size = 3L, n_perm = 50L, seed = 3L
    )
    el <- igraph::as_edgelist(ig, names = FALSE)
    b <- motif_enrichment_rs(el[, 1], el[, 2], ct,
        size = 3L, n_perm = 50L, seed = 3L
    )
    expect_identical(attr(a, "n_instances"), attr(b, "n_instances"))
    expect_setequal(a$motif_id, b$motif_id)
})
