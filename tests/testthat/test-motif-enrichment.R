# `motif_enrichment_rs()` enumerates once and folds every permutation into
# per-class statistics inside Rust. These tests pin the invariants that are
# cheap to state and expensive to lose: instance conservation, both-tail
# p-value bounds, that the exact canonical form keeps structurally different
# colorings apart, and that the level ordering comes from R rather than from
# Rust's byte collation.

.cycle <- function(n) list(from = seq_len(n), to = c(seq_len(n)[-1], 1L))

test_that("observed counts sum to the instance count", {
    cy <- .cycle(8)
    for (sz in 2:4) {
        r <- motif_enrichment_rs(cy$from, cy$to, rep(c("A", "B"), 4),
            size = sz, n_perm = 20L, seed = 1L
        )
        expect_equal(sum(r$observed), attr(r, "n_instances"))
    }
})

test_that("a permutation moves instances between classes but never creates them", {
    cy <- .cycle(10)
    r <- motif_enrichment_rs(cy$from, cy$to, rep(c("A", "B", "C", "A", "B"), 2),
        size = 3L, n_perm = 300L, seed = 3L
    )
    expect_equal(sum(r$expected), attr(r, "n_instances"), tolerance = 1e-9)
})

test_that("both tails are reported and neither can be exactly zero", {
    cy <- .cycle(9)
    r <- motif_enrichment_rs(cy$from, cy$to, rep(c("A", "B", "C"), 3),
        size = 3L, n_perm = 99L, seed = 5L
    )
    expect_true(all(r$p_enrich > 0 & r$p_enrich <= 1))
    expect_true(all(r$p_deplete > 0 & r$p_deplete <= 1))
    expect_gte(min(c(r$p_enrich, r$p_deplete)), 1 / 100)
    expect_true(all(r$p_adj >= 0 & r$p_adj <= 1))
})

test_that("a wedge's centre is distinguished from its ends", {
    # C6 labelled A B A B A B: three wedges centred on A (ends B,B) and three
    # centred on B (ends A,A). These are different motifs, not one class of six.
    cy <- .cycle(6)
    r <- motif_enrichment_rs(cy$from, cy$to, rep(c("A", "B"), 3),
        size = 3L, n_perm = 19L, seed = 1L
    )
    seen <- r[r$observed > 0, ]
    expect_setequal(seen$motif_id, c("size3_open_A-B-B", "size3_open_B-A-A"))
    expect_true(all(seen$observed == 3))
})

test_that("alternating and adjacent 4-cycle colourings are different classes", {
    # the collision the (degree, colour) sort produced: both reduce to degrees
    # (2,2,2,2) and sorted colours A-A-B-B, but they are not isomorphic
    cy <- .cycle(4)
    alt <- motif_enrichment_rs(cy$from, cy$to, c("A", "B", "A", "B"),
        size = 4L, n_perm = 9L, seed = 1L
    )
    adj <- motif_enrichment_rs(cy$from, cy$to, c("A", "A", "B", "B"),
        size = 4L, n_perm = 9L, seed = 1L
    )
    id_alt <- alt$motif_id[alt$observed > 0]
    id_adj <- adj$motif_id[adj$observed > 0]
    expect_length(id_alt, 1L)
    expect_length(id_adj, 1L)
    expect_false(identical(id_alt, id_adj))
})

test_that("a single cell type collapses to one class per realized topology", {
    cy <- .cycle(5)
    r <- motif_enrichment_rs(cy$from, cy$to, rep("A", 5),
        size = 3L, n_perm = 25L, seed = 1L
    )
    expect_equal(nrow(r), 1L)
    expect_equal(r$observed, 5)
    # with one label the null cannot move anything
    expect_equal(r$expected, 5, tolerance = 1e-9)
    expect_lt(r$sd_null, 1e-9)
})

test_that("the stratified null cannot exchange labels across strata", {
    # two disjoint triangles, each a pure single-label stratum
    from <- c(1L, 2L, 3L, 4L, 5L, 6L)
    to <- c(2L, 3L, 1L, 5L, 6L, 4L)
    r <- motif_enrichment_rs(from, to, c("A", "A", "A", "B", "B", "B"),
        size = 3L, n_perm = 50L, seed = 2L,
        strata = c(1, 1, 1, 2, 2, 2)
    )
    expect_equal(r$expected, r$observed, tolerance = 1e-9)
    expect_true(all(r$sd_null < 1e-9))
})

test_that("results depend on the seed but not on run order", {
    cy <- .cycle(12)
    ct <- rep(c("A", "B", "C"), 4)
    a <- motif_enrichment_rs(cy$from, cy$to, ct, size = 3L, n_perm = 60L, seed = 11L)
    b <- motif_enrichment_rs(cy$from, cy$to, ct, size = 3L, n_perm = 60L, seed = 11L)
    d <- motif_enrichment_rs(cy$from, cy$to, ct, size = 3L, n_perm = 60L, seed = 12L)
    expect_identical(a$expected, b$expected)
    expect_false(identical(a$expected, d$expected))
})

test_that("cell type level order is taken from R, not from Rust byte order", {
    # names whose locale collation differs from byte order; supplying a factor
    # with an explicit level order must be honoured verbatim
    cy <- .cycle(6)
    ct <- factor(rep(c("cd8_T", "CD8+ T cell"), 3),
        levels = c("cd8_T", "CD8+ T cell")
    )
    r <- motif_enrichment_rs(cy$from, cy$to, ct, size = 2L, n_perm = 9L, seed = 1L)
    expect_true(all(unlist(r$color_tuple) %in% levels(ct)))
    # reversing the declared level order must not change the counts, only labels
    ct2 <- factor(as.character(ct), levels = rev(levels(ct)))
    r2 <- motif_enrichment_rs(cy$from, cy$to, ct2, size = 2L, n_perm = 9L, seed = 1L)
    expect_equal(sum(r$observed), sum(r2$observed))
})

test_that("bad input is rejected with a clear message", {
    cy <- .cycle(4)
    expect_error(
        motif_enrichment_rs(cy$from, cy$to, rep("A", 4), size = 5L),
        "size must be 2, 3 or 4"
    )
    expect_error(
        motif_enrichment_rs(cy$from, cy$to, rep("A", 4), n_perm = 0L),
        "positive integer"
    )
    expect_error(
        motif_enrichment_rs(cy$from, cy$to[-1], rep("A", 4)),
        "equal length"
    )
    expect_error(
        motif_enrichment_rs(cy$from, cy$to, rep("A", 3)),
        "cell_type covers only 3 nodes"
    )
    expect_error(
        motif_enrichment_rs(c(0L, 1L), c(1L, 2L), rep("A", 3)),
        "1-based"
    )
})
