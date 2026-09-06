# The conditional null exists to separate a genuine multicellular niche from an
# echo of pairwise attraction. These are its two acceptance criteria: it must
# NOT call a purely pairwise pattern a higher-order motif, and it must still
# detect structure the pairwise table cannot see. A third test pins the mixing
# diagnostic, because a chain held too cold returns "nothing is significant"
# and that failure is indistinguishable from a clean result without it.

.ring <- function(n) list(from = seq_len(n), to = c(seq_len(n)[-1], 1L))

test_that("a purely pairwise pattern is not called a higher-order motif", {
    # Labels are drawn as a first-order Markov chain along the ring, strongly
    # favouring alternation. That is the one construction that is *exactly*
    # pairwise: x[i+1] depends on x[i] and, given x[i], on nothing earlier. So
    # the pairwise edge-type table carries all the information there is, and a
    # null that preserves it should find no higher-order signal.
    #
    # Two constructions that look similar are wrong here and were tried first.
    # A perfect alternation pins the assignment outright -- with A-A and B-B
    # both zero the chain cannot move, and "not significant" would be an
    # artefact of a frozen chain. Alternation with a fraction of labels flipped
    # is worse than that: it correlates positions i and i+2, which the pairwise
    # table does not constrain, so it contains genuine higher-order structure
    # and the conditional null is right to flag it.
    set.seed(11)
    n <- 500
    g <- .ring(n)
    ct <- character(n)
    ct[1] <- "A"
    for (i in 2:n) {
        ct[i] <- if (runif(1) < 0.8) {
            if (ct[i - 1] == "A") "B" else "A"
        } else {
            ct[i - 1]
        }
    }

    marg <- motif_enrichment_rs(g$from, g$to, ct,
        size = 4L, n_perm = 300L, seed = 4L, null = "label"
    )
    cond <- motif_enrichment_rs(g$from, g$to, ct,
        size = 4L, n_perm = 300L, seed = 4L, null = "conditional"
    )
    # the comparison is only meaningful if the chain actually explored
    expect_gt(attr(cond, "cond_moved"), 0.25)
    expect_lt(attr(cond, "cond_dev"), 0.2)

    alt <- "size4_path_A-B-A-B"
    zm <- marg$z[marg$motif_id == alt]
    zc <- cond$z[cond$motif_id == alt]
    expect_length(zm, 1L)
    expect_length(zc, 1L)
    # the marginal null is fooled by the pairwise attraction
    expect_gt(zm, 5)
    # holding that attraction fixed removes essentially all of the signal
    expect_lt(abs(zc), 3)
})

test_that("structure the pairwise table cannot see is still detected", {
    # Same pairwise composition, different placement: all heterotypic edges sit
    # inside triangles while the 3-paths stay pure. A null that only preserves
    # edge-type counts may move them into the paths, so the triangle counts must
    # still stand out.
    from <- integer(0)
    to <- integer(0)
    ct <- character(0)
    id <- 0L
    for (i in 1:80) {
        from <- c(from, id + 1L, id + 2L, id + 3L)
        to <- c(to, id + 2L, id + 3L, id + 1L)
        ct <- c(ct, "A", "B", "C")
        id <- id + 3L
    }
    for (i in 1:80) {
        lab <- c("A", "B", "C")[(i %% 3) + 1L]
        from <- c(from, id + 1L, id + 2L)
        to <- c(to, id + 2L, id + 3L)
        ct <- c(ct, lab, lab, lab)
        id <- id + 3L
    }
    r <- motif_enrichment_rs(from, to, ct,
        size = 3L, n_perm = 200L, seed = 2L, null = "conditional"
    )
    tri <- r[r$motif_id == "size3_closed_A-B-C", ]
    expect_equal(nrow(tri), 1L)
    expect_gt(tri$z, 5)
    expect_lt(tri$p_adj, 0.05)
})

test_that("a chain that fails to mix is reported, not silently trusted", {
    from <- integer(0)
    to <- integer(0)
    ct <- character(0)
    id <- 0L
    for (i in 1:60) {
        from <- c(from, id + 1L, id + 2L, id + 3L)
        to <- c(to, id + 2L, id + 3L, id + 1L)
        ct <- c(ct, "A", "B", "C")
        id <- id + 3L
    }
    for (i in 1:60) {
        lab <- c("A", "B", "C")[(i %% 3) + 1L]
        from <- c(from, id + 1L, id + 2L)
        to <- c(to, id + 2L, id + 3L)
        ct <- c(ct, lab, lab, lab)
        id <- id + 3L
    }
    # far too cold: the chain barely accepts anything
    expect_warning(
        motif_enrichment_rs(from, to, ct,
            size = 3L, n_perm = 100L, seed = 2L,
            null = "conditional", cond_temp = 0.2
        ),
        "mixed poorly"
    )
    # the default temperature is fine on the same data
    expect_no_warning(
        motif_enrichment_rs(from, to, ct,
            size = 3L, n_perm = 100L, seed = 2L, null = "conditional"
        )
    )
})

test_that("conditional runs report their diagnostics and marginal runs do not", {
    g <- .ring(120)
    ct <- rep(c("A", "B", "C"), 40)
    cond <- motif_enrichment_rs(g$from, g$to, ct,
        size = 3L, n_perm = 60L, seed = 1L, null = "conditional"
    )
    for (a in c("cond_dev", "cond_accept", "cond_moved", "cond_moves_per_draw")) {
        expect_true(is.finite(attr(cond, a)), info = a)
    }
    expect_identical(attr(cond, "null"), "conditional")

    marg <- motif_enrichment_rs(g$from, g$to, ct,
        size = 3L, n_perm = 60L, seed = 1L, null = "label"
    )
    expect_true(is.na(attr(marg, "cond_dev")))
    expect_identical(attr(marg, "null"), "label")
})

test_that("the stratified null demands strata", {
    g <- .ring(10)
    expect_error(
        motif_enrichment_rs(g$from, g$to, rep(c("A", "B"), 5),
            size = 3L, null = "stratified"
        ),
        "needs a strata vector"
    )
})
