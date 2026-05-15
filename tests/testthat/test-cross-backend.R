test_that("rust and igraph backends produce identical catalogs at sizes 2/3/4", {
  skip_if_no_smotif()
  for (seed in c(1L, 7L, 23L)) {
    sg <- fixture_random(N = 40L, seed = seed)
    for (size in c(2L, 3L, 4L)) {
      m_ig <- smotif::find_motifs(sg, size = size, backend = "igraph")
      m_rs <- smotif::find_motifs(sg, size = size, backend = "rust")
      expect_equal(canon_catalog(m_ig$catalog),
                   canon_catalog(m_rs$catalog),
                   info = sprintf("seed=%d size=%d", seed, size))
    }
  }
})

test_that("incidence (motif_id, cell_id, orbit_id) frequencies match across backends", {
  skip_if_no_smotif()
  sg <- fixture_random(N = 30L, seed = 5L)
  for (size in c(2L, 3L, 4L)) {
    m_ig <- smotif::find_motifs(sg, size = size, backend = "igraph")
    m_rs <- smotif::find_motifs(sg, size = size, backend = "rust")
    expect_equal(canon_incidence(m_ig$incidence),
                 canon_incidence(m_rs$incidence),
                 info = sprintf("size=%d", size))
  }
})

test_that("uncolored mode collapses motifs to topology only", {
  skip_if_no_smotif()
  sg <- fixture_random(N = 40L, seed = 11L)
  m <- smotif::find_motifs(sg, size = 3L, colored = FALSE, backend = "rust")
  expect_true(all(is.na(m$catalog$color_tuple)))
  expect_setequal(m$catalog$canonical_iso_class, c("triangle", "wedge"))
  # total = total connected 3-subgraphs
  m_ig <- smotif::find_motifs(sg, size = 3L, colored = FALSE, backend = "igraph")
  expect_equal(sum(m$catalog$count), sum(m_ig$catalog$count))
})

test_that("anchored_on restricts enumeration in both backends and matches", {
  skip_if_no_smotif()
  sg <- fixture_random(N = 50L, seed = 19L)
  anchors <- sg$nodes$cell_id[1:5]
  m_ig <- smotif::find_motifs(sg, size = 3L, anchored_on = anchors,
                              backend = "igraph", colored = FALSE)
  m_rs <- smotif::find_motifs(sg, size = 3L, anchored_on = anchors,
                              backend = "rust", colored = FALSE)
  expect_equal(canon_catalog(m_ig$catalog), canon_catalog(m_rs$catalog))
  # every instance must touch an anchor
  data.table::setDT(m_rs$incidence)
  by_inst <- m_rs$incidence[, .(touches = any(cell_id %in% anchors)),
                            by = "instance_id"]
  expect_true(all(by_inst$touches))
})

test_that("instance_meta tables exist for the right classes", {
  skip_if_no_smotif()
  sg <- fixture_random(N = 30L, seed = 3L)
  m2 <- smotif::find_motifs(sg, size = 2L, backend = "rust")
  m3 <- smotif::find_motifs(sg, size = 3L, backend = "rust")
  m4 <- smotif::find_motifs(sg, size = 4L, backend = "rust")
  expect_true("edge" %in% names(m2$instance_meta))
  expect_setequal(intersect(names(m3$instance_meta),
                            c("triangle", "wedge")),
                  intersect(c("triangle", "wedge"),
                            unique(m3$catalog$canonical_iso_class)))
  # at least one size-4 class should appear
  expect_true(any(c("path", "claw", "cycle", "paw", "diamond", "K4")
                  %in% names(m4$instance_meta)))
})

test_that("Motifs object meta records the rust backend", {
  skip_if_no_smotif()
  sg <- fixture_random(N = 20L)
  m <- smotif::find_motifs(sg, size = 3L, backend = "rust")
  expect_equal(m$meta$backend, "rust")
  expect_equal(m$meta$size, 3L)
  expect_true(isTRUE(m$meta$colored))
})
