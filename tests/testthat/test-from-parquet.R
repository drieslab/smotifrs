test_that("find_motifs_from_parquet matches the SpatialGraph round-trip", {
  skip_if_no_smotif()
  testthat::skip_if_not_installed("arrow")

  sg <- fixture_random(N = 40L, seed = 4L)
  td <- withr::local_tempdir()
  np <- file.path(td, "nodes.parquet")
  ep <- file.path(td, "edges.parquet")
  smotif::write_spatial_graph(sg, np, ep)

  for (size in c(2L, 3L, 4L)) {
    m_via_sg <- smotif::find_motifs(sg, size = size, backend = "rust")
    m_direct <- smotifrs::find_motifs_from_parquet(np, ep, size = size)
    expect_equal(canon_catalog(m_via_sg$catalog),
                 canon_catalog(m_direct$catalog),
                 info = sprintf("size=%d", size))
  }
})

test_that("find_motifs_from_parquet honors anchored_on", {
  skip_if_no_smotif()
  testthat::skip_if_not_installed("arrow")
  sg <- fixture_random(N = 30L, seed = 8L)
  td <- withr::local_tempdir()
  np <- file.path(td, "nodes.parquet")
  ep <- file.path(td, "edges.parquet")
  smotif::write_spatial_graph(sg, np, ep)

  anchors <- sg$nodes$cell_id[1:3]
  m <- smotifrs::find_motifs_from_parquet(np, ep, size = 3L,
                                          anchored_on = anchors,
                                          colored = FALSE)
  data.table::setDT(m$incidence)
  by_inst <- m$incidence[, .(touches = any(cell_id %in% anchors)),
                         by = "instance_id"]
  expect_true(all(by_inst$touches))
})

test_that("find_motifs_from_parquet errors clearly on bad paths", {
  expect_error(
    smotifrs::find_motifs_from_parquet("/no/such/nodes.parquet",
                                       "/no/such/edges.parquet",
                                       size = 2L),
    "nodes_path not found"
  )
})
