# Benchmark behind vignettes/rust-backend.Rmd section 4.
#
# Writes vignettes/bench/results.csv and vignettes/bench/meta.dcf, which the
# vignette reads instead of re-running this (~15 min, ~6 GB). Run from the
# package root against an installed smotifrs:
#
#   Rscript data-raw/bench-rust-backend.R
#
# CI runs it via .github/workflows/bench.yml (manual dispatch or release).
# BENCH_REPS sets repetitions per case; the median is recorded.

library(smotif)
library(smotifrs)

reps <- as.integer(Sys.getenv("BENCH_REPS", "3"))
out_dir <- file.path("vignettes", "bench")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_synthetic <- function(n, k = 6L, seed = 1L) {
  set.seed(seed)
  coords <- cbind(x = stats::runif(n, 0, 1000),
                  y = stats::runif(n, 0, 1000))
  types <- sample(c("T", "B", "M", "E", "X"), n, replace = TRUE)
  build_spatial_graph(
    coords, types, "synth",
    cell_id = paste0("c", seq_len(n)),
    method = "knn", k = k
  )
}

peak_mb_from_gc <- function(g) {
  # Robust across R versions: locate the "max used" column and read
  # the Mb counterpart immediately after it.
  i <- which(colnames(g) == "max used")
  if (!length(i) || ncol(g) < i + 1L) return(NA_real_)
  sum(g[, i + 1L])
}

bench_one <- function(fn) {
  invisible(gc(reset = TRUE, full = TRUE))
  t1 <- Sys.time()
  # Keep the result alive until t2 so its allocation is captured by gc().
  result_holder <- fn()
  t2 <- Sys.time()
  g <- gc(verbose = FALSE)
  rm(result_holder)
  c(elapsed   = as.numeric(difftime(t2, t1, units = "secs")),
    r_peak_mb = peak_mb_from_gc(g))
}

bench_median <- function(fn) {
  runs <- vapply(seq_len(reps), function(i) bench_one(fn), numeric(2))
  apply(runs, 1L, stats::median)
}

sizes <- c(small = 2500L, medium = 10000L, large = 30000L)
out <- list()
for (sz_name in names(sizes)) {
  n <- sizes[[sz_name]]
  message("building synthetic graph: ", n, " cells")
  sg_n <- make_synthetic(n)
  for (motif_size in c(3L, 4L)) {
    for (backend in c("igraph", "rust")) {
      label <- sprintf("%s/N=%d/size=%d/%s", sz_name, n, motif_size, backend)
      message("  ", label)
      m <- bench_median(
        function() find_motifs(sg_n, size = motif_size, backend = backend)
      )
      out[[length(out) + 1L]] <- data.frame(
        label = label, elapsed = m[["elapsed"]], r_peak_mb = m[["r_peak_mb"]],
        backend = backend, N = n, motif_size = motif_size, sz_name = sz_name
      )
    }
  }
}
bench <- do.call(rbind, out)
utils::write.csv(bench, file.path(out_dir, "results.csv"), row.names = FALSE)

cpu_model <- function() {
  if (file.exists("/proc/cpuinfo")) {
    l <- grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1L]
    return(trimws(sub("^model name\\s*:", "", l)))
  }
  if (Sys.info()[["sysname"]] == "Darwin") {
    return(system2("sysctl", c("-n", "machdep.cpu.brand_string"),
                   stdout = TRUE))
  }
  NA_character_
}

sha <- Sys.getenv("GITHUB_SHA")
if (!nzchar(sha)) sha <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)

meta <- data.frame(
  commit   = substr(sha, 1L, 7L),
  date     = format(Sys.Date()),
  machine  = if (nzchar(Sys.getenv("GITHUB_ACTIONS"))) {
    paste("GitHub Actions", Sys.getenv("RUNNER_OS"), "runner")
  } else {
    "local"
  },
  cpu      = cpu_model(),
  cores    = parallel::detectCores(),
  reps     = reps,
  R        = paste(R.version$major, R.version$minor, sep = "."),
  rustc    = system2("rustc", "--version", stdout = TRUE),
  smotif   = as.character(utils::packageVersion("smotif")),
  smotifrs = as.character(utils::packageVersion("smotifrs"))
)
write.dcf(meta, file.path(out_dir, "meta.dcf"))
message("wrote ", out_dir)
