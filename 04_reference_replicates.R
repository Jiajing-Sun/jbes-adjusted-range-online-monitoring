# Resumable reference simulation. Each replication has its own L'Ecuyer stream;
# changing the number of workers does not change any saved reference draw.
reference_replicates <- function(fun, n_rep, n_cores, seed, directory, metadata) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  manifest <- file.path(directory, 'manifest.rds')
  signature <- list(
    n_rep = n_rep, seed = seed, metadata = metadata,
    implementation = tools::md5sum(file.path(ROOT, c('02_methods.R', '03_affine_extrema.R', '07_aersn_backend.R')))
  )
  if (file.exists(manifest)) stopifnot(identical(readRDS(manifest), signature)) else saveRDS(signature, manifest)
  RNGkind("L'Ecuyer-CMRG");set.seed(seed)
  streams <- vector('list', n_rep);streams[[1L]] <- .Random.seed
  if (n_rep > 1L) for (i in 2:n_rep) streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
  one_reference <- function(i) {
    assign('.Random.seed', streams[[i]], envir = .GlobalEnv);fun(i)
  }
  ans <- vector('list', n_rep)
  for (start in seq.int(1L, n_rep, by = 100L)) {
    ids <- seq.int(start, min(start + 99L, n_rep));path <- file.path(directory, sprintf('rep_%05d_%05d.rds', min(ids), max(ids)))
    if (file.exists(path)) batch <- readRDS(path) else {
      batch <- par_lapply_safe(ids, one_reference, n_cores = n_cores, seed = seed)
      tmp <- paste0(path, '.tmp');saveRDS(batch, tmp);stopifnot(file.rename(tmp, path))
    }
    stopifnot(length(batch) == length(ids));ans[ids] <- batch
    msg('Reference checkpoint: ', max(ids), ' / ', n_rep)
  }
  ans
}
