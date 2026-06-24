#' Configure threading for terra, OpenMP, GDAL, and BLAS
#'
#' This helper centralizes parallel configuration for spatial workflows that use
#' terra (and underlying GDAL/OpenMP/BLAS). It:
#'  - sets environment variables OMP_NUM_THREADS, OMP_THREAD_LIMIT, GDAL_NUM_THREADS
#'  - optionally limits BLAS threads via RhpcBLASctl (if installed)
#'  - sets options terra.cores and mc.cores for convenience
#'
#' Notes:
#'  - Keep BLAS threads low (default 1) to avoid oversubscription when terra or GDAL
#'    use multiple threads.
#'  - GDAL accepts either an integer or the special string "ALL_CPUS" on some systems.
#'
#' @param n_terra Integer: target number of cores for terra operations (>= 1).
#' @param omp_threads Integer or NULL: OpenMP threads. If NULL, set to n_terra.
#' @param gdal_threads Integer, "ALL_CPUS", or NULL: GDAL threads. If NULL, set to n_terra.
#' @param blas_threads Integer: number of BLAS threads to allow (default 1).
#'
#' @return Invisibly returns a named list with the configured thread counts.
#' @export
set_terra_threads <- function(
    n_terra,                  # default cores for terra functions
    omp_threads   = NULL,     # OpenMP threads (NULL -> n_terra)
    gdal_threads  = NULL,     # GDAL threads (NULL -> n_terra)
    blas_threads  = 1         # keep BLAS low to avoid oversubscription
) {
  stopifnot(is.numeric(n_terra), length(n_terra) == 1, n_terra >= 1)
  n_terra <- as.integer(n_terra)

  # Detect available logical cores; if NA, continue without capping
  avail_logical <- parallel::detectCores(logical = TRUE)

  # Cap requested n_terra by available logical cores when possible
  if (!is.na(avail_logical)) n_terra <- min(n_terra, avail_logical)

  # Default omp/gdal to n_terra if NULL
  if (is.null(omp_threads))  omp_threads  <- n_terra
  if (is.null(gdal_threads)) gdal_threads <- n_terra

  # If avail_logical provided, cap omp and gdal (unless gdal is set to "ALL_CPUS")
  if (!is.na(avail_logical)) {
    omp_threads <- min(as.integer(omp_threads), avail_logical)
    if (!identical(gdal_threads, "ALL_CPUS")) {
      gdal_threads <- min(as.integer(gdal_threads), avail_logical)
    }
  }

  # Set OpenMP environment variables
  Sys.setenv(
    OMP_NUM_THREADS  = as.character(omp_threads),
    OMP_THREAD_LIMIT = as.character(omp_threads)
  )

  # Set GDAL threading (string accepted)
  Sys.setenv(GDAL_NUM_THREADS = as.character(gdal_threads))

  # Configure BLAS threading if RhpcBLASctl is available; otherwise set env vars
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    # Limit BLAS worker threads to avoid nested parallelism
    tryCatch({
      RhpcBLASctl::blas_set_num_threads(as.integer(blas_threads))
      # Optionally align OpenMP threads visible to BLAS libraries
      suppressWarnings(RhpcBLASctl::omp_set_num_threads(as.integer(omp_threads)))
    }, error = function(e) {
      warning("RhpcBLASctl present but failed to set threads: ", conditionMessage(e))
    })
  }

  # Also set common BLAS env vars for packages that inspect them
  Sys.setenv(
    MKL_NUM_THREADS        = as.character(blas_threads),
    OPENBLAS_NUM_THREADS   = as.character(blas_threads),
    VECLIB_MAXIMUM_THREADS = as.character(blas_threads),
    NUMEXPR_NUM_THREADS    = as.character(blas_threads)
  )

  # Set convenient R options many code paths check
  options(
    mc.cores    = n_terra,
    terra.cores = n_terra
  )

  # Report configuration to user
  message(sprintf(
    "Configured: terra.cores=%d | OMP_NUM_THREADS=%s | GDAL_NUM_THREADS=%s | BLAS threads=%s",
    getOption("terra.cores"),
    Sys.getenv("OMP_NUM_THREADS"),
    Sys.getenv("GDAL_NUM_THREADS"),
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) "unknown")
    } else {
      "unknown"
    }
  ))

  invisible(list(
    terra_cores = n_terra,
    omp_threads = omp_threads,
    gdal_threads = gdal_threads,
    blas_threads = blas_threads
  ))
}

# Convenience wrapper functions that use options("terra.cores") by default.
# These thin wrappers reduce repetitive boilerplate when calling terra functions.

# Apply a function over layers/cells using terra::app with default cores
t_app <- function(x, fun, ..., cores = getOption("terra.cores", 1)) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install 'terra' to use t_app()")
  terra::app(x, fun, ..., cores = cores)
}

# Focal operation wrapper using terra::focal
t_focal <- function(x, w, fun, ..., cores = getOption("terra.cores", 1)) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install 'terra' to use t_focal()")
  terra::focal(x, w = w, fun = fun, ..., cores = cores)
}

# Temporal aggregate (tapp) wrapper
t_tapp <- function(x, index, fun, ..., cores = getOption("terra.cores", 1)) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install 'terra' to use t_tapp()")
  terra::tapp(x, index = index, fun = fun, ..., cores = cores)
}

# Predict wrapper that forwards cores to terra::predict
t_predict <- function(x, model, ..., cores = getOption("terra.cores", 1)) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install 'terra' to use t_predict()")
  terra::predict(x, model, ..., cores = cores)
}

# Helper to construct GDAL-aware terra write options for multithreaded raster output.
# comp: compression algorithm (e.g., "LZW", "DEFLATE"), threads: number or "ALL_CPUS"
terra_wopt_mt <- function(comp = "LZW", threads = Sys.getenv("GDAL_NUM_THREADS", "ALL_CPUS")) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install 'terra' to use terra_wopt_mt()")
  predictor <- if (comp %in% c("LZW", "DEFLATE")) 2 else 1
  terra::wopt(gdal = c(
    paste0("COMPRESS=", comp),
    paste0("PREDICTOR=", predictor),
    paste0("NUM_THREADS=", threads)
  ))
}
