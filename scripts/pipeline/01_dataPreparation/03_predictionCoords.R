# Create or load a regular prediction grid (centers) inside the kelp (tareskog) mask
#
# Behavior:
#  - Builds a regular grid of points with spacing `grid_res` (meters) covering the
#    union of the tareskog polygons, retains only points that fall inside the
#    tareskog geometry, and saves the result as a parquet table for quick reload.
#  - On subsequent runs, reads the parquet and reconstructs an sf object from WKT.
#
# Requirements / assumptions:
#  - Objects `tare` (sf polygons) and `study_region` (sf polygon) exist in the
#    calling environment and share a projected CRS with meter units.
#  - `path_for_processed_data` is defined and writable.
#  - The `arrow` package is available to write/read parquet files (CSV fallback not implemented).
#
# Output:
#  - pred_grid : sf points (centers) inside tareskog, saved/loaded from pred_grid.parquet
#
# Optional diagnostic plot produced when `checks` is TRUE.
# ------------------------------------------------------------------------------

# Parameters
grid_res <- 100   # grid spacing in meters

# Dependency checks
req <- c("sf", "dplyr", "arrow", "ggplot2")
missing_pkgs <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))

# Ensure environment variables exist
if (!exists("path_for_processed_data")) stop("path_for_processed_data must be defined (output directory).")
if (!exists("tare")) stop("`tare` (kelp polygons) must be present in the environment.")
if (!exists("study_region")) stop("`study_region` must be present in the environment.")

dir.create(path_for_processed_data, recursive = TRUE, showWarnings = FALSE)
parquet_path <- file.path(path_for_processed_data, "pred_grid.parquet")

if (!file.exists(parquet_path)) {
  message("No pred_grid parquet found; building grid at ", grid_res, " m resolution...")

  # 1) Union tareskog polygons into a single geometry to define extent
  tare_union <- sf::st_union(tare)

  # 2) Create grid of cell centers covering the union geometry
  grid_pts <- sf::st_make_grid(
    tare_union,
    cellsize = grid_res,
    what = "centers"
  )

  # 3) Keep only points that intersect the tare_union (inside kelp polygons)
  inside_idx <- sf::st_intersects(grid_pts, tare_union, sparse = FALSE)[, 1]
  grid_pts <- grid_pts[inside_idx]

  # 4) Convert to sf and attach CRS of study_region (assumes projected CRS in meters)
  pred_grid <- sf::st_sf(geometry = grid_pts, crs = sf::st_crs(study_region))

  # 5) Add coordinate columns for convenience (useful for model prediction)
  coords <- sf::st_coordinates(pred_grid)
  pred_grid$x <- coords[, 1]
  pred_grid$y <- coords[, 2]

  # 6) Persist as parquet: store all attribute columns and a WKT geometry column
  pred_grid_parquet <- pred_grid %>%
    dplyr::mutate(geometry_wkt = sf::st_as_text(sf::st_geometry(.))) %>%
    sf::st_drop_geometry()

  arrow::write_parquet(pred_grid_parquet, parquet_path)
  message("Saved pred grid to: ", parquet_path)

} else {
  message("Loading pred_grid from: ", parquet_path)
  pred_grid_parquet <- arrow::read_parquet(parquet_path)

  # Reconstruct sf object from WKT geometry column
  if (!"geometry_wkt" %in% names(pred_grid_parquet)) stop("parquet missing 'geometry_wkt' column")
  pred_grid <- sf::st_as_sf(pred_grid_parquet, wkt = "geometry_wkt", crs = sf::st_crs(study_region))

  # Recompute coordinate columns if absent
  if (!all(c("x", "y") %in% names(pred_grid))) {
    coords <- sf::st_coordinates(pred_grid)
    pred_grid$x <- coords[, 1]
    pred_grid$y <- coords[, 2]
  }
}

# Optional diagnostic plot (only when `checks` exists and is TRUE)
if (exists("checks") && isTRUE(checks)) {
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = study_region, fill = "grey90", color = "grey50") +
    ggplot2::geom_sf(data = pred_grid, color = "blue", size = 0.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Regular prediction grid (centers inside tareskog)",
      subtitle = paste0("Grid resolution: ", grid_res, " m")
    ) -> p_pred_grid

  print(p_pred_grid)
}

# The object `pred_grid` (sf) is now available in the environment.
