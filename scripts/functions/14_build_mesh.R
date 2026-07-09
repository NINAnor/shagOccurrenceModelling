create_fmesher_mesh_from_locations <- function(points_sf,
                                               study_region_sf,
                                               target_crs = NULL,
                                               max_edge = c(Inf, Inf),
                                               cutoff = 0,
                                               offset = c(0, 0)) {

  if (!inherits(points_sf, "sf"))
    stop("points_sf must be an sf object.")
  if (!inherits(study_region_sf, "sf"))
    stop("study_region_sf must be an sf object.")
  if (nrow(points_sf) == 0)
    stop("points_sf is empty.")
  if (nrow(study_region_sf) == 0)
    stop("study_region_sf is empty.")

  ## ---- Choose CRS ----
  if (!is.null(target_crs)) {
    crs_use <- sf::st_crs(target_crs)
    if (is.na(crs_use))
      stop("target_crs not recognized.")
  } else {
    if (is.na(sf::st_crs(points_sf)))
      stop("points_sf has no CRS. Provide target_crs.")

    crs_use <- sf::st_crs(points_sf)

    if (sf::st_is_longlat(crs_use))
      stop(
        "points_sf is in lon/lat. Supply a projected CRS (e.g. UTM) via target_crs."
      )
  }

  ## ---- Reproject ----
  pts <- sf::st_transform(points_sf, crs_use)
  region <- sf::st_transform(study_region_sf, crs_use)

  ## ---- Valid study region ----
  region_union <- region |>
    sf::st_union() |>
    sf::st_make_valid()

  ## ---- Keep points inside ----
  inside <- lengths(sf::st_intersects(pts, region_union)) > 0
  pts_in <- pts[inside, ]

  if (nrow(pts_in) < 3)
    stop("Need at least 3 points inside the study region.")

  ## ---- Boundary ----
  boundary <- fmesher::fm_as_segm(region_union)

  ## ---- Mesh ----
  mesh <- fmesher::fm_mesh_2d(
    loc      = sf::st_coordinates(pts_in),
    boundary = boundary,
    max.edge = max_edge,
    cutoff   = cutoff,
    offset   = offset,
    crs      = crs_use
  )

  list(
    mesh = mesh,
    used_points = pts_in,
    study_region = region_union,
    crs = crs_use
  )
}
