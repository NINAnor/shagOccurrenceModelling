#' Create a two-dimensional finite element mesh for spatial SPDE models
#'
#' Constructs a constrained Delaunay triangulation suitable for spatial
#' Gaussian random fields fitted with INLA or inlabru. Observation locations
#' are first reprojected into a projected coordinate reference system (CRS),
#' clipped to the study region, and used to generate a mesh bounded by the
#' supplied study area polygon.
#'
#' The mesh forms the basis of the SPDE approximation to the Matérn Gaussian
#' random field (Lindgren et al., 2011).
#'
#' @param points_sf An `sf` object containing point observations.
#'
#' @param study_region_sf An `sf` polygon defining the spatial domain over
#'   which the latent Gaussian field is estimated.
#'
#' @param target_crs Optional projected CRS used for mesh construction.
#'   If `NULL`, the CRS of `points_sf` is used. Geographic (longitude/latitude)
#'   coordinate systems are not permitted because mesh distances must be in
#'   linear units (typically metres).
#'
#' @param max_edge Numeric vector of length two specifying the maximum triangle
#'   edge length inside and outside the observation domain.
#'
#' @param cutoff Minimum allowed distance between mesh vertices. Observation
#'   points closer than this distance are merged before mesh generation.
#'
#' @param offset Numeric vector of length two specifying the extension of the
#'   mesh beyond the study region. A non-zero offset reduces boundary effects
#'   when estimating the spatial field.
#'
#' @details
#' The function performs the following steps:
#'
#' \enumerate{
#'   \item Validates inputs.
#'   \item Reprojects all spatial objects into a projected CRS.
#'   \item Unions and validates the study-region geometry.
#'   \item Removes observations outside the study region.
#'   \item Simplifies the study-region boundary to improve mesh stability.
#'   \item Generates a constrained triangulation using
#'         `fmesher::fm_mesh_2d()`.
#' }
#'
#' The resulting mesh can be used directly in
#' `INLA::inla.spde2.pcmatern()` or
#' `bru_mapper_spde()`.
#'
#' @return A list containing:
#'
#' \describe{
#'   \item{mesh}{Finite element mesh (`fm_mesh_2d`).}
#'   \item{used_points}{Observation points retained after clipping.}
#'   \item{study_region}{Validated study-region polygon.}
#'   \item{crs}{Projected coordinate reference system.}
#' }
#'
#' @references
#' Lindgren F, Rue H, Lindström J (2011).
#' An explicit link between Gaussian fields and Gaussian Markov random fields:
#' the stochastic partial differential equation approach.
#' *Journal of the Royal Statistical Society: Series B*,
#' **73**, 423--498.
#'
#' @export
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

  region_simple <-
    region_union |>
    sf::st_simplify(dTolerance = 500)


  ## ---- Boundary ----
  boundary <- fmesher::fm_as_segm(region_simple)

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

# choose_pc_priors_binomial <- function(coords,
#                                       range_fraction = 100, #0.2,
#                                       p_range = 0.5,
#                                       s0 = 1,
#                                       p_sigma = 0.01) {
#
#   coords <- as.matrix(coords)
#
#   max_dist <- max(dist(coords))
#   min_dist <- min(dist(coords))
#
#   r0 <- range_fraction * min_dist #range_fraction * max_dist
#   if(r0 > 3000) r0 <- 3000
#
#   list(
#     prior.range = c(r0, p_range),
#     prior.sigma = c(s0, p_sigma),
#     max.distance = max_dist,
#     min_dist = min_dist
#   )
# }

#' Choose penalised complexity priors for an SPDE model
#'
#' Computes Penalised Complexity (PC) priors for the Matérn spatial random
#' field used in INLA/inlabru SPDE models.
#'
#' The prior range is derived from the spatial extent of the study region,
#' producing weakly informative priors that scale naturally with the size
#' of the study area.
#'
#' @param study_region An `sf` polygon defining the modelling domain.
#'
#' @param range_fraction Positive number controlling the prior spatial range.
#' The reference range is calculated as
#'
#' `extent / range_fraction`
#'
#' where *extent* is the maximum width or height of the study region.
#'
#' Smaller values produce larger prior spatial ranges.
#'
#' @param p_range Probability assigned to
#'
#' `P(range < r0) = p_range`
#'
#' where `r0` is the computed reference range.
#'
#' @param s0 Reference marginal standard deviation of the latent Gaussian
#' field.
#'
#' @param p_sigma Probability assigned to
#'
#' `P(sigma > s0) = p_sigma`.
#'
#' @details
#'
#' PC priors shrink the spatial model towards the simpler base model
#' (no spatial variation) unless supported by the data
#' (Simpson et al., 2017).
#'
#' Rather than deriving the prior range from the spacing between observations,
#' this implementation scales the prior using the study-region extent.
#' This makes the prior independent of sampling density while remaining
#' proportional to the ecological domain being modelled.
#'
#' For example, if the study region is 60 km across and
#'
#' `range_fraction = 20`
#'
#' then
#'
#' `r0 = 3000 m`
#'
#' and the prior becomes
#'
#' `P(range < 3000 m) = 0.5`.
#'
#' @return A list containing
#'
#' \describe{
#'   \item{prior.range}{Vector passed directly to
#'   `INLA::inla.spde2.pcmatern()`.}
#'
#'   \item{prior.sigma}{Vector passed directly to
#'   `INLA::inla.spde2.pcmatern()`.}
#' }
#'
#' @references
#' Simpson D, Rue H, Riebler A, Martins TG, Sørbye SH (2017).
#' Penalising model component complexity:
#' A principled, practical approach to constructing priors.
#' *Statistical Science*,
#' **32**, 1--28.
#'
#' Fuglstad GA, Simpson D, Lindgren F, Rue H (2019).
#' Constructing priors that penalize the complexity of Gaussian random fields.
#' *Journal of the American Statistical Association*,
#' **114**, 445--452.
#'
#' @export
choose_pc_priors_binomial <- function(study_region,
                                      range_fraction = 20,
                                      p_range = 0.5,
                                      s0 = 1,
                                      p_sigma = 0.01) {

  bb <- sf::st_bbox(study_region)

  extent <- max(
    bb["xmax"] - bb["xmin"],
    bb["ymax"] - bb["ymin"]
  )

  r0 <- extent / range_fraction

  list(
    prior.range = c(r0, p_range),
    prior.sigma = c(s0, p_sigma)
  )
}

