study_region_details <- function(study_region){
  study_region <- sf::st_buffer(colony, dist = 33000 + 10000)  # 33 km study buffer
  # Ensure study_region is a single (unioned) polygon and obtain CRS WKT
  study_union <- sf::st_union(study_region)
  study_crs_wkt <- sf::st_crs(study_region)$wkt
  if (is.na(study_crs_wkt) || study_crs_wkt == "") stop("study_region must have a valid projected CRS (with meter units)")

  # Convert study region to terra vect for raster operations

  study_union_sf <- sf::st_sf(
    geometry = study_union
  )

  study_vect <- terra::vect(study_union_sf)

  #study_vect <- terra::vect(study_union)

  # Create template raster aligned to study region extent and CRS
  template <- terra::rast(
    extent = terra::ext(study_vect),
    resolution = 10,
    crs = study_crs_wkt
  )

  ret <- list(study_region = study_region,
              study_union = study_union,
              study_vect = study_vect,
              template = template)

  return(ret)
}



#' Create a 10m (or user-specified resolution) raster layer from vector or raster input
#'
#' This helper loads, reprojects, crops/masks, and (if needed) rasterises input
#' spatial data to a common raster grid defined by a study region. It supports
#' both input terra SpatRaster objects and sf vector objects.
#'
#' Behaviour summary:
#'  - If input_data is a SpatRaster: reproject -> crop -> mask -> resample to template.
#'  - If input_data is an sf object: transform -> intersect with study region ->
#'    rasterise to template -> mask.
#'  - The returned SpatRaster is aligned to the study_region extent and CRS.
#'
#' Inputs:
#'  - input_data   : terra SpatRaster OR sf object (vector)
#'  - study_region : sf polygon defining the study extent (CRS used for output)
#'  - resolution   : raster resolution in meters (default 10)
#'  - field        : field/value used for rasterisation (vector inputs). For
#'                   categorical fields supply the column name or numeric index.
#'  - background   : background value for rasterised areas (default 0)
#'  - raster_name  : name assigned to the output raster layer
#'  - raster_method: resampling method for raster reprojection/resample (default "bilinear")
#'
#' Returns:
#'  - a terra SpatRaster aligned to the study_region with one layer named raster_name
#'
#' @export
create_10m_raster <- function(input_data,
                              study_region_list,
                              resolution = 10,
                              field = 1,
                              background = 0,
                              raster_name = "layer",
                              raster_method = "bilinear") {

  # ------------------------
  # Dependency checks
  # ------------------------
  needed <- c("sf", "terra")
  missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Install missing packages: ", paste(missing_pkgs, collapse = ", "))

  # ------------------------
  # Validate inputs
  # ------------------------
  if (!inherits(study_region_list, "list")) stop("study_region must be an sf object")
  if (!(inherits(input_data, "SpatRaster") || inherits(input_data, "sf"))) {
    stop("input_data must be either a terra SpatRaster or an sf object")
  }



  # Ensure study_region is a single (unioned) polygon and obtain CRS WKT
  study_union <- sf::st_union(study_region_list$study_union)
   study_crs_wkt <- sf::st_crs(study_region_list$study_region)$wkt
  # if (is.na(study_crs_wkt) || study_crs_wkt == "") stop("study_region must have a valid projected CRS (with meter units)")
  #
  # # Convert study region to terra vect for raster operations
  #
  # study_union_sf <- sf::st_sf(
  #   geometry = study_union
  # )
  #
  # study_vect <- terra::vect(study_union_sf)

  study_vect <- study_region_list$study_vect

  # Create template raster aligned to study region extent and CRS
  template <- study_region_list$template
  # template <- terra::rast(
  #   extent = terra::ext(study_vect),
  #   resolution = resolution,
  #   crs = study_crs_wkt
  # )

  # ------------------------
  # Case A: input is a SpatRaster
  # ------------------------
  if (inherits(input_data, "SpatRaster")) {

    # If input raster and template have different CRS, reproject to study CRS
    # Use terra::project; target CRS can be provided as WKT
    input_proj <- tryCatch(
      terra::project(input_data, study_crs_wkt, method = raster_method),
      error = function(e) stop("Failed to reproject input SpatRaster: ", conditionMessage(e))
    )

    # Crop and mask to study extent
    input_crop <- terra::crop(input_proj,  study_vect)
    input_mask <- terra::mask(input_crop, study_vect)

    # Resample to the template resolution & alignment
    out_raster <- terra::resample(input_mask, template, method = raster_method)

  } else if (inherits(input_data, "sf")) {

    # ------------------------
    # Case B: input is vector (sf)
    # ------------------------
    # Reproject vector to study CRS
    input_proj <- sf::st_transform(input_data, crs = study_crs_wkt)

    # Clip to study region; intersection may produce empty result
    input_clip <- sf::st_intersection(input_proj, study_union)

    if (nrow(input_clip) == 0) {
      # No features inside study region: return a raster filled with background
      out_raster <- template
      terra::values(out_raster) <- background
      names(out_raster) <- raster_name
      return(out_raster)
    }

    # Convert to terra vector and rasterise.
    # field may be a column name (character) or numeric value/index.
    vrt <- terra::vect(input_clip)

    # If field is a character column name and not present, error
    if (is.character(field) && !(field %in% names(input_clip))) {
      stop("Field '", field, "' not found in input sf object")
    }

    # Rasterize to the template: value = field, background outside polygons
    out_raster <- terra::rasterize(vrt, template, field = field, background = background, touches = TRUE)

    # Ensure masking to study region (in case rasterize spilled outside)
    out_raster <- terra::mask(out_raster, study_vect)

  } else {
    stop("Unhandled input_data type. Provide either a terra SpatRaster or an sf object.")
  }

  # ------------------------
  # Finalise: set layer name and return
  # ------------------------
  names(out_raster) <- raster_name
  out_raster
}
