# ============================================================
# Function: create_10m_raster
#
# Purpose:
#   Load, project, crop, mask, and rasterise spatial data
#   to a common 10 m raster grid within the study region.
#
# Supports:
#   - vector polygon/line/point data (.shp)
#   - raster data (.nc, .tif, etc.)
#
# Inputs:
#   input_data   = sf object OR terra SpatRaster
#   study_region = sf polygon defining study extent
#   resolution   = raster resolution in meters
#   field        = value to rasterise for vector data
#   background   = background value outside polygons
#   raster_name  = output raster layer name
#
# Output:
#   10 m SpatRaster aligned to study region
# ============================================================

create_10m_raster <- function(input_data,
                              study_region,
                              resolution = 10,
                              field = 1,
                              background = 0,
                              raster_name = "layer",
                              raster_method = "bilinear") {
  
  library(sf)
  library(terra)
  
  # ----------------------------------------------------------
  # Study-region CRS
  # ----------------------------------------------------------
  
  study_crs <- sf::st_crs(study_region)$wkt
  
  # ----------------------------------------------------------
  # Union study region
  # ----------------------------------------------------------
  
  study_union <- sf::st_union(study_region)
  
  study_vect <- terra::vect(study_union)
  
  # ----------------------------------------------------------
  # Create 10 m template raster
  # ----------------------------------------------------------
  
  template_10m <- terra::rast(
    extent = terra::ext(study_vect),
    resolution = resolution,
    crs = study_crs
  )
  
  # ==========================================================
  # CASE 1: INPUT IS SPATRASTER
  # ==========================================================
  
  if (inherits(input_data, "SpatRaster")) {
    
    # --------------------------------------------------------
    # Reproject raster
    # --------------------------------------------------------
    
    input_proj <- terra::project(
      input_data,
      study_crs,
      method = raster_method
    )
    
    # --------------------------------------------------------
    # Crop and mask
    # --------------------------------------------------------
    
    input_crop <- terra::crop(
      input_proj,
      study_vect
    )
    
    input_crop <- terra::mask(
      input_crop,
      study_vect
    )
    
    # --------------------------------------------------------
    # Resample to 10 m
    # --------------------------------------------------------
    
    out_raster <- terra::resample(
      input_crop,
      template_10m,
      method = raster_method
    )
  }
  
  # ==========================================================
  # CASE 2: INPUT IS SF OBJECT
  # ==========================================================
  
  else if (inherits(input_data, "sf")) {
    
    # --------------------------------------------------------
    # Reproject vector
    # --------------------------------------------------------
    
    input_proj <- sf::st_transform(
      input_data,
      crs = study_crs
    )
    
    # --------------------------------------------------------
    # Crop to study region
    # --------------------------------------------------------
    
    input_crop <- sf::st_intersection(
      input_proj,
      study_union
    )
    
    # --------------------------------------------------------
    # Rasterise
    # --------------------------------------------------------
    
    out_raster <- terra::rasterize(
      terra::vect(input_crop),
      template_10m,
      field = field,
      background = background,
      touches = TRUE
    )
    
    # --------------------------------------------------------
    # Mask to study region
    # --------------------------------------------------------
    
    out_raster <- terra::mask(
      out_raster,
      study_vect
    )
  }
  
  # ==========================================================
  # INVALID INPUT
  # ==========================================================
  
  else {
    
    stop(
      "input_data must be either:
       - terra SpatRaster
       - sf object"
    )
  }
  
  # ----------------------------------------------------------
  # Set raster name
  # ----------------------------------------------------------
  
  names(out_raster) <- raster_name
  
  return(out_raster)
}