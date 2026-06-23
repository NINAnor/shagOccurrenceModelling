if (!file.exists(file.path(path_for_processed_data, "raster_stack_2023.parquet"))) {
  
  raster_results_23 <- create_harvest_dist_rasters(
    harvest = SK_35km_22_23,
    colony = colony,
    resolution = 10
  )
  
  save_raster_results(
    raster_results = raster_results_23,
    year = "2023",
    path_for_processed_data = path_for_processed_data
  )
  
} else {
  
  raster_results_23 <- load_raster_results(
    year = "2023",
    path_for_processed_data = path_for_processed_data,
    crs_reference = colony
  )
}


if (!file.exists(file.path(path_for_processed_data, "raster_stack_2024.parquet"))) {
  
  raster_results_24 <- create_harvest_dist_rasters(
    harvest = SK_35km_23_24,
    colony = colony,
    resolution = 10
  )
  
  save_raster_results(
    raster_results = raster_results_24,
    year = "2024",
    path_for_processed_data = path_for_processed_data
  )
  
} else {
  
  raster_results_24 <- load_raster_results(
    year = "2024",
    path_for_processed_data = path_for_processed_data,
    crs_reference = colony
  )
}