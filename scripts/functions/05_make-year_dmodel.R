# ============================================================
# Generate full modelling dataset for one year
# ============================================================

make_year_dmodel <- function(dloc_year,
                             tare,
                             colony,
                             harvest,
                             year,
                             n_per_point = 5,
                             min_dist_m = 500,
                             ID_buffer_m = 10000, 
                             group_by_ID = TRUE) {
  
  library(dplyr)
  library(sf)
  library(ggplot2)
  
  # ----------------------------------------------------------
  # 1. Restrict dive locations to kelp habitat within buffer
  # ----------------------------------------------------------
  
  dloctare <- process_dloctare(
    dloc = dloc_year,
    tare = tare,
    colony = colony,
    year = year
  )
  
  dloctare_resolved <- dloctare$dloctare_resolved %>%
    dplyr::mutate(
      naturtype = as.factor(naturtype)
    )
  
  
  # ----------------------------------------------------------
  # 2. Generate available random points
  # ----------------------------------------------------------
  
  avail_results <- generate_available_points(
    dloc = dloctare_resolved,
    study_region = dloctare$SK33km,
    tare = tare,
    colony = colony,
    n_per_point = n_per_point,
    min_dist_m = min_dist_m,
    id_buffer_m = ID_buffer_m,
    group_by_id = group_by_ID
  )
  
  combined_data <- avail_results$combined_data
  
  
  # ----------------------------------------------------------
  # 3. Add harvest overlap variables
  # ----------------------------------------------------------
  
  dmodel_year <- add_harvest_overlap(
    dmodel = combined_data,
    harvest = harvest %>% dplyr::filter(year == !!year),
    radii_m = c(200, 50),
    time_col = "GIS_TIME",
    datetime_col = "UTC_datetime",
    id_col = "seq"
  )
  
  
  # ----------------------------------------------------------
  # 4. Return objects for checking
  # ----------------------------------------------------------
  
  return(list(
    dloctare = dloctare,
    dloctare_resolved = dloctare_resolved,
    avail_results = avail_results,
    dmodel = dmodel_year
  ))
}