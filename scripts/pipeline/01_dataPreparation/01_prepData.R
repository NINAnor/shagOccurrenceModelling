# ============================================================
# Build used-available modelling dataset
# ============================================================
# Spatial reference:
# All spatial objects are assumed to be in WGS 84 / UTM zone 32N.
#
# Main outputs:
#   dmodel123 = processed model data for 2023
#   dmodel124 = processed model data for 2024
#   dmodel    = combined 2023 + 2024 model dataset
# ============================================================


# ------------------------------------------------------------
# 1. Load spatial and tracking data
# ------------------------------------------------------------

# Colony location
load("../colony")
study_region <- sf::st_buffer(colony, dist = 33000)

# Kelp forest polygons
load("../tare")

# Dive locations from 2024
# Object name after loading: dloc
load("../dloc24")

# Dive locations from 2023
# Object name after loading: dloc23
load("../dloc23")

# Kelp harvest information for 2022/2023,
# restricted to within 35 km of Sklinna
load("../SK_35_22_23")

# Kelp harvest information for 2023/2024,
# restricted to within 35 km of Sklinna
load("../SK_35_23_24")


if(!file.exists(file.path(path_for_processed_data, "dmodel.RData"))){
# ------------------------------------------------------------
# 2. Process 2023 dive locations
# ------------------------------------------------------------

result23 <- make_year_dmodel(
  dloc_year   = dloc23,
  tare        = tare,
  colony      = colony,
  harvest     = SK_35km_22_23,
  year        = "2023",
  n_per_point = n_per_point,
  min_dist_m  = 1000,
  ID_buffer_m = ID_buffer_m, 
  group_by_ID = TRUE
)

# Extract the processed 2023 model dataset
dmodel123 <- result23$dmodel


# ------------------------------------------------------------
# 3. Process 2024 dive locations
# ------------------------------------------------------------

result24 <- make_year_dmodel(
  dloc_year   = dloc,
  tare        = tare,
  colony      = colony,
  harvest     = SK_35km_23_24,
  year        = "2024",
  n_per_point = n_per_point,
  min_dist_m  = 1000,
  ID_buffer_m = ID_buffer_m, 
  group_by_ID = TRUE
)

# Extract the processed 2024 model dataset
dmodel124 <- result24$dmodel


# ------------------------------------------------------------
# 4. Diagnostic plots
# ------------------------------------------------------------

# 2023: dive locations clipped to kelp habitat
result23$dloctare$plot

# 2023: observed and available points
result23$avail_results$plot

# 2024: dive locations clipped to kelp habitat
result24$dloctare$plot

# 2024: observed and available points
result24$avail_results$plot


# ------------------------------------------------------------
# 5. Combine 2023 and 2024 datasets
# ------------------------------------------------------------

dmodel <- finalise_dmodel(
  dmodel23     = dmodel123,
  dmodel24     = dmodel124,
  datetime_col = "UTC_datetime"
)

save(
  dmodel,
  file = file.path(path_for_processed_data, "dmodel.RData")
)
} else {
  load(file.path(path_for_processed_data, "dmodel.RData"))
}
# ------------------------------------------------------------
# 6. Inspect final dataset
# ------------------------------------------------------------

if(checks){
summary(dmodel)
str(dmodel)

table(dmodel$Year, dmodel$random)
table(dmodel$Year, dmodel$harvested50m)
table(dmodel$Year, dmodel$harvested50m, dmodel$random)
}

# ------------------------------------------------------------
# Load GEBCO bathymetry
# ------------------------------------------------------------
if(!file.exists(file.path(path_for_processed_data, "env_stack_10m.parquet"))){
gebco <- terra::rast(
  "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_multiple_Sites/Environmental Covariates/gebco_2019_n85.0_s50.0_w0.0_e40.0.nc"
)

gebco_10m <- create_10m_raster(
  input_data = gebco,
  study_region = study_region,
  resolution = 10,
  raster_name = "bathymetry"
)

plot(gebco_10m)

# ------------------------------------------------------------
# Load hard-substrate shapefile
# ------------------------------------------------------------
hard_substr <- sf::st_read(
  "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_multiple_Sites/Environmental Covariates/Substrate_Rock_hard_substr.shp"
)

hard_substr_10m <- create_10m_raster(
  input_data = hard_substr,
  study_region = study_region,
  resolution = 10,
  field = 1,
  background = 0,
  raster_name = "hard_substrate"
)

plot(hard_substr_10m)

# ------------------------------------------------------------
# Load substrate_sediment shapefile
# ------------------------------------------------------------
sediment_substr <- sf::st_read(
  "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_multiple_Sites/Environmental Covariates/Substrate_Sediment.shp"
)

sediment_substr_10m <- create_10m_raster(
  input_data = sediment_substr,
  study_region = study_region,
  resolution = 10,
  field = 1,
  background = 0,
  raster_name = "sediment_substrate"
)

plot(sediment_substr_10m)
# ------------------------------------------------------------
# Load kelp forest shapefile
# ------------------------------------------------------------

tareskog <- sf::st_read(
  "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_multiple_Sites/Environmental Covariates/Tareskog_merge.shp"
)

tareskog_10m <- create_10m_raster(
  input_data = tareskog,
  study_region = study_region,
  resolution = 10,
  field = 1,
  background = 0,
  raster_name = "tareskog"
)

plot(tareskog_10m)


# ------------------------------------------------------------
# Combine the three 10 m rasters into one raster stack
# ------------------------------------------------------------

env_stack_10m <- c(
  gebco_10m,
  hard_substr_10m,
  sediment_substr_10m,
  tareskog_10m
)

names(env_stack_10m) <- c(
  "bathymetry",
  "hard_substrate",
  "sediment_substrate",
  "tareskog"
)

# Check
env_stack_10m
plot(env_stack_10m)


# ------------------------------------------------------------
# Convert raster stack to dataframe
# ------------------------------------------------------------

env_stack_df <- terra::as.data.frame(
  env_stack_10m,
  xy = TRUE,
  na.rm = FALSE
)


# ------------------------------------------------------------
# Save as parquet file
# ------------------------------------------------------------

arrow::write_parquet(
  env_stack_df,
  file.path(path_for_processed_data, "env_stack_10m.parquet")
)
} else {
  env_stack_df <- arrow::read_parquet(
    file.path(path_for_processed_data, "env_stack_10m.parquet")
  )
  
  env_stack_10m <- terra::rast(
    env_stack_df,
    type = "xyz",
    crs = sf::st_crs(study_region)$wkt
  )
  
  names(env_stack_10m) <- c(
    "bathymetry",
    "hard_substrate",
    "sediment_substrate",
    "tareskog"
  )
  
  #plot(env_stack_10m)
}


