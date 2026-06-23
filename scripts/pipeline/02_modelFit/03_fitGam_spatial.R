# ------------------------------------------------------------
# Add coordinates to model data
# ------------------------------------------------------------

model_data <- dmodel %>%
  #sf::st_drop_geometry() %>%
  dplyr::filter(sex == pars_to_use$sex,
                Year == pars_to_use$Year,
                harvesttimeSK24 == pars_to_use$harvesttimeSK24)  %>%
  dplyr::select(ID, Logger_model, random, Year, sex, harvesttimeSK24)

# ------------------------------------------------------------
# Extract raster values at model_data point locations
# ------------------------------------------------------------

# env_stack_10m should contain:
#     "bathymetry", "hard_substrate", "sediment_substrate", "tareskog"

raster_values <- terra::extract(
  env_stack_10m,
  terra::vect(model_data)
)

# Extract raster values at prediction-grid points
raster_values_pred_grid <- terra::extract(
  env_stack_10m,
  terra::vect(pred_grid)
) %>%
  dplyr::select(-ID)



if(pars_to_use$Year == "2023"){
raster_values_harvest_regions <- terra::extract(
  raster_results_23$raster_stack,
  terra::vect(model_data)
)


harvestRegionsPredGrid <- terra::extract(
  raster_results_23$raster_stack,
  terra::vect(pred_grid)
)%>%
  dplyr::select(-ID)

harvestRegions <- raster_results_23$harvest50_buffer

} else {
  raster_values_harvest_regions <- terra::extract(
    raster_results_24$raster_stack,
    terra::vect(model_data)
  )
  
  harvestRegionsPredGrid <- terra::extract(
    raster_results_24$raster_stack,
    terra::vect(pred_grid)
  )%>%
    dplyr::select(-ID)
  
  harvestRegions <- raster_results_24$harvest50_buffer
}

# Remove terra ID column
raster_values <- raster_values %>%
  dplyr::select(-ID)

raster_values_harvest_regions <- raster_values_harvest_regions %>%
  dplyr::select(-ID)

# Bind extracted raster values back to model_data
model_data_raster <- model_data %>%
  dplyr::bind_cols(raster_values,
                   raster_values_harvest_regions)

# Remove terra ID column and bind values back to pred_grid
pred_grid_env <- pred_grid  %>%
  dplyr::bind_cols(
    raster_values_pred_grid,
    harvestRegionsPredGrid
  )

# Check
# head(model_data_raster)
# 
# summary(model_data_raster[, c(
#   "bathymetry",
#   "hard_substrate",
#   "tareskog"
# )])

if(checks){
eda <- run_eda_by_year(
  dmodel = model_data_raster,
  years = pars_to_use$Year,
  response_col = "random",
  cov_cols = c("bathymetry", "DistNest", "hard_substrate", "sediment_substrate", "tareskog", "harvested50m", "harvested200m", "Logger_model", "harvesttimeSK24", "Year", "sex"),
  factor_vars = c("hard_substrate", "sediment_substrate", "tareskog", "harvested50m", "harvested200m", "Logger_model", "harvesttimeSK24", "Year", "sex")
)

yr <- as.character(pars_to_use$Year)
eda[[yr]]$missing_values
eda[[yr]]$frequency_tables
eda[[yr]]$cross_tables
eda[[yr]]$plots$distnest_hist
eda[[yr]]$plots$pairs
}


model_data <- model_data_raster %>%
  drop_na()

coords <- sf::st_coordinates(model_data)

model_data <-  model_data %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(
    x = coords[, 1],
    y = coords[, 2],
    random = as.numeric(as.character(random)),
    Logger_model = factor(Logger_model),
    harvested50m = factor(harvested50m),
    tareskog = factor(tareskog),
    hard_substrate = factor(hard_substrate),
    sediment_substrate = factor(sediment_substrate),
    ID = factor(ID),
    bathymetry_sc = as.numeric(scale(bathymetry)),
    DistNest_sc = as.numeric(scale(DistNest))
  )


gam_harvest_spatial <- mgcv::gam(
  random ~ harvested50m +
    s(DistNest_sc, k = 2) +
    s(bathymetry_sc, k = 2) +
    harvested50m +
    hard_substrate + 
    sediment_substrate + 
    tareskog +
    s(x, y, bs = "tp", k = 20) +
    s(ID, bs = "re"),
  data = model_data,
  family = binomial,
  method = "REML"
)

# ------------------------------------------------------------
# Convert SpatRaster to prediction dataframe
# ------------------------------------------------------------

pred_df <- terra::as.data.frame(
  pred_grid_env,
  xy = TRUE,
  na.rm = TRUE
)

# ------------------------------------------------------------
# Use the same DistNest scaling as the model data
# ------------------------------------------------------------

dist_center <- mean(model_data$DistNest, na.rm = TRUE)
dist_scale  <- sd(model_data$DistNest, na.rm = TRUE)

bathymetry_center <- mean(model_data$bathymetry, na.rm = TRUE)
bathymetry_scale <- sd(model_data$bathymetry, na.rm = TRUE)

pred_df <- pred_df %>%
  dplyr::mutate(
    harvested50m = factor(
      harvested50m,
      levels = levels(model_data$harvested50m)
    ),

    tareskog = factor(tareskog,
                      levels = levels(model_data$tareskog)
                      ),
    hard_substrate = factor(hard_substrate,
                            levels = levels(model_data$hard_substrate)
                            ),
    sediment_substrate = factor(sediment_substrate,
                                levels = levels(model_data$sediment_substrate)
                                ),
    
    DistNest_sc = as.numeric(
      (DistNest - dist_center) / dist_scale
    ),
    
    bathymetry_sc = as.numeric(
      (bathymetry - bathymetry_center) / bathymetry_scale
    ),
    
    ID = factor(
      levels(model_data$ID)[1],
      levels = levels(model_data$ID)
    )
  )

pred_df$pred_prob <- predict(
  gam_harvest_spatial,
  newdata = pred_df,
  type = "response",
  exclude = "s(ID)"
)

pred_raster <- terra::rast(
  pred_df[, c("x", "y", "pred_prob")],
  type = "xyz",
  crs = sf::st_crs(study_region)$wkt
)

names(pred_raster) <- "pred_prob"

library(sf)
library(ggplot2)
library(terra)

# ------------------------------------------------------------
# Keep only harvest buffers inside study region
# ------------------------------------------------------------

harvestRegions_sf <- sf::st_as_sf(harvestRegions)

# Match CRS if needed
if (sf::st_crs(harvestRegions_sf) != sf::st_crs(study_region)) {
  
  harvestRegions_sf <- sf::st_transform(
    harvestRegions_sf,
    sf::st_crs(study_region)
  )
}

# Clip harvest buffers to study region
harvestRegions_sf <- sf::st_intersection(
  harvestRegions_sf,
  sf::st_union(study_region)
)

# ------------------------------------------------------------
# Convert raster to plotting dataframe
# ------------------------------------------------------------

pred_df_plot <- terra::as.data.frame(
  pred_raster,
  xy = TRUE,
  na.rm = TRUE
)

# Rename prediction column if needed
colnames(pred_df_plot)[3] <- "pred_prob"

# ------------------------------------------------------------
# Plot prediction raster + clipped harvest regions
# ------------------------------------------------------------

ggplot() +
  
  geom_raster(
    data = pred_df_plot,
    aes(x = x,
        y = y,
        fill = pred_prob)
  ) +
  
  geom_sf(
    data = harvestRegions_sf,
    fill = NA,
    color = "red",
    linewidth = 0.5
  ) +
  
  geom_sf(
    data = study_region,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  
  coord_sf(
    crs = sf::st_crs(study_region)
  ) +
  
  scale_fill_viridis_c(
    name = "Predicted probability",
    na.value = NA
  ) +
  
  labs(
    title = "Predicted habitat suitability",
    subtitle = "Harvest buffers clipped to study region",
    x = "Easting",
    y = "Northing"
  ) +
  
  theme_minimal()

if(checks){
summary(gam_harvest_spatial)

library(ggeffects)

pred <- ggeffects::ggpredict(
  gam_harvest_spatial,
  terms = c(
    "harvested50m"
  )
)

plot(pred)


gam.check(gam_harvest_spatial)
}
