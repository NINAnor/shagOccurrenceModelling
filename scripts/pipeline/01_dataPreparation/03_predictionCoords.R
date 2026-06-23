library(sf)
library(dplyr)
library(ggplot2)

# Grid spacing in meters
grid_res <- 100

if(!file.exists(file.path(path_for_processed_data, "pred_grid.parquet"))){
  
  # --------------------------------------------------
  # 1. Union all tareskog polygons
  # --------------------------------------------------
  
  tare_union <- tare %>%
    sf::st_union()
  
# Create regular grid
pred_grid <- sf::st_make_grid(
  tare_union,
  cellsize = grid_res,
  what = "centers"
)

# Keep only points inside tare region
pred_grid <- pred_grid[
  sf::st_intersects(pred_grid, tare_union, sparse = FALSE)
]

# Convert to sf object
pred_grid <- sf::st_sf(
  geometry = pred_grid,
  crs = sf::st_crs(study_region)
)

# Add coordinates if needed for GAM prediction
coords <- sf::st_coordinates(pred_grid)

# Save sf object as parquet
pred_grid_parquet <- pred_grid %>%
  dplyr::mutate(
    geometry_wkt = sf::st_as_text(sf::st_geometry(.))
  ) %>%
  sf::st_drop_geometry()

arrow::write_parquet(
  pred_grid_parquet,
  file.path(path_for_processed_data, "pred_grid.parquet")
)
} else {

pred_grid <- arrow::read_parquet(
  file.path(path_for_processed_data, "pred_grid.parquet")
)

pred_grid <- pred_grid %>%
  sf::st_as_sf(
    wkt = "geometry_wkt",
    crs = sf::st_crs(study_region)
  )
}

if(checks){
pred_grid$x <- coords[,1]
pred_grid$y <- coords[,2]

# Plot
ggplot() +
  geom_sf(data = study_region,
          fill = "grey90",
          color = "grey50") +
  geom_sf(data = pred_grid,
          color = "blue",
          size = 0.2) +
  theme_minimal() +
  labs(
    title = "Regular prediction grid",
    subtitle = paste0("Grid resolution: ", grid_res, " m")
  )
}