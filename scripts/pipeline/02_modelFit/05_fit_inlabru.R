library(sf)
library(terra)
library(fmesher)
library(INLA)
library(inlabru)



if (pars_to_use$Year == "2023") {
  if(!file.exists(file.path(path_for_processed_data, "cov_stack_proj_2023.tif"))){
  # create combined SpatRaster pointing to on-disk files (this avoids big in-memory copies)
  cov_stack <- terra::rast(c(file.path(path_for_processed_data, "raster_stack_2023.tif"),
                             file.path(path_for_processed_data, "env_stack_10m.tif")))

  cov_stack <- cov_stack |>
    tidyterra::select(DistNest, DistHarvest, bathymetry)

  cov_stack_proj <- terra::project(cov_stack,
                                   paste0("EPSG:",modelCrs),
                                   filename = file.path(path_for_processed_data, "cov_stack_proj_2023.tif"),
                                   overwrite = TRUE)
  }

  cov_stack <- terra::rast(file.path(path_for_processed_data, "cov_stack_proj_2023.tif"))

  # Harvest information (two datasets covering different year windows)
  load("../SK_35_22_23")    # defines 'SK_35km_22_23'

  harvest_region <- SK_35km_22_23 %>%
    sf::st_transform(modelCrs) %>%
    select(geometry) %>%
    sf::st_crop(., study_region %>% sf::st_transform(modelCrs))

} else {

  if(!file.exists(file.path(path_for_processed_data, "cov_stack_proj_2024.tif"))){
    # create combined SpatRaster pointing to on-disk files (this avoids big in-memory copies)
    cov_stack <- terra::rast(c(file.path(path_for_processed_data, "raster_stack_2024.tif"),
                               file.path(path_for_processed_data, "env_stack_10m.tif")))

    cov_stack <- cov_stack |>
      tidyterra::select(DistNest, DistHarvest, bathymetry)

    cov_stack_proj <- terra::project(cov_stack,
                                     paste0("EPSG:",modelCrs),
                                     filename = file.path(path_for_processed_data, "cov_stack_proj_2024.tif"),
                                     overwrite = TRUE)
  }

  cov_stack <- terra::rast(file.path(path_for_processed_data, "cov_stack_proj_2024.tif"))

  load("../SK_35_23_24")    # defines 'SK_35km_23_24'
  harvest_region <- SK_35km_23_24 %>%
    sf::st_transform(modelCrs) %>%
    select(geometry) %>%
    sf::st_crop(., study_region %>% sf::st_transform(modelCrs))

}


##--------------------------------------------------
## Presence data
##--------------------------------------------------

pres_sf <- dmodel |>
  dplyr::filter(
    Year == pars_to_use$Year,
    sex == pars_to_use$sex,
    harvesttimeSK24 == pars_to_use$harvesttimeSK24,
    random == 1
  )

# Keep the first row for each exact (X,Y) coordinate
pres_sf_unique <- pres_sf %>%
  dplyr::mutate(
    X = sf::st_coordinates(.)[,1],
    Y = sf::st_coordinates(.)[,2]
  ) %>%
  dplyr::distinct(X, Y, .keep_all = TRUE) %>%
  dplyr::select(-X, -Y)

# verify
sum(duplicated(sf::st_coordinates(pres_sf_unique)))
nrow(pres_sf); nrow(pres_sf_unique)

library(sf)

keep <- rep(TRUE, nrow(pres_sf_unique))

for (i in seq_len(nrow(pres_sf_unique) - 1)) {

  if (!keep[i]) next

  d <- st_distance(
    pres_sf_unique[i, ],
    pres_sf_unique[(i + 1):nrow(pres_sf_unique), ]
  )

  close <- which(as.numeric(d) < 10)

  if (length(close) > 0) {
    keep[(i + 1):nrow(pres_sf_unique)][close] <- FALSE
  }
}

pres_sf_5km <- pres_sf_unique[keep, ]
nrow(pres_sf_5km)

mesh_builder <- create_fmesher_mesh_from_locations(
  points_sf       = pres_sf_5km,
  study_region_sf = study_region,
  target_crs      = modelCrs,
  offset          = c(2000, 10000),
  max_edge        = c(3000, 10000),
  cutoff          = 1000
)

mesh      <- mesh_builder$mesh
region_sf <- mesh_builder$study_region
pres_sf   <- mesh_builder$used_points


spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(3000, 0.5),  # P(range < 1 km) = 0.05
  prior.sigma = c(0.5, 0.01)    # P(sigma > 0.5) = 0.01 (strong shrinkage)
)

cmp <-
  ~
  Intercept(1) +
  DistNest(main = DistNest,
           model = "linear") +
  DistHarvest(main = DistHarvest,
              model = "linear") +
  bathymetry(main = bathymetry,
             model = "linear") +
  spatial(
    geometry,
    model = spde
  )

cov_env <- list(
  DistNest = cov_stack$DistNest,
  DistHarvest = cov_stack$DistHarvest,
  bathymetry = cov_stack$bathymetry
)

# pres_sf$DistNest     <- terra::extract(cov_env$DistNest, pres_sf)[,2]
# pres_sf$DistHarvest  <- terra::extract(cov_env$DistHarvest, pres_sf)[,2]
# pres_sf$bathymetry   <- terra::extract(cov_env$bathymetry, pres_sf)[,2]
#
# fit <- bru(
#   components = cmp,
#   bru_obs(
#     formula = geometry ~ .,
#     family  = "cp",
#     data    = pres_sf,
#     domain  = list(geometry = mesh),
#     samplers = list(geometry = region_sf)
#   )
# )

#2. fit with informative prior on the intercept (use expected_intercept_m2)
expected_intercept_m2 <- log(nrow(pres_sf) / as.numeric(st_area(region_sf)))

pres_sf <- pres_sf %>%
  dplyr::select(random, geometry)


fit <- bru(
  components = cmp,

  bru_obs(
    formula  = geometry ~ .,
    family   = "cp",
    data     = pres_sf,

    domain   = list(
      geometry = mesh
    ),

    samplers = list(
      geometry = region_sf
    ),

    data_extra = cov_env
  ),
  options = bru_options(
   # inla = list(
      control.fixed = list(
        mean.intercept = expected_intercept_m2,
        prec.intercept = 1,
        mean = 0,
        prec = 1
     # )
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE
    ),
    bru_verbose = TRUE
  )
)


cmp_nospatial <-
  ~
  Intercept(1) +
  DistNest(
    main = DistNest,
    model = "linear"
  ) +
  DistHarvest(
    main = DistHarvest,
    model = "linear"
  ) +
  bathymetry(
    main = bathymetry,
    model = "linear"
  )

fit_nospatial <- bru(
  components = cmp_nospatial,

  bru_obs(
    formula = geometry ~ .,
    family = "cp",
    data = pres_sf,

    domain = list(
      geometry = mesh
    ),

    samplers = list(
      geometry = region_sf
    ),

    data_extra = cov_env
  ),

  options = bru_options(
    control.fixed = list(
      mean.intercept = expected_intercept_m2,
      prec.intercept = 1,
      mean = 0,
      prec = 1
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE
    ),
    bru_verbose = TRUE
  )
)

allModels <- list(fit,
                  fit_nospatial)

waicEsts <-lapply(allModels, function(x){
  x$waic$waic
})

minWAIC <- waicEsts %>%
  do.call("c", .) %>%
  which.min


bestModel <- allModels[[minWAIC]]

px <- fm_pixels(mesh,
                dim = c(500, 500),
                mask = region_sf)

pred_data <- px

pred_data$DistNest <- terra::extract(
  cov_stack$DistNest,
  terra::vect(pred_data)
)[, 2]

pred_data$DistHarvest <- terra::extract(
  cov_stack$DistHarvest,
  terra::vect(pred_data)
)[, 2]

pred_data$bathymetry <- terra::extract(
  cov_stack$bathymetry,
  terra::vect(pred_data)
)[, 2]

pred_dh <- predict(
  bestModel,
  pred_data,
  ~ exp(
    Intercept +
      DistHarvest
  )
)

pred_allCovs <- predict(
  bestModel,
  pred_data,
  ~ exp(
    Intercept +
      DistNest +
      DistHarvest +
      bathymetry
  )
)





plot1 <- ggplot(pred_dh) +
  geom_sf(aes(colour = mean)) +
  geom_sf(
    data = harvest_region,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_sf(
    data = pres_sf,
    colour = "#D55E00",
    size = 2
  ) +
  scale_colour_viridis_c(
    option = "cividis"
  ) +
  theme_bw()+
  labs(subtitle = paste("Year:", pars_to_use$Year, "| Sex: ", pars_to_use$sex, "| Harvest:", pars_to_use$harvesttimeSK24),
       title = paste("Distance to harvest prediction: ", pars_to_use$method))


plot2 <- ggplot(pred_allCovs) +
  geom_sf(aes(colour = mean)) +
  geom_sf(
    data = harvest_region,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_sf(
    data = pres_sf,
    colour = "#D55E00",
    size = 2
  ) +
  scale_colour_viridis_c(
    option = "cividis"
  ) +
  theme_bw()+
  labs(subtitle = paste("Year:", pars_to_use$Year, "| Sex: ", pars_to_use$sex, "| Harvest:", pars_to_use$harvesttimeSK24),
       title = paste("Full model predictions: ", pars_to_use$method))


