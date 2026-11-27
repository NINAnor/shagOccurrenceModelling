# Spatio-temporal LGCP of the daily harvest effect by sex, with a bird-level
# bootstrap (see vignettes/supplementary_inla.qmd).
#
# Requires in the environment (from masterScript / 01_prepData.R):
#   colony, dloc23, dloc, SK_35km_22_23, SK_35km_23_24, env_parquet_path,
#   modelCrs, result_dir, version_label
# Settings (defaults if not set): st70_B = 15, st70_threshold = 0.6,
#   st70_min_harvest_dives = 20, st70_min_birds = 3.

st70_B <- if (exists("st70_B")) st70_B else 15
st70_threshold <- if (exists("st70_threshold")) st70_threshold else 0.6
st70_min_harvest_dives <- if (exists("st70_min_harvest_dives")) {
  st70_min_harvest_dives
} else {
  20
}
st70_min_birds <- if (exists("st70_min_birds")) st70_min_birds else 3

library(INLA)
library(inlabru)
library(ggplot2)

# GEBCO is read from a local copy: GDAL's netcdf driver can fail on network drives
gebco_nc <- normalizePath(
  file.path(
    here::here(),
    "..",
    "..",
    "Analyses_Nina_multiple_Sites",
    "Environmental Covariates",
    "gebco_2019_n85.0_s50.0_w0.0_e40.0.nc"
  ),
  winslash = "/"
)
gebco_local <- file.path(tempdir(), "gebco.nc")
file.copy(gebco_nc, gebco_local, overwrite = TRUE)

env_stack_10m <- terra::rast(env_parquet_path)
cov35 <- build_inla_grid(
  colony = colony,
  tareskog = env_stack_10m$tareskog,
  gebco_path = gebco_local,
  modelCrs = modelCrs,
  buffer_m = 35000,
  res = 50
)
grid <- cov35$grid

dives_all <- prepare_dives(
  rbind(
    dloc23[, c("ID", "sex", "Logger_model", "UTC_datetime", "Year")],
    dloc[, c("ID", "sex", "Logger_model", "UTC_datetime", "Year")]
  ),
  modelCrs = modelCrs
)

# Pooled harvest: all cuts Sept 2022 - July 2024
hall <- dplyr::bind_rows(
  dplyr::select(SK_35km_22_23, GIS_TIME, year, geometry),
  dplyr::select(SK_35km_23_24, GIS_TIME, year, geometry)
)
hall <- hall[
  !duplicated(data.frame(
    t = hall$GIS_TIME,
    g = sf::st_as_text(sf::st_geometry(hall))
  )),
]

# Mesh: dive locations as nodes, 35 km study region as boundary
mesh_b <- create_fmesher_mesh_from_locations(
  points_sf = dives_all,
  study_region_sf = cov35$region,
  target_crs = modelCrs,
  offset = c(500, 1000),
  max_edge = c(15000, 30000),
  cutoff = 500
)
mesh <- mesh_b$mesh

# Dive covariates at exact 50 m cells (dives on non-water cells are dropped)
dv <- terra::extract(
  grid[[c("kelp", "depth", "dist", "water")]],
  cbind(dives_all$X, dives_all$Y)
)
keep <- !is.na(dv$water) & stats::complete.cases(dv)
dsf0 <- sf::st_sf(
  data.frame(
    kelp = as.numeric(dv$kelp[keep] > 0),
    depth = dv$depth[keep],
    dist = dv$dist[keep]
  ),
  geometry = sf::st_geometry(dives_all)[keep]
)
dsf0$depth2 <- dsf0$depth^2
dmeta <- sf::st_drop_geometry(dives_all)[
  keep,
  c("ID", "sex", "Year", "UTC_datetime")
]

# prep70$brk (distance-spline bins) and prep70$tmap (day -> time index) are
# fixed here and reused for every bootstrap resample
prep70 <- st70_prepare(dsf0, dmeta, grid, cov35$template, hall)

st70_out_dir <- file.path(result_dir, version_label, "inla_st70")
st70_boot_dir <- file.path(st70_out_dir, "bootstrap")
dir.create(st70_boot_dir, recursive = TRUE, showWarnings = FALSE)

# Mesh plot
p_mesh <- ggplot() +
  fmesher::geom_fm(
    data = mesh,
    alpha = 0,
    linewidth = 0.15,
    colour = "grey30"
  ) +
  geom_point(
    data = data.frame(X = mesh$loc[, 1], Y = mesh$loc[, 2]),
    aes(X, Y),
    size = 0.4,
    colour = "blue"
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "Mesh: dive locations as nodes, 35 km study region boundary",
    subtitle = sprintf("%d nodes", mesh$n)
  ) +
  theme_void()
ggsave(
  file.path(st70_out_dir, "mesh.png"),
  p_mesh,
  width = 6.5,
  height = 6,
  dpi = 150
)

# Full-data fit
res70s <- st70_fit_sex(prep70, mesh)
utils::write.csv(
  res70s$curves,
  file.path(st70_out_dir, "harvest_effect_daily_by_sex.csv"),
  row.names = FALSE
)
utils::write.csv(
  res70s$fit$summary.fixed,
  file.path(st70_out_dir, "fixed_effects.csv")
)
utils::write.csv(
  res70s$fit$summary.hyperpar,
  file.path(st70_out_dir, "hyperparameters.csv")
)
ggsave(
  file.path(st70_out_dir, "harvest_effect_daily_by_sex.png"),
  st70_plot_sex(res70s$curves),
  width = 10,
  height = 6.5,
  dpi = 150
)

# Bird bootstrap (resumes: existing resample files are skipped)
st70_bootstrap(
  dsf0,
  dmeta,
  grid,
  cov35$template,
  hall,
  mesh,
  prep70,
  ids = seq_len(st70_B),
  save_dir = st70_boot_dir
)

dm_t <- dmeta
dm_t$ID <- as.character(dm_t$ID)
st70_summary <- st70_bootstrap_summary(
  res70s$curves,
  prep70$dsf,
  dm_t,
  st70_boot_dir,
  threshold = st70_threshold,
  min_harvest_dives = st70_min_harvest_dives,
  min_birds = st70_min_birds
)
utils::write.csv(
  st70_summary$results,
  file.path(st70_out_dir, "daily_harvest_effect_by_sex_bootstrap.csv"),
  row.names = FALSE
)
utils::write.csv(
  st70_summary$zones,
  file.path(st70_out_dir, "robust_time_zones.csv"),
  row.names = FALSE
)
ggsave(
  file.path(st70_out_dir, "harvest_effect_daily_by_sex_bootstrap.png"),
  st70_summary$plot,
  width = 10,
  height = 6.5,
  dpi = 150
)
print(as.data.frame(st70_summary$zones))

# Data summary: dives per day in harvested vs not harvested cells, per year and
# sex (figure and tables saved with the st70 results; also shown in the Plots pane)
st70_counts <- st70_dive_counts(
  prep70, dmeta,
  n_raw = nrow(dloc23) + nrow(dloc), n_unique = nrow(dives_all),
  out_dir = st70_out_dir
)
