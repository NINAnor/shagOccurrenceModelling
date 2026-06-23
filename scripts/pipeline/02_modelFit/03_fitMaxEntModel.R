library(dplyr)
library(sf)
library(maxnet)
library(pROC)

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------

maxent_data <- dmodel %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(
    presence = as.numeric(as.character(used)),   # 1 = used, 0 = available
    sex = factor(sex),
    Year = factor(Year),
    Logger_model = factor(Logger_model),
    harvesttimeSK24 = factor(harvesttimeSK24),
    harvested50m = factor(harvested50m),
    DistNest_sc = as.numeric(scale(DistNest))
  ) %>%
  dplyr::select(
    presence,
    sex,
    Year,
    Logger_model,
    harvesttimeSK24,
    harvested50m,
    DistNest_sc
  ) %>%
  tidyr::drop_na()

# ------------------------------------------------------------
# Fit MaxEnt model
# ------------------------------------------------------------

x <- maxent_data %>%
  dplyr::select(
    sex,
    Year,
    Logger_model,
    harvesttimeSK24,
    harvested50m,
    DistNest_sc
  )

y <- maxent_data$presence

mx <- maxnet::maxnet(
  p = y,
  data = x,
  f = maxnet::maxnet.formula(
    p = y,
    data = x,
    classes = "lq"
  )
)

summary(mx)