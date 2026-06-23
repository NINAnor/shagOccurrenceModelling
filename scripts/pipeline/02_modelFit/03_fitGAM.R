model_data <- dmodel %>%
  
  sf::st_drop_geometry() %>%
  
  dplyr::mutate(
    
    random = as.numeric(as.character(random)),
    
    sex = factor(sex),
    
    Logger_model = factor(Logger_model),
    
    Year = factor(Year),
    
    harvested50m = factor(harvested50m),
    
    harvesttimeSK24 = factor(harvesttimeSK24),
    ID = factor(ID),
    
    DistNest_sc = as.numeric(scale(DistNest))
  )


gam_harvest <- mgcv::gam(
  
  random ~
    
    #sex +
    
   # Logger_model +
    

   sex * Year * harvesttimeSK24 * harvested50m +
    
    s(DistNest_sc) +
    +
    
    s(ID, bs = "re"),
  
  data = model_data,
  
  family = binomial,
  
  method = "REML"
)

summary(gam_harvest)

library(ggeffects)

pred <- ggeffects::ggpredict(
  gam_harvest,
  terms = c(
    "harvested50m",
    "harvesttimeSK24",
    "sex",
    "Year"
  )
)

plot(pred)
