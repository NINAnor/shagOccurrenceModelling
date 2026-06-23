# ============================================================
# Finalise model dataset after harvest overlap processing
#
# Inputs:
#   dmodel23, dmodel24 = sf objects after add_harvest_overlap()
#
# Returns:
#   Combined and cleaned sf object ready for modelling
# ============================================================

finalise_dmodel <- function(dmodel23,
                            dmodel24,
                            datetime_col = "UTC_datetime") {
  
  # ----------------------------------------------------------
  # Combine years
  # ----------------------------------------------------------
  
  dmodel <- dplyr::bind_rows(dmodel23, dmodel24)
  
  
  # ----------------------------------------------------------
  # Convert harvest timing classes to binary harvested variables
  #
  # after     = 1
  # before    = 0
  # noharvest = 0
  # ----------------------------------------------------------
  
  dmodel <- dmodel %>%
    dplyr::mutate(
      harvested50m = dplyr::case_when(
        afterharvest50m == "after" ~ "1",
        afterharvest50m %in% c("before", "noharvest") ~ "0",
        TRUE ~ NA_character_
      ),
      
      harvested200m = dplyr::case_when(
        afterharvest200m == "after" ~ "1",
        afterharvest200m %in% c("before", "noharvest") ~ "0",
        TRUE ~ NA_character_
      ),
      
      harvested50m = as.factor(harvested50m),
      harvested200m = as.factor(harvested200m)
    )
  
  
  # ----------------------------------------------------------
  # Standardise column names
  #
  # This handles cases where previous joins created .x names.
  # If clean names already exist, they are left unchanged.
  # ----------------------------------------------------------
  
  if ("Logger_model.x" %in% names(dmodel)) {
    dmodel$Logger_model <- dmodel$Logger_model.x
  }
  
  if ("sex.x" %in% names(dmodel)) {
    dmodel$sex <- dmodel$sex.x
  }
  
  if ("UTC_datetime.x" %in% names(dmodel)) {
    dmodel$UTC_datetime <- dmodel$UTC_datetime.x
  }
  
  if ("DistNest.x" %in% names(dmodel)) {
    dmodel$DistNest <- dmodel$DistNest.x
  }
  
  if ("random.x" %in% names(dmodel)) {
    dmodel$random <- dmodel$random.x
  }
  
  if ("Year.x" %in% names(dmodel)) {
    dmodel$Year <- dmodel$Year.x
  }
  
  if ("ID.x" %in% names(dmodel)) {
    dmodel$ID <- dmodel$ID.x
  }
  
  
  # ----------------------------------------------------------
  # Remove redundant .x columns if present
  # ----------------------------------------------------------
  
  dmodel <- dmodel %>%
    dplyr::select(
      -dplyr::any_of(c(
        "Year.x",
        "sex.x",
        "Logger_model.x",
        "UTC_datetime.x",
        "random.x",
        "DistNest.x",
        "ID.x"
      ))
    )
  
  
  # ----------------------------------------------------------
  # Ensure key modelling variables have correct classes
  # ----------------------------------------------------------
  
  dmodel <- dmodel %>%
    dplyr::mutate(
      Year = as.factor(Year),
      random = as.factor(used),
      afterharvest50m = as.factor(afterharvest50m),
      afterharvest200m = as.factor(afterharvest200m)
    )
  
  
  # ----------------------------------------------------------
  # Classify each point relative to June harvest period
  #
  # Uses year-specific dates:
  #   2023: 17 June to 27 June
  #   2024: 17 June to 27 June
  # ----------------------------------------------------------
  
  dmodel <- dmodel %>%
    dplyr::mutate(
      harvest_start = as.POSIXct(
        paste0(Year, "-06-17 11:00:00"),
        tz = "UTC"
      ),
      harvest_end = as.POSIXct(
        paste0(Year, "-06-27 16:00:00"),
        tz = "UTC"
      ),
      
      harvesttimeSK24 = dplyr::case_when(
        .data[[datetime_col]] < harvest_start ~ "beforeJuneHarvest",
        .data[[datetime_col]] >= harvest_start &
          .data[[datetime_col]] <= harvest_end ~ "duringJuneHarvest",
        .data[[datetime_col]] > harvest_end ~ "afterJuneHarvest",
        TRUE ~ NA_character_
      ),
      
      harvesttimeSK24 = as.factor(harvesttimeSK24)
    ) %>%
    dplyr::select(
      -harvest_start,
      -harvest_end
    )
  
  
  # ----------------------------------------------------------
  # Return final modelling dataset
  # ----------------------------------------------------------
  
  return(dmodel)
}