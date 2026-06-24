#' Finalise modelling dataset after harvest-overlap annotation
#'
#' Combine year-specific datasets, standardise column names created by joins,
#' derive binary "harvested" flags from harvest-timing categories, ensure
#' correct data types, and classify each observation relative to the June
#' harvest window.
#'
#' Inputs:
#'  - dmodel23, dmodel24 : sf objects (results of add_harvest_overlap or similar)
#'  - datetime_col       : name of the observation datetime column (POSIXct expected)
#'  - tz                 : timezone used when constructing harvest-period POSIXct
#'
#' Returns:
#'  - an sf object (combined and cleaned) ready for modelling
#'
finalise_dmodel <- function(dmodel23,
                            dmodel24,
                            datetime_col = "UTC_datetime",
                            tz = "UTC") {
  # ----------------------------
  # Namespace checks
  # ----------------------------
  pkgs <- c("dplyr", "sf", "tidyr")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Install missing packages: ", paste(missing_pkgs, collapse = ", "))

  `%>%` <- dplyr::`%>%`

  # ----------------------------
  # 1) Combine year datasets
  # ----------------------------
  dmodel <- dplyr::bind_rows(dmodel23, dmodel24)

  # ----------------------------
  # 2) Create binary harvested indicators from timing classes
  #    - after => "1"
  #    - before/noharvest => "0"
  #    - otherwise NA
  # ----------------------------
  if ("afterharvest50m" %in% names(dmodel)) {
    dmodel <- dmodel %>%
      dplyr::mutate(
        harvested50m = dplyr::case_when(
          afterharvest50m == "after" ~ "1",
          afterharvest50m %in% c("before", "noharvest") ~ "0",
          TRUE ~ NA_character_
        )
      )
  } else {
    dmodel$harvested50m <- NA_character_
  }

  if ("afterharvest200m" %in% names(dmodel)) {
    dmodel <- dmodel %>%
      dplyr::mutate(
        harvested200m = dplyr::case_when(
          afterharvest200m == "after" ~ "1",
          afterharvest200m %in% c("before", "noharvest") ~ "0",
          TRUE ~ NA_character_
        )
      )
  } else {
    dmodel$harvested200m <- NA_character_
  }

  # Convert harvested flags to factors with consistent levels if they exist
  if ("harvested50m" %in% names(dmodel)) dmodel$harvested50m <- factor(dmodel$harvested50m, levels = c("0", "1"))
  if ("harvested200m" %in% names(dmodel)) dmodel$harvested200m <- factor(dmodel$harvested200m, levels = c("0", "1"))

  # ----------------------------
  # 3) Standardise columns that may have been suffixed by joins
  #    Use coalesce to prefer unsuffixed names, fall back to .x or .y variants.
  # ----------------------------
  coalesce_if_present <- function(df, target, variants) {
    present_variants <- variants[variants %in% names(df)]
    if (length(present_variants) == 0) return(df)
    df[[target]] <- dplyr::coalesce(!!!rlang::syms(c(target, present_variants)))
    df
  }

  dmodel <- coalesce_if_present(dmodel, "Logger_model", c("Logger_model.x", "Logger_model.y"))
  dmodel <- coalesce_if_present(dmodel, "sex", c("sex.x", "sex.y"))
  dmodel <- coalesce_if_present(dmodel, datetime_col, c(paste0(datetime_col, ".x"), paste0(datetime_col, ".y")))
  dmodel <- coalesce_if_present(dmodel, "DistNest", c("DistNest.x", "DistNest.y"))
  dmodel <- coalesce_if_present(dmodel, "random", c("random.x", "random.y"))
  dmodel <- coalesce_if_present(dmodel, "Year", c("Year.x", "Year.y"))
  dmodel <- coalesce_if_present(dmodel, "ID", c("ID.x", "ID.y"))

  # ----------------------------
  # 4) Remove common '.x' / '.y' artefact columns created by joins
  # ----------------------------
  remove_suffix_cols <- function(df) {
    drop <- grep("\\.(x|y)$", names(df), value = TRUE)
    if (length(drop) > 0) df <- df %>% dplyr::select(-dplyr::any_of(drop))
    df
  }
  dmodel <- remove_suffix_cols(dmodel)

  # ----------------------------
  # 5) Ensure key modelling variables have correct classes
  # ----------------------------
  if ("Year" %in% names(dmodel)) dmodel$Year <- as.factor(as.character(dmodel$Year))

  # 'random' indicates whether point is background/used in some pipelines.
  # If not present, create it from 'used' if available.
  if (!"random" %in% names(dmodel) && "used" %in% names(dmodel)) {
    dmodel$random <- factor(dmodel$used)
  } else if ("random" %in% names(dmodel)) {
    dmodel$random <- as.factor(dmodel$random)
  }

  if ("afterharvest50m" %in% names(dmodel)) dmodel$afterharvest50m <- as.factor(dmodel$afterharvest50m)
  if ("afterharvest200m" %in% names(dmodel)) dmodel$afterharvest200m <- as.factor(dmodel$afterharvest200m)

  # ----------------------------
  # 6) Ensure datetime column is POSIXct
  # ----------------------------
  if (!datetime_col %in% names(dmodel)) {
    warning("Datetime column '", datetime_col, "' not found; harvest-time classification will be skipped.")
    return(dmodel)
  }
  if (!inherits(dmodel[[datetime_col]], "POSIXt")) {
    # Try to coerce; let NA result if parsing fails
    dmodel[[datetime_col]] <- as.POSIXct(dmodel[[datetime_col]], tz = tz)
  }

  # ----------------------------
  # 7) Classify each point relative to the June harvest window
  #    Uses fixed window (17 Jun 11:00 to 27 Jun 16:00) per year (as in original script).
  # ----------------------------
  # Build harvest start/end per-row using the Year column (fall back to observed year if Year missing)
  if (!"Year" %in% names(dmodel)) {
    # derive Year from datetime if missing
    dmodel$Year <- as.factor(format(dmodel[[datetime_col]], "%Y"))
  }

  # Create POSIXct boundaries based on Year
  dmodel <- dmodel %>%
    dplyr::mutate(
      .harvest_start = as.POSIXct(paste0(as.character(Year), "-06-17 11:00:00"), tz = tz),
      .harvest_end   = as.POSIXct(paste0(as.character(Year), "-06-27 16:00:00"), tz = tz),

      harvesttimeSK24 = dplyr::case_when(
        .data[[datetime_col]] < .harvest_start ~ "beforeJuneHarvest",
        .data[[datetime_col]] >= .harvest_start & .data[[datetime_col]] <= .harvest_end ~ "duringJuneHarvest",
        .data[[datetime_col]] > .harvest_end ~ "afterJuneHarvest",
        TRUE ~ NA_character_
      ),
      harvesttimeSK24 = as.factor(harvesttimeSK24)
    ) %>%
    dplyr::select(-.harvest_start, -.harvest_end)

  # ----------------------------
  # 8) Return final dataset
  # ----------------------------
  dmodel
}
