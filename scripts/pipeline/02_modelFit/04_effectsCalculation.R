# ============================================================
# Create year-specific cross-tabulation of harvest effects
# from GAM coefficients
#
# Model:
# random ~ Logger_model +
#          sex * Year * harvesttimeSK24 * harvested50m
#
# Interpretation:
#   Negative = avoidance
#   Positive = preference
#   Near zero = neutral
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)

# ------------------------------------------------------------
# Extract coefficient table
# ------------------------------------------------------------

coef_table <- summary(gam_harvest)$p.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term")

# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------

get_coef <- function(term_name) {
  
  val <- coef_table %>%
    dplyr::filter(term == term_name) %>%
    dplyr::pull(Estimate)
  
  if (length(val) == 0) return(0)
  
  return(val)
}

# ============================================================
# Main coefficients
# ============================================================

H  <- get_coef("harvested50m1")

B  <- get_coef(
  "harvesttimeSK24beforeJuneHarvest:harvested50m1"
)

D  <- get_coef(
  "harvesttimeSK24duringJuneHarvest:harvested50m1"
)

M  <- get_coef(
  "sexMale:harvested50m1"
)

Y  <- get_coef(
  "Year2024:harvested50m1"
)

MB <- get_coef(
  "sexMale:harvesttimeSK24beforeJuneHarvest:harvested50m1"
)

MD <- get_coef(
  "sexMale:harvesttimeSK24duringJuneHarvest:harvested50m1"
)

YB <- get_coef(
  "Year2024:harvesttimeSK24beforeJuneHarvest:harvested50m1"
)

YD <- get_coef(
  "Year2024:harvesttimeSK24duringJuneHarvest:harvested50m1"
)

MY <- get_coef(
  "sexMale:Year2024:harvested50m1"
)

MYB <- get_coef(
  "sexMale:Year2024:harvesttimeSK24beforeJuneHarvest:harvested50m1"
)

MYD <- get_coef(
  "sexMale:Year2024:harvesttimeSK24duringJuneHarvest:harvested50m1"
)

# ============================================================
# Calculate net harvest effects
# ============================================================

effect_table <- tibble::tribble(
  ~Year, ~Sex, ~HarvestPeriod, ~Effect,
  
  # ----------------------------------------------------------
  # FEMALE 2023
  # ----------------------------------------------------------
  
  "2023", "Female", "Before harvest",
  H + B,
  
  "2023", "Female", "During harvest",
  H + D,
  
  "2023", "Female", "After harvest",
  H,
  
  # ----------------------------------------------------------
  # MALE 2023
  # ----------------------------------------------------------
  
  "2023", "Male", "Before harvest",
  H + B + M + MB,
  
  "2023", "Male", "During harvest",
  H + D + M + MD,
  
  "2023", "Male", "After harvest",
  H + M,
  
  # ----------------------------------------------------------
  # FEMALE 2024
  # ----------------------------------------------------------
  
  "2024", "Female", "Before harvest",
  H + B + Y + YB,
  
  "2024", "Female", "During harvest",
  H + D + Y + YD,
  
  "2024", "Female", "After harvest",
  H + Y,
  
  # ----------------------------------------------------------
  # MALE 2024
  # ----------------------------------------------------------
  
  "2024", "Male", "Before harvest",
  H + B + M + Y + MB + YB + MY + MYB,
  
  "2024", "Male", "During harvest",
  H + D + M + Y + MD + YD + MY + MYD,
  
  "2024", "Male", "After harvest",
  H + M + Y + MY
)

# ============================================================
# Biological interpretation
# ============================================================

effect_table <- effect_table %>%
  
  dplyr::mutate(
    
    Interpretation = dplyr::case_when(
      
      Effect < -0.3 ~ "Strong avoidance",
      
      Effect >= -0.3 &
        Effect < -0.1 ~ "Weak avoidance",
      
      Effect >= -0.1 &
        Effect <= 0.1 ~ "Neutral",
      
      Effect > 0.1 ~ "Preference"
    )
  )

# ============================================================
# Create cross-tabulated table
# ============================================================

effect_table_wide <- effect_table %>%
  
  dplyr::mutate(
    Effect = round(Effect, 3)
  ) %>%
  
  tidyr::pivot_wider(
    names_from = HarvestPeriod,
    values_from = c(
      Effect,
      Interpretation
    )
  )

# ============================================================
# Split by year
# ============================================================

table_2023 <- effect_table_wide %>%
  dplyr::filter(Year == "2023")

table_2024 <- effect_table_wide %>%
  dplyr::filter(Year == "2024")

# ============================================================
# Print tables
# ============================================================

cat("\n====================\n")
cat("2023 EFFECT TABLE\n")
cat("====================\n")

print(table_2023)

cat("\n====================\n")
cat("2024 EFFECT TABLE\n")
cat("====================\n")

print(table_2024)