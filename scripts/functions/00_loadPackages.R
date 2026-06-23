# This script load all the packages needed to fit the models

# Get a list of installed packages
installedPackages <- installed.packages()


# Now we bring int he rest of the necessary packages
necessaryPackages <- c("sf",
                       "terra",
                       "dplyr",
                       "lubridate",
                       "data.table",
                      # "chillR",
                       #"diveMove",
                       "geosphere",
                       "ggplot2",
                       "adehabitatHR",
                       "ggspatial",
                      # "gamm4",
                       "pROC",
                       "car",
                       "ggiraphExtra",
                       "ggeffects",
                       "patchwork",
                       "tidyr",
                       "rlang",
                       "purrr",
                       "GGally",
                       "arrow",
                       "biomod2",
                       "cito",
                       "mda",
                       "earth",
                       "smoothr",
                       "gam"
)
uninstalledPackages <- necessaryPackages[!(necessaryPackages %in% row.names(installedPackages))]


load_required_packages <- function(pkgs) {
  
  missing_pkgs <- pkgs[
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nInstall them once manually before running parallel jobs:\n",
      "install.packages(c(",
      paste(sprintf('\"%s\"', missing_pkgs), collapse = ", "),
      "))"
    )
  }
  
  invisible(
    lapply(
      pkgs,
      library,
      character.only = TRUE
    )
  )
}
#install.packages(uninstalledPackages)
load_required_packages(necessaryPackages)

# library(sf)
# library(terra)
# library(dplyr)
# library(lubridate)
# library(data.table)
# library(chillR)
# library(diveMove)
# library(geosphere)
# library(ggplot2)
# library(sp)
# library(adehabitatHR)
# library(ggspatial)
# library(purrr)
# library(rlang)
# library(GGally)
# library(tidyr)
# library(arrow)
# library(smoothr)
# 
# # for modelling
# 
# library(gamm4)
# library(mgcv)
# library(pROC)
# library(sf)
# library(ggplot2)
# library(car)
# library(ggiraphExtra)
# library(ggeffects)
# library(patchwork)
# library(tidyr)
# library(biomod2)

