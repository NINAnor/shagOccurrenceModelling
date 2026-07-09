run_model <- function(param_row) {

  assign(
    "param_row",
    param_row,
    envir = .GlobalEnv
  )

  source("scripts/masterScript.R")
}

for(i in 1:12){
  run_model(i)
}
