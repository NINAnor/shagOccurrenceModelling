for (i in 1:12) {
tryCatch({

  options(warn = 2)
  on.exit(options(warn = 0), add = TRUE)

  assign("param_row", i, envir = .GlobalEnv)
  source("scripts/masterScript.R")

}, error = function(e) {
  message("Model ", i, ": ", conditionMessage(e))
})

}
