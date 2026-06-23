eda <- run_eda_by_year(
  dmodel = dmodel,
  years = c("2023", "2024"),
  response_col = "used"
)


eda[["2023"]]$missing_values
eda[["2023"]]$frequency_tables
eda[["2023"]]$cross_tables
eda[["2023"]]$plots$distnest_hist
eda[["2023"]]$plots$pairs

eda[["2024"]]$plots$harvest50
eda[["2024"]]$plots$pairs
eda[["2024"]]$missing_values
eda[["2024"]]$frequency_tables
eda[["2024"]]$cross_tables
eda[["2024"]]$plots$distnest_hist
