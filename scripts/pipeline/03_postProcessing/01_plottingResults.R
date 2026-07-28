library(readxl)
library(dplyr)
library(tidyr)
library(scales)
library(purrr)
library(ggplot2)
library(ggrepel)

file <- "scripts/data/modelOutputs.xlsx"

# Get all sheet names
sheet_names <- excel_sheets(file)

# Read and combine all sheets
modelOutputs <- map_dfr(sheet_names, function(sheet) {

  read_excel(file, sheet = sheet) %>%
    fill(Year, Gender, .direction = "down")%>%
    dplyr::mutate(Mean  = as.numeric(Mean),
                  Lower = as.numeric(Lower),
                  Upper = as.numeric(Upper),
                  Year = factor(Year),
                  Gender = factor(Gender),
                  Response = if_else(Mean > 0, "Preference", "Avoidance"),
                  Significant = factor(Significant, levels = c("YES", "NO")),
                  HarvestRegime = factor(
                    `Harvest Regime`,
                    levels = c(
                      "Before June Harvest",
                      "During June Harvest",
                      "After June Harvest"
                    )
                  )
    ) %>%
    mutate(Model = sheet)

})


pd <- position_dodge(width = 0.6)

p <- ggplot(modelOutputs,
       aes(x = Mean,
           y = HarvestRegime,
           colour = Response,
           shape = Significant,
           group = Model)) +

  geom_vline(xintercept = 0,
             linetype = 2,
             colour = "grey40") +

  geom_errorbarh(
    aes(xmin = Lower,
        xmax = Upper,
        linetype = Model),
    height = 0.15,
    linewidth = 0.8,
    position = pd
  ) +

  geom_point(
    aes(fill = Model),
    size = 3.5,
    position = pd
  ) +

  facet_grid(Gender ~ Year) +

  scale_colour_manual(values = c(
    Preference = "#1f78b4",
    Avoidance = "#d73027"
  )) +

  scale_shape_manual(values = c(
    YES = 16,
    NO = 1
  )) +

  scale_fill_brewer(palette = "Dark2") +

  scale_x_continuous(
    trans = pseudo_log_trans(base = 10)
  ) +

  labs(
    x = "Selection coefficient",
    y = NULL,
    fill = "Model",
    colour = "Response",
    shape = "Significant",
    title = "Shag response to kelp harvesting"
  ) +

  theme_bw(base_size = 13) +

  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    axis.text.y = element_text(size = 11)
  )

output_dir <- "R:/Prosjekter/Toppskarv/Analyses_Nina_KelpHarvest/kwaku/results/version_0.1.2/publication_results"

# PDF (vector graphics; ideal for journals)
ggsave(
  filename = file.path(output_dir, "shag_harvest_selection.pdf"),
  plot = p,
  width = 15,
  height = 10,
  units = "in",
  device = cairo_pdf
)

# PNG (high-resolution)
ggsave(
  filename = file.path(output_dir, "shag_harvest_selection.png"),
  plot = p,
  width = 15,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
