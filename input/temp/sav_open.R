library(readr)
library(dplyr)
library(leaflet)
library(htmlwidgets)
library(here)
library(tidyr)

setwd("C:/Ford_BA_FPM25")
file.create(".here")

df <- read_csv(
  here("input", "temp", "regionmask_full.csv"),
  show_col_types = FALSE
)


# Keep only cells belonging to a US state
us <- df %>%
  filter(state != "NotState")

nrow(us)

#distribution of column
summary_dist <-
  df %>%
  mutate(
    epa_region  = as.character(epa_region),
    nca4_region = as.character(nca4_region),
    mvm_region  = as.character(mvm_region)
  ) %>%
  pivot_longer(
    cols = c(state, epa_region, nca4_region, mvm_region),
    names_to = "variable",
    values_to = "category"
  ) %>%
  count(variable, category) %>%
  arrange(variable, desc(n))

summary_dist

#plot
pal <- colorFactor(
  palette = "Set3",
  domain = us$state
)

m <- leaflet(us) %>%
  
  # global background
  addProviderTiles(providers$CartoDB.Positron) %>%
  setView(lng = 0, lat = 25, zoom = 2) %>%
  
  #points with values
  addCircleMarkers(
    lng = ~lon,
    lat = ~lat,
    radius = 3,
    stroke = FALSE,
    fillOpacity = 0.7,
    color = ~pal(state),
    
    popup = ~paste0(
      "<b>State:</b> ", state, "<br>",
      "<b>Lat:</b> ", round(lat,3), "<br>",
      "<b>Lon:</b> ", round(lon,3)
    )
  ) %>%
  
  addLegend(
    "bottomright",
    pal = pal,
    values = us$state,
    title = "US State",
    opacity = 1
  )

m

