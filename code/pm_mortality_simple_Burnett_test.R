# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### GEMM Mortality (Burnett et al. 2018 PNAS) ###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

library(tidyverse)
library(readxl)
library(janitor)
library(countrycode)
library(here)
library(purrr)
library(stringr)
library(readr)
library(dplyr)
library(tibble)
library(sf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################### 1) LOAD GRID DATA (PM, fire PM, population) #########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grid <- read_csv(
  here("output", "pop_pm_with_countries.csv"),
  show_col_types = FALSE
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################### 2) GEMM PARAMETERS θ, α, μ, ν  ######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gemm_xlsx <- here("input", "Burnett_et_al_2018", "GEMM Calculator (PNAS)_ab.xlsx")

gemm_params <- read_excel(
  gemm_xlsx,
  sheet = "GEMM fit parameters",
  range = "B11:G120",
  col_names = c("age_group","theta","se","alpha","mu","nu")
) %>%
  filter(age_group == "25+") %>%
  mutate(across(theta:nu, as.numeric)) %>%
  slice(1)

theta <- gemm_params$theta
alpha <- gemm_params$alpha
mu    <- gemm_params$mu
nu    <- gemm_params$nu

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################## 3) BUILD BASELINE MORTALITY Y_ct (FROM NCD+LRI) ########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ------------------------------------------------------------------------------
# Step 3.1 Read country population (denominator)
# ------------------------------------------------------------------------------

gbd_pop <- read_excel(
  gemm_xlsx,
  sheet = "GBD 2015 assessment",
  range = "B7:D400",
  col_names = c("country","pop_all","pop_25p")
) %>%
  mutate(
    pop_all = as.numeric(pop_all),
    pop_25p = as.numeric(pop_25p)
  )

# ------------------------------------------------------------------------------
# Step 3.2 Read NCD+LRI baseline deaths (numerator)
# ------------------------------------------------------------------------------

ncd_lri_raw <- read_excel(
  gemm_xlsx,
  sheet = "NCD+ LRI baseline",
  col_names = FALSE
)

# Convert to a character matrix for robust header detection
mat <- as.matrix(ncd_lri_raw)

clean_cell <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\u00A0", " ")  # non-breaking spaces
  x <- str_squish(x)                      # trim + collapse whitespace
  x
}

mat_chr <- apply(mat, c(1, 2), clean_cell)

# Find the header row containing the label "Location" anywhere
loc_hits <- which(mat_chr == "Location", arr.ind = TRUE)
if (nrow(loc_hits) == 0) stop("Could not find a header cell labeled 'Location'.")

header_row <- loc_hits[1, "row"]

# Read the header labels from that row
hdr <- mat_chr[header_row, ]

# Identify the key columns by header text (NOT by position)
location_col <- which(hdr == "Location")[1]
under_col    <- which(str_detect(hdr, "^under\\s*25$"))[1]
over_col     <- which(str_detect(hdr, "^over\\s*25$"))[1]

if (any(is.na(c(location_col, under_col, over_col)))) {
  stop("Could not locate one of the required header columns: Location / under 25 / over 25.")
}

# Slice data rows below the header
dat <- ncd_lri_raw[(header_row + 1):nrow(ncd_lri_raw), ]

# Extract country names from the Location column
country_vec <- clean_cell(dat[[location_col]])

# Keep real countries only (drop blanks + GLOBAL)
keep <- !is.na(country_vec) & country_vec != "" & country_vec != "GLOBAL"
dat <- dat[keep, ]
country_vec <- country_vec[keep]

# Parse numeric deaths
ncd_lri <- tibble(
  country = country_vec,
  deaths_under25 = parse_number(clean_cell(dat[[under_col]])),
  deaths_over25  = parse_number(clean_cell(dat[[over_col]]))
) %>%
  mutate(D_ct = deaths_under25 + deaths_over25)

# Quick sanity check: Afghanistan should be ~40842 and ~163765
ncd_lri %>% filter(country == "Afghanistan")

# ------------------------------------------------------------------------------
# Step 3.3 Construct baseline mortality rate Y_ct
# ------------------------------------------------------------------------------

Y <- gbd_pop %>%
  left_join(ncd_lri, by = "country") %>%
  mutate(
    Y_ct = D_ct / pop_all,
    iso3 = countrycode(country, "country.name", "iso3c")
  ) %>%
  select(iso3, Y_ct)

print(summary(Y$Y_ct))

# ------------------------------------------------------------------------------
# Step 3.4 Assign Y_ct to each grid cell
# ------------------------------------------------------------------------------

grid <- grid %>%
  left_join(Y, by = c("country_code_iso3" = "iso3"))

cat("Grid cells with valid Y_ct:", sum(!is.na(grid$Y_ct)), "\n")
cat("Grid cells missing Y_ct:", sum(is.na(grid$Y_ct)), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################## 4) GEMM ATTRIBUTABLE FRACTION (WITH cf UNCERTAINTY) ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(123)                           # for reproducibility
ndraw <- 300
cf_draws <- runif(ndraw, 2.4, 5.9)    # draw cf once globally

gemm_af_fast <- function(pm, cf = cf_draws) {
  if (is.na(pm)) return(NA_real_)
  z  <- pmax(0, pm - cf)
  g  <- log(1 + z / alpha) / (1 + exp((mu - z) / nu))
  rr <- exp(theta * g)
  mean((rr - 1) / rr)
}

# 1) Build AF lookup for unique PM values (huge speedup when PM repeats on a grid)
pm_vals <- grid$pm_2000
pm_unique <- sort(unique(pm_vals[!is.na(pm_vals)]))

af_lookup <- tibble(pm_2000 = pm_unique) %>%
  mutate(AF_2000 = map_dbl(pm_2000, gemm_af_fast))

# 2) Join AF back to grid and compute mortalities in ONE mutate
grid <- grid %>%
  left_join(af_lookup, by = "pm_2000") %>%
  mutate(
    pm_mort_2000  = pop_tot_2009 * Y_ct * AF_2000,
    # fire fraction (protect against divide-by-zero)
    fpm_mort_2000 = pm_mort_2000 * dplyr::if_else(pm_2000 > 0, fpm_2000 / pm_2000, NA_real_)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################### 7) USA TOTAL #######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

usa <- grid %>%
  filter(country_code_iso3 == "USA") %>%
  summarise(
    total_PM_deaths   = sum(pm_mort_2000, na.rm = TRUE),
    total_fire_deaths = sum(fpm_mort_2000, na.rm = TRUE)
  )

print(usa)

cat("\nFinished — GEMM mortality calculated\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################### 7) PLOT #######################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ============================================================
# Unified diagnostic maps (PM + Fire PM)
# >0 = green, 0 = red, NA = black
# ============================================================

library(dplyr)
library(leaflet)
library(scales)
library(htmltools)

# ------------------------------------------------------------
# 1. Detect lon/lat columns
# ------------------------------------------------------------

lon_candidates <- c("lon", "longitude", "x", "X", "Long", "LONGITUDE")
lat_candidates <- c("lat", "latitude", "y", "Y", "Lat", "LATITUDE")

lon_col <- intersect(lon_candidates, names(grid))[1]
lat_col <- intersect(lat_candidates, names(grid))[1]

if (is.na(lon_col) || is.na(lat_col)) {
  stop("Cannot find lon/lat columns in grid.")
}

cat("Longitude:", lon_col, "\n")
cat("Latitude:", lat_col, "\n")

# ------------------------------------------------------------
# 2. Aggregate grid to 1° cells (fast)
# ------------------------------------------------------------

cell_deg <- 1

grid_map <- grid %>%
  mutate(
    lon = .data[[lon_col]],
    lat = .data[[lat_col]]
  ) %>%
  filter(is.finite(lon), is.finite(lat)) %>%
  mutate(
    lon_bin = floor(lon / cell_deg) * cell_deg + cell_deg/2,
    lat_bin = floor(lat / cell_deg) * cell_deg + cell_deg/2
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarise(
    pm_mort_2000  = sum(pm_mort_2000,  na.rm = FALSE),
    fpm_mort_2000 = sum(fpm_mort_2000, na.rm = FALSE),
    pop_tot_2009  = sum(pop_tot_2009,  na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

cat("Number of map cells:", nrow(grid_map), "\n")

# ------------------------------------------------------------
# 3. Status classification (shared color logic)
# ------------------------------------------------------------

status_color <- function(x){
  ifelse(is.na(x), "yellow",
         ifelse(x == 0, "red", "green"))
}

# PM status
grid_map <- grid_map %>%
  mutate(
    pm_status = case_when(
      is.na(pm_mort_2000) ~ "NA",
      pm_mort_2000 == 0   ~ "0",
      pm_mort_2000 > 0    ~ ">0"
    ),
    pm_color = status_color(pm_mort_2000)
  )

# ------------------------------------------------------------
# 4. PM mortality map
# ------------------------------------------------------------

m_pm <- leaflet(grid_map) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(
    lng = ~lon_bin,
    lat = ~lat_bin,
    radius = 3,
    stroke = FALSE,
    fillOpacity = 0.8,
    color = ~pm_color,
    popup = ~HTML(paste0(
      "<b>PM status:</b> ", pm_status, "<br/>",
      "<b>PM deaths:</b> ", pm_mort_2000, "<br/>",
      "<b>Population:</b> ", comma(pop_tot_2009)
    ))
  ) %>%
  addLegend(
    position = "bottomright",
    colors = c("green","red","yellow"),
    labels = c(">0 (positive deaths)",
               "=0 (computed but none)",
               "NA (cannot compute)"),
    title = "PM mortality",
    opacity = 1
  )

m_pm

# ------------------------------------------------------------
# 5. Fire PM classification
# ------------------------------------------------------------

grid_map <- grid_map %>%
  mutate(
    fpm_status = case_when(
      is.na(fpm_mort_2000) ~ "NA",
      fpm_mort_2000 == 0   ~ "0",
      fpm_mort_2000 > 0    ~ ">0"
    ),
    fpm_color = status_color(fpm_mort_2000)
  )

# ------------------------------------------------------------
# 6. Fire PM mortality map
# ------------------------------------------------------------

m_fpm <- leaflet(grid_map) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(
    lng = ~lon_bin,
    lat = ~lat_bin,
    radius = 3,
    stroke = FALSE,
    fillOpacity = 0.8,
    color = ~fpm_color,
    popup = ~HTML(paste0(
      "<b>Fire PM status:</b> ", fpm_status, "<br/>",
      "<b>Fire deaths:</b> ", fpm_mort_2000, "<br/>",
      "<b>Population:</b> ", comma(pop_tot_2009)
    ))
  ) %>%
  addLegend(
    position = "bottomright",
    colors = c("green","red","yellow"),
    labels = c(">0 (positive deaths)",
               "=0 (computed but none)",
               "NA (cannot compute)"),
    title = "Fire PM mortality",
    opacity = 1
  )

m_fpm
