# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################GEMM Mortality (Burnett et al. 2018 PNAS)##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

library(tidyverse)
library(readxl)
library(here)
library(countrycode)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #########################1) Load grid##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid <- read_csv(here("output", "pop_pm_country_death.csv"), show_col_types = FALSE)

years <- 2001:2010 # select years

stopifnot(all(paste0("death_rate_", years) %in% names(grid)))
stopifnot(all(paste0("pop_tot_",   years) %in% names(grid)))
stopifnot(all(c("pm_2000","fpm_2000","country_code_iso3") %in% names(grid)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################2) Load GEMM parameters (25+)##########################
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
# ) ##########################3) GEMM: compute AF and diagnostic RR (mean RR) for a cell##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gemm_af_rr_cell <- function(pm, cf_vec, theta, alpha, mu, nu) {
  if (is.na(pm)) return(c(AF = NA_real_, RR = NA_real_))
  z  <- pmax(0, pm - cf_vec)
  g  <- log(1 + z / alpha) / (1 + exp((mu - z) / nu))
  rr <- exp(theta * g)
  c(
    AF = mean((rr - 1) / rr),
    RR = mean(rr)
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################4) Cell-specific counterfactual draws + AF_2000 + RR_2000###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)
ndraw  <- 300
cf_min <- 2.4
cf_max <- 5.9

cf_mat <- matrix(
  runif(nrow(grid) * ndraw, cf_min, cf_max),
  nrow = nrow(grid),
  ncol = ndraw
)

af_rr_list <- purrr::map(
  seq_len(nrow(grid)),
  \(i) gemm_af_rr_cell(grid$pm_2000[i], cf_mat[i, ], theta, alpha, mu, nu)
)

af_rr_df <- bind_rows(lapply(af_rr_list, as.list))
grid$AF_2000 <- af_rr_df$AF
grid$RR_2000 <- af_rr_df$RR

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################5) Mortality 2001–2010 (pm_2000 fixed)##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (yr in years) {
  pop_col <- paste0("pop_tot_", yr)
  dr_col  <- paste0("death_rate_", yr)
  pm_out  <- paste0("pm_mort_", yr)
  fpm_out <- paste0("fpm_mort_", yr)
  
  grid[[pm_out]] <- grid[[pop_col]] * grid[[dr_col]] * grid$AF_2000
  grid[[fpm_out]] <- grid[[pm_out]] * if_else(grid$pm_2000 > 0, grid$fpm_2000 / grid$pm_2000, NA_real_)
}

# rename final dataset
Burnett_mortality <- grid
rm(grid)

#=================================== save results===================================
#write_csv(Burnett_mortality, here("output", "Burnett_mortality.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################6a) USA totals (2001–2010 sum)##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
usa_total <- Burnett_mortality %>%
  filter(country_code_iso3 == "USA") %>%
  summarise(
    total_PM_deaths_2001_2010 = sum(rowSums(across(all_of(paste0("pm_mort_",  years))),  na.rm = TRUE), na.rm = TRUE),
    total_fire_PM_deaths_2001_2010 = sum(rowSums(across(all_of(paste0("fpm_mort_", years))), na.rm = TRUE), na.rm = TRUE)
  )

print(usa_total)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################6b) USA totals BY YEAR##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
usa_yearly <- Burnett_mortality %>%
  filter(country_code_iso3 == "USA") %>%
  summarise(
    across(all_of(paste0("pm_mort_",  years)),  ~sum(.x, na.rm = TRUE)),
    across(all_of(paste0("fpm_mort_", years)), ~sum(.x, na.rm = TRUE))
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "deaths"
  ) %>%
  mutate(
    type = if_else(str_detect(variable, "fpm"), "Fire PM", "Total PM"),
    year = as.integer(str_extract(variable, "\\d{4}"))
  ) %>%
  select(year, type, deaths) %>%
  arrange(year, type)

print(usa_yearly)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################10) GLOBAL totals BY YEAR##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

global_yearly <- Burnett_mortality %>%
  filter(!is.na(country_code_iso3), country_code_iso3 != "-99") %>%
  summarise(
    across(all_of(paste0("pm_mort_",  years)),  ~sum(.x, na.rm = TRUE)),
    across(all_of(paste0("fpm_mort_", years)), ~sum(.x, na.rm = TRUE))
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "deaths"
  ) %>%
  mutate(
    type = if_else(str_detect(variable, "fpm"), "Fire PM", "Total PM"),
    year = as.integer(str_extract(variable, "\\d{4}"))
  ) %>%
  select(year, type, deaths) %>%
  arrange(year, type)

cat("\n================ GLOBAL YEARLY DEATHS (2001–2010) ================\n")
print(global_yearly %>% mutate(deaths = round(deaths, 0)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################7) Is RR identical within each country?##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rr_check <- Burnett_mortality %>%
  filter(!is.na(RR_2000), !is.na(country_code_iso3)) %>%
  group_by(country_code_iso3) %>%
  summarise(
    n_cells = n(),
    rr_min  = min(RR_2000),
    rr_max  = max(RR_2000),
    rr_same_within_country = (rr_max - rr_min) < 1e-12,
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

print(rr_check)

cat("\nCounts (RR identical within country?):\n")
print(rr_check %>% count(rr_same_within_country))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################8) Are RR the same ACROSS countries?##########################
#     (using country-level mean RR)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rr_country_mean <- Burnett_mortality %>%
  filter(!is.na(RR_2000), !is.na(country_code_iso3)) %>%
  group_by(country_code_iso3) %>%
  summarise(
    n_cells = n(),
    rr_mean = mean(RR_2000),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

print(rr_country_mean)

# TRUE if all countries have (numerically) the same mean RR
all_equal_across_countries <-
  (max(rr_country_mean$rr_mean) - min(rr_country_mean$rr_mean)) < 1e-12

cat("\nAre country mean RR values identical across countries?\n")
print(all_equal_across_countries)

cat("\nRange of country mean RR (min, max):\n")
print(range(rr_country_mean$rr_mean, na.rm = TRUE))

cat("\nNumber of unique country mean RR values (rounded to 8 decimals):\n")
print(n_distinct(round(rr_country_mean$rr_mean, 8)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##########################9) Overall mortality per country (ALL NUMERIC OUTPUT)##########################
#     - deaths by country × year (2001–2010): total PM + fire PM
#     - rates per 100,000 (numeric)
#     - plus 2001–2010 totals per country (numeric)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
options(scipen = 999)         # avoid scientific notation in printing

# 9a) country × year totals (deaths + rates)
country_mortality_long <- Burnett_mortality %>%
  filter(!is.na(country_code_iso3), country_code_iso3 != "-99") %>%
  group_by(country_code_iso3) %>%
  summarise(
    across(all_of(paste0("pop_tot_", years)),  ~sum(.x, na.rm = TRUE), .names = "{.col}"),
    across(all_of(paste0("pm_mort_",  years)), ~sum(.x, na.rm = TRUE), .names = "{.col}"),
    across(all_of(paste0("fpm_mort_", years)), ~sum(.x, na.rm = TRUE), .names = "{.col}"),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = matches("^(pop_tot_|pm_mort_|fpm_mort_)"),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(str_extract(var, "\\d{4}")),
    metric = case_when(
      str_detect(var, "^pop_tot_")  ~ "pop",
      str_detect(var, "^pm_mort_")  ~ "pm_deaths",
      str_detect(var, "^fpm_mort_") ~ "fire_deaths"
    )
  ) %>%
  select(country_code_iso3, year, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(
    pm_rate_per100k   = (pm_deaths   / pop) * 1e5,
    fire_rate_per100k = (fire_deaths / pop) * 1e5,
    country_name      = countrycode(country_code_iso3, "iso3c", "country.name")
  ) %>%
  select(country_code_iso3, country_name, year,
         pop, pm_deaths, fire_deaths, pm_rate_per100k, fire_rate_per100k) %>%
  arrange(country_code_iso3, year)

print(country_mortality_long %>% mutate(across(where(is.numeric), ~round(.x, 2))))


### THE END 

