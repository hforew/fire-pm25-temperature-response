# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FPM–GMT RELATIONSHIP COMPARISON: Compare beta^(c) (per-capita fire PM2.5 change per 1°C GMT)
##                        from multiple regression specifications
##
## Goal: For each country c and regressions specification, compare betas from the following
##    pooled FE regression across Pierce/Park/Zhao data
##    regression with Pierce data
##    regression with Park data
## Inputs:
##   fpm_gmt_regression_coefs_classic.csv
##   fpm_gmt_regression_coefs_jules.csv
##   fpm_gmt_regression_coefs_ssib4.csv
##   fpm_gmt_regression_coefs_classic_cli.csv
##   fpm_gmt_regression_coefs_jules_cli.csv
##   fpm_gmt_regression_coefs_ssib4_cli.csv
##   fpm_gmt_regression_coefs_FE_all.csv
##   fpm_gmt_regression_coefs_pierce.csv
## Outputs: 
##    beta_select.tex
##    beta_summary.tex
## Execution Order:
##   Execute after: fpm_gmt_regression_park.R, fpm_gmt_regression_pierce.R, fpm_gmt_regression_FE_all.R
##   Execute before: NA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment to start fresh
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)       # for portable file paths relative to project root
library(tidyverse)  # for data manipulation (dplyr, tidyr, purrr) and ggplot2
library(tibble)     # for tidy data frames (part of tidyverse, loaded explicitly for clarity)
library(xtable)     # for LaTeX table output

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# betas from Park et al models (CLASSIC, JULES, SSiB4)
beta_classic <- read_csv(here("output", "fpm_gmt_regression_coefs_classic.csv"))
beta_jules   <- read_csv(here("output", "fpm_gmt_regression_coefs_jules.csv"))
beta_ssib4   <- read_csv(here("output", "fpm_gmt_regression_coefs_ssib4.csv"))

# betas from Park et al climate-attributable regressions (CLASSIC, JULES, SSiB4)
beta_classic_cli <- read_csv(here("output", "fpm_gmt_regression_coefs_classic_cli.csv"))
beta_jules_cli   <- read_csv(here("output", "fpm_gmt_regression_coefs_jules_cli.csv"))
beta_ssib4_cli   <- read_csv(here("output", "fpm_gmt_regression_coefs_ssib4_cli.csv"))

head(beta_classic)
head(beta_jules)
head(beta_ssib4)

# beta from FE regression (Park/Pierce/Zhao data)
beta_fe     <- read_csv(here("output", "fpm_gmt_regression_coefs_FE_all.csv"))
beta_fe_cli <- read_csv(here("output", "fpm_gmt_regression_coefs_FE_all_cli.csv"))
head(beta_fe)
head(beta_fe_cli)

# beta from Pirce et al
beta_pierce <- read_csv(here("output", "fpm_gmt_regression_coefs_pierce.csv"))
head(beta_pierce)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Join and Compare betas ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# one row per country; beta_* = estimate, p_* = p-value, se_* = std error for each spec
beta_comparison <- beta_fe |>
  select(country_code_iso3, country_name, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
  rename(beta_fe = estimate_beta_c, p_fe = p.value_beta_c, se_fe = std.error_beta_c) |>          # FE (all data)
  left_join(
    beta_pierce |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
      rename(beta_cesm = estimate_beta_c, p_cesm = p.value_beta_c, se_cesm = std.error_beta_c), # Pierce et al
    by = "country_code_iso3"
  ) |>
  left_join(
    beta_classic |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
      rename(beta_classic = estimate_beta_c, p_classic = p.value_beta_c, se_classic = std.error_beta_c), # Park CLASSIC
    by = "country_code_iso3"
  ) |>
  left_join(
    beta_jules |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
      rename(beta_jules = estimate_beta_c, p_jules = p.value_beta_c, se_jules = std.error_beta_c),   # Park JULES
    by = "country_code_iso3"
  ) |>
  left_join(
    beta_ssib4 |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
      rename(beta_ssib4 = estimate_beta_c, p_ssib4 = p.value_beta_c, se_ssib4 = std.error_beta_c),   # Park SSiB4
    by = "country_code_iso3"
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))                   # round all numeric to 4dp

# one row per country; climate-attributable (cli) Park betas alongside FE and CESM benchmarks
beta_comparison_cli <- beta_fe_cli |>
  select(country_code_iso3, country_name, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
  rename(beta_fe = estimate_beta_c, p_fe = p.value_beta_c, se_fe = std.error_beta_c) |>          # FE cli (all data)
  left_join(
    beta_pierce |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
      rename(beta_cesm = estimate_beta_c, p_cesm = p.value_beta_c, se_cesm = std.error_beta_c),    # Pierce et al
    by = "country_code_iso3"
  ) |>
  left_join(
    beta_classic_cli |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
      rename(beta_classic_cli = estimate_beta_c, p_classic_cli = p.value_beta_c, se_classic_cli = std.error_beta_c), # Park CLASSIC cli
    by = "country_code_iso3"
  ) |>
  left_join(
    beta_jules_cli |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
      rename(beta_jules_cli = estimate_beta_c, p_jules_cli = p.value_beta_c, se_jules_cli = std.error_beta_c),     # Park JULES cli
    by = "country_code_iso3"
  ) |>
  left_join(
    beta_ssib4_cli |>
      select(country_code_iso3, estimate_beta_c, p.value_beta_c, std.error_beta_c) |>
      rename(beta_ssib4_cli = estimate_beta_c, p_ssib4_cli = p.value_beta_c, se_ssib4_cli = std.error_beta_c),     # Park SSiB4 cli
    by = "country_code_iso3"
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))                   # round all numeric to 4dp


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Compare select betas ##############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 11 countries of interest; order preserved via factor levels
countries <- c("global","CHN", "IND", "USA", "IDN", "PAK", "NGA", "BRA", "BDG", "RUS", "MEX")

# total fire PM betas (obsclim fire PM minus no-fire baseline)
beta_comparison_select <- beta_comparison |>
  filter(country_code_iso3 %in% countries) |>
  select(-country_name, -starts_with("p_")) |>                          # drop name and p-values; use SE
  mutate(country_code_iso3 = factor(country_code_iso3, levels = countries)) |> # convert to factor so arrange() sorts by countries order, not alphabetically
  arrange(country_code_iso3)

beta_comparison_select

# climate-attributable fire PM betas (obsclim minus counterclim)
beta_comparison_select_cli <- beta_comparison_cli |>
  filter(country_code_iso3 %in% countries) |>
  select(-country_name, -starts_with("p_")) |>                          # drop name and p-values; use SE
  mutate(country_code_iso3 = factor(country_code_iso3, levels = countries)) |>
  arrange(country_code_iso3)

beta_comparison_select
beta_comparison_select_cli

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Summary statistics by model #############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# total fire PM betas (obsclim fire PM minus no-fire baseline)
beta_summary <- beta_comparison |>
  select(country_code_iso3, starts_with("beta_"), starts_with("p_")) |>
  pivot_longer(
    cols = -country_code_iso3,
    names_to  = c(".value", "model"),
    names_pattern = "^(beta|p)_(.*)"     # splits e.g. "beta_fe" --> .value="beta", model="fe"
  ) |>
  group_by(model) |>
  summarise(
    pct_beta_pos  = round(mean(beta > 0,      na.rm = TRUE) * 100, 2), # % countries with beta > 0
    pct_beta_sig5 = round(mean(p    < 0.05,   na.rm = TRUE) * 100, 2), # % sig. at 5%
    beta_min      = round(min(beta,            na.rm = TRUE),       2),
    beta_p25      = round(quantile(beta, 0.25, na.rm = TRUE),       2),
    beta_median   = round(median(beta,         na.rm = TRUE),       2),
    beta_p75      = round(quantile(beta, 0.75, na.rm = TRUE),       2),
    beta_max      = round(max(beta,            na.rm = TRUE),       2)
  ) |>
  mutate(model = recode(model,
    fe      = "All models",
    cesm    = "CESM",
    classic = "CLASSIC",
    jules   = "JULES",
    ssib4   = "SSiB4"
  )) |>
  mutate(model = factor(model,
    levels = c("All models", "CESM", "CLASSIC", "JULES", "SSiB4"))) |>
  arrange(model)

beta_summary

# climate-attributable fire PM betas (obsclim minus counterclim)
beta_summary_cli <- beta_comparison_cli |>
  select(country_code_iso3, starts_with("beta_"), starts_with("p_")) |>
  pivot_longer(
    cols = -country_code_iso3,
    names_to  = c(".value", "model"),
    names_pattern = "^(beta|p)_(.*)"     # splits e.g. "beta_classic_cli" --> .value="beta", model="classic_cli"
  ) |>
  group_by(model) |>
  summarise(
    pct_beta_pos  = round(mean(beta > 0,      na.rm = TRUE) * 100, 2),
    pct_beta_sig5 = round(mean(p    < 0.05,   na.rm = TRUE) * 100, 2),
    beta_min      = round(min(beta,            na.rm = TRUE),       2),
    beta_p25      = round(quantile(beta, 0.25, na.rm = TRUE),       2),
    beta_median   = round(median(beta,         na.rm = TRUE),       2),
    beta_p75      = round(quantile(beta, 0.75, na.rm = TRUE),       2),
    beta_max      = round(max(beta,            na.rm = TRUE),       2)
  ) |>
  mutate(model = recode(model,
    fe          = "All models",
    cesm        = "CESM",
    classic_cli = "CLASSIC",
    jules_cli   = "JULES",
    ssib4_cli   = "SSiB4"
  )) |>
  mutate(model = factor(model,
    levels = c("All models", "CESM", "CLASSIC", "JULES", "SSiB4"))) |>
  arrange(model)

beta_summary
beta_summary_cli

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ LaTeX Tables #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LaTeX math column headers; rendered correctly via sanitize.colnames.function = identity
colnames(beta_summary) <- c("Model",
  "\\% $\\beta_c > 0$", "\\% $p < .05$",
  "Min", "P25", "Median", "P75", "Max")

print(
  xtable(
    beta_summary,
    caption = paste0("Summary statistics of country-level $\\beta_c$ estimates",
                     " across regression specifications."),
    label   = "tab:beta_summary",                         # \ref{tab:beta_summary} in LaTeX
    # length 9: 1 mandatory row-name slot + 8 data cols; 0dp for % cols, 2dp for beta cols
    digits  = c(0, 0, 0, 0, 2, 2, 2, 2, 2)
  ),
  booktabs                   = TRUE,                      # \toprule/\midrule/\bottomrule rules
  caption.placement          = "top",                     # caption above table
  include.rownames           = FALSE,                     # suppress R row indices (1,2,3...)
  sanitize.colnames.function = identity,                  # pass col names to LaTeX as-is (no escaping)
  sanitize.text.function     = identity,                  # pass cell values to LaTeX as-is (no escaping)
  file = here("output", "latex", "beta_summary.tex")
)

# Double-header table for beta_comparison_select.
# xtable does not natively support spanning headers, so the two-row header is
# injected as raw LaTeX via add.to.row; include.colnames = FALSE suppresses
# the default single-row header.
latex_beta_select <- beta_comparison_select |>
  mutate(country_code_iso3 = as.character(country_code_iso3))  # drop factor for xtable

# Row 1: model group labels spanning 2 cols each (beta + p)
# Row 2: Country | beta | p | beta | p | ...
# \cmidrule(lr){} draws partial rules under each group in row 1
header_row <- paste0(
  " & \\multicolumn{2}{c}{FE} & \\multicolumn{2}{c}{CESM} & ",
  "\\multicolumn{2}{c}{CLASSIC} & \\multicolumn{2}{c}{JULES} & ",
  "\\multicolumn{2}{c}{SSiB4} \\\\\n",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}\n",
  "Country & $\\beta$ & SE & $\\beta$ & SE & $\\beta$ & SE & ",
  "$\\beta$ & SE & $\\beta$ & SE \\\\\n",
  "\\midrule\n"
)

print(
  xtable(
    latex_beta_select,
    caption = paste0("$\\beta_c$ estimates and p-values across",
                     " specifications for largest countries."),
    label  = "tab:beta_select",                              # \ref{tab:beta_select} in LaTeX
    # length 12: 1 row-name slot + 11 data cols; 0dp for country, 3dp for all numeric
    digits = c(0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
  ),
  booktabs                   = TRUE,                         # \toprule/\midrule/\bottomrule rules
  caption.placement          = "top",                        # caption above table
  include.rownames           = FALSE,                        # suppress R row indices (1,2,3...)
  include.colnames           = FALSE,                        # suppressed; replaced by add.to.row header
  add.to.row                 = list(pos = list(0),           # inject header before first data row
                                    command = header_row),
  sanitize.text.function     = identity,                     # pass cell values to LaTeX as-is (no escaping)
  file = here("output", "latex", "beta_select.tex")
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Save output ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(beta_comparison, here("output", "fpm_gmt_beta_comparison.csv"))

# THE END
