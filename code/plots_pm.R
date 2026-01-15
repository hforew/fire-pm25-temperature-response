###########################
## Plot select countries ##
###########################


# Remove all objects from the environment
rm(list = ls())

############ Packages #####################################################

library(here)
library(tidyverse)
library(ggplot2)
library(maps)
library(tibble)


############ Import #####################################################

# Import annual average PM2.5 data
annual_ave <- read_csv(here("output", "annual_ave_pm25.csv"))
month4to9_ave <- read_csv(here("output", "month4to9_ave_pm25.csv"))


# Preview the data
head(annual_ave)
str(annual_ave)

############ filter data select countries #####################################################

annual_ave_USA <- annual_ave %>%
  filter(country_code_iso3 == "USA")

month4to9_ave_USA <- month4to9_ave %>%
  filter(country_code_iso3 == "USA")

############ Totals - plot select countries on map #####################################################

# 1. Get world map data
world_map <- map_data("world")

# USA - BASELINE 
world_usa_pm <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = pm_2000)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "PM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average PM2.5 - USA (2000)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_usa_pm.png"), world_usa_pm, width = 10, height = 6)

world_usa_fpm <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2000)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "FPM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA (2000)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_usa_fpm.png"), world_usa_fpm, width = 10, height = 6)

world_usa_fpm_mon4to9 <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = month4to9_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2000)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "FPM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Apr. to Sep. Average FPM2.5 - USA (2000)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_usa_fpm_mon4to9.png"), world_usa_fpm_mon4to9, width = 10, height = 6)

# USA - 2050 RCP45 

world_usa_pm_2050_45 <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = pm_2050_45)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "PM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average PM2.5 - USA (2050) RCP4.5",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_usa_pm_2050_45.png"), world_usa_pm_2050_45, width = 10, height = 6)

world_usa_fpm_2050_45 <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2050_45)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "FPM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA (2050) RCP4.5",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_usa_fpm_2050_45.png"), world_usa_fpm_2050_45, width = 10, height = 6)

world_usa_fpm_mon4to9_2050_45 <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = month4to9_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2050_45)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "FPM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Apr. to Sep. Average FPM2.5 - USA (2050) RCP4.5",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_usa_fpm_mon4to9_2050_45.png"), world_usa_fpm_mon4to9_2050_45, width = 10, height = 6)



# India
annual_ave_IND <- annual_ave %>%
  filter(country_code_iso3 == "IND")

world_ind_pm <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_IND,
            aes(x = lon, y = lat, fill = pm_2000)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "PM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average PM2.5 - India (2000)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_ind_pm.png"), world_ind_pm, width = 10, height = 6)

# China
annual_ave_CHN <- annual_ave %>%
  filter(country_code_iso3 == "CHN")

world_chn_pm <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_CHN,
            aes(x = lon, y = lat, fill = pm_2000)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "PM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average PM2.5 - China (2000)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_chn_pm.png"), world_chn_pm, width = 10, height = 6)

# Russia
annual_ave_RUS <- annual_ave %>%
  filter(country_code_iso3 == "RUS")

world_rus_pm <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_RUS,
            aes(x = lon, y = lat, fill = pm_2000)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "PM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average PM2.5 - Russia (2000)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_rus_pm.png"), world_rus_pm, width = 10, height = 6)

# Global - only mapped countries (exclude ocean/unmapped areas)
annual_ave_mapped <- annual_ave %>%
  filter(!is.na(country_code_iso3))

world_pm <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_mapped,
            aes(x = lon, y = lat, fill = pm_2000)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "PM2.5 (μg/m³)",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average PM2.5 - Global (2000)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_pm.png"), world_pm, width = 12, height = 6)


print("All plots saved to images folder")

############ changes - plot select countries on map #####################################################

############### change relative to base year (2000) ############### 

# USA map only (not world) - RCP4.5: 2050 - 2000
usa_fpm_2050_45_base_chg <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2050_45_base_chg)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA RCP4.5: 2050 - 2000",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2050_45_base_chg.png"), usa_fpm_2050_45_base_chg, width = 10, height = 6)

# USA map only (not world) - RCP8.5: 2050 - 2000
usa_fpm_2050_85_base_chg <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2050_85_base_chg)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA RCP8.5: 2050 - 2000",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2050_85_base_chg.png"), usa_fpm_2050_85_base_chg, width = 10, height = 6)

# USA map only (not world) - RCP4.5: 2100 - 2000
usa_fpm_2100_45_base_chg <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2100_45_base_chg)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA RCP4.5: 2100 - 2000",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2100_45_base_chg.png"), usa_fpm_2100_45_base_chg, width = 10, height = 6)

# USA map only (not world) - RCP8.5: 2100 - 2000
usa_fpm_2100_85_base_chg <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2100_85_base_chg)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA RCP8.5: 2100 - 2000",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2100_85_base_chg.png"), usa_fpm_2100_85_base_chg, width = 10, height = 6)


############### within a year, change between RCPs ########################################

# USA - 2050: RCP8.5 - 4.5
world_usa_fpm_2050_rcp_chg <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2050_rcp_chg)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA 2050: RCP8.5 - RCP4.5",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_usa_fpm_2050_rcp_chg.png"), world_usa_fpm_2050_rcp_chg, width = 10, height = 6)

# USA map only not world  - 2050: RCP8.5 - 4.5
usa_fpm_2050_rcp_chg <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2050_rcp_chg)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    breaks = scales::pretty_breaks(n = 6),
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA 2050: RCP8.5 - RCP4.5",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2050_rcp_chg.png"), usa_fpm_2050_rcp_chg, width = 10, height = 6)


# USA - 2100: RCP8.5 - 4.5
world_usa_fpm_2100_rcp_chg <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2100_rcp_chg)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA 2100: RCP8.5 - RCP4.5",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_usa_fpm_2100_rcp_chg.png"), world_usa_fpm_2100_rcp_chg, width = 10, height = 6)

# USA map only not world  - 2100: RCP8.5 - 4.5
usa_fpm_2100_rcp_chg <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2100_rcp_chg)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    breaks = scales::pretty_breaks(n = 6),
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA 2100: RCP8.5 - RCP4.5",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2100_rcp_chg.png"), usa_fpm_2100_rcp_chg, width = 10, height = 6)

# Global - 2050: RCP8.5 - 4.5 - (exclude ocean/unmapped areas)

world_fpm_2050_rcp_chg <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_mapped,
            aes(x = lon, y = lat, fill = fpm_2050_rcp_chg)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - Global 2050: RCP8.5 - RCP4.5",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_fpm_2050_rcp_chg.png"), world_fpm_2050_rcp_chg, width = 10, height = 6)


# Global - 2100: RCP8.5 - 4.5 - (exclude ocean/unmapped areas)

world_fpm_2100_rcp_chg <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = annual_ave_mapped,
            aes(x = lon, y = lat, fill = fpm_2100_rcp_chg)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - Global 2100: RCP8.5 - RCP4.5",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "world_fpm_2100_rcp_chg.png"), world_fpm_2100_rcp_chg, width = 10, height = 6)


############### within a year, change between RCPs WITH human influence ########################################

# USA map only not world  - 2050: RCP8.5 - 4.5 (human influence)
usa_fpm_2050_rcp_chg_hi <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2050_rcp_chg_hi)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    breaks = scales::pretty_breaks(n = 6),
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA 2050: RCP8.5 - RCP4.5 (human influence)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2050_rcp_chg_hi.png"), usa_fpm_2050_rcp_chg_hi, width = 10, height = 6)

# USA map only not world  - 2100: RCP8.5 - 4.5 (human influence)
usa_fpm_2100_rcp_chg_hi <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2100_rcp_chg_hi)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    breaks = scales::pretty_breaks(n = 6),
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA 2100: RCP8.5 - RCP4.5 (human influence)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2100_rcp_chg_hi.png"), usa_fpm_2100_rcp_chg_hi, width = 10, height = 6)

############### change relative to base year (2000) (human influence) ############### 
# for comparison with Ford et al 2018

# USA map only (not world) - RCP4.5: 2050 - 2000 (human influence)
usa_fpm_2050_45_base_chg_hi <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2050_45_base_chg_hi)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA RCP4.5: 2050 - 2000 (human influence)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2050_45_base_chg_hi.png"), usa_fpm_2050_45_base_chg_hi, width = 10, height = 6)

# USA map only (not world) - RCP8.5: 2050 - 2000 (human influence)
usa_fpm_2050_85_base_chg_hi <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2050_85_base_chg_hi)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA RCP8.5: 2050 - 2000 (human influence)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2050_85_base_chg_hi.png"), usa_fpm_2050_85_base_chg_hi, width = 10, height = 6)

# USA map only (not world) - RCP4.5: 2100 - 2000 (human influence)
usa_fpm_2100_45_base_chg_hi <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2100_45_base_chg_hi)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA RCP4.5: 2100 - 2000 (human influence)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2100_45_base_chg_hi.png"), usa_fpm_2100_45_base_chg_hi, width = 10, height = 6)

# USA map only (not world) - RCP8.5: 2100 - 2000 (human influence)
usa_fpm_2100_85_base_chg_hi <- ggplot() +
  geom_tile(data = annual_ave_USA,
            aes(x = lon, y = lat, fill = fpm_2100_85_base_chg_hi)) +
  scale_fill_gradient2(
    low = "blue",        # Negative values
    mid = "white",       # Zero
    high = "darkred",    # Positive values
    midpoint = 0,
    name = "Change in FPM2.5 (μg/m³)",
    na.value = NA
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Annual Average FPM2.5 - USA RCP8.5: 2100 - 2000 (human influence)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images", "usa_fpm_2100_85_base_chg_hi.png"), usa_fpm_2100_85_base_chg_hi, width = 10, height = 6)






