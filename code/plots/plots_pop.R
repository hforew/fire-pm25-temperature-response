# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOT POPULATION DATA ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ggplot2)
library(maps)
library(tibble)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import pop data
pop_pm_country <- read_csv(here("output", "pop_pm_with_countries.csv"))


# Preview the data
head(pop_pm_country)
str(pop_pm_country)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ filter data select countries #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pop_pm_country_USA <- pop_pm_country %>%
  filter(country_code_iso3 == "USA")

# India
pop_pm_country_IND <- pop_pm_country %>%
  filter(country_code_iso3 == "IND")


# China
pop_pm_country_CHN <- pop_pm_country %>%
  filter(country_code_iso3 == "CHN")

# Russia
pop_pm_country_RUS <- pop_pm_country %>%
  filter(country_code_iso3 == "RUS")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Totals - plot select countries on global map #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get world map data
world_map <- map_data("world")

# 2009 population USA 
world_usa_pop_2009 <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = pop_pm_country_USA,
            aes(x = lon, y = lat, fill = pop_tot_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Total",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Total Population -- USA (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "world_usa_pop_2009.png"), world_usa_pop_2009, width = 10, height = 6)

world_usa_pop_den_2009 <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = pop_pm_country_USA,
            aes(x = lon, y = lat, fill = pop_dens_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Density",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Population Density -- USA (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "world_usa_pop_den_2009.png"), world_usa_pop_den_2009, width = 10, height = 6)


# Global - only mapped countries (exclude ocean/unmapped areas)
pop_pm_country_mapped <- pop_pm_country %>%
  filter(!is.na(country_code_iso3))

world_pop_2009 <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", size = 0.1) +
  geom_tile(data = pop_pm_country_mapped,
            aes(x = lon, y = lat, fill = pop_tot_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Total Population",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Total Population -- World (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "world_pop_2009.png"), world_pop_2009, width = 12, height = 6)


print("All plots saved to images folder")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Totals - specific country maps #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# USA - Population Total (2009)
usa_pop_2009 <- ggplot() +
  geom_tile(data = pop_pm_country_USA,
            aes(x = lon, y = lat, fill = pop_tot_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Total",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Total Population - USA (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "usa_pop_2009.png"), usa_pop_2009, width = 10, height = 6)

# USA - Population Density (2009)
usa_pop_den_2009 <- ggplot() +
  geom_tile(data = pop_pm_country_USA,
            aes(x = lon, y = lat, fill = pop_dens_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Density",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Population Density - USA (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "usa_pop_den_2009.png"), usa_pop_den_2009, width = 10, height = 6)

# India - Population Total (2009)
ind_pop_2009 <- ggplot() +
  geom_tile(data = pop_pm_country_IND,
            aes(x = lon, y = lat, fill = pop_tot_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Total",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Total Population - India (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "ind_pop_2009.png"), ind_pop_2009, width = 10, height = 6)

# India - Population Density (2009)
ind_pop_den_2009 <- ggplot() +
  geom_tile(data = pop_pm_country_IND,
            aes(x = lon, y = lat, fill = pop_dens_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Density",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Population Density - India (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "ind_pop_den_2009.png"), ind_pop_den_2009, width = 10, height = 6)

# China - Population Total (2009)
chn_pop_2009 <- ggplot() +
  geom_tile(data = pop_pm_country_CHN,
            aes(x = lon, y = lat, fill = pop_tot_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Total",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Total Population - China (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "chn_pop_2009.png"), chn_pop_2009, width = 10, height = 6)

# China - Population Density (2009)
chn_pop_den_2009 <- ggplot() +
  geom_tile(data = pop_pm_country_CHN,
            aes(x = lon, y = lat, fill = pop_dens_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Density",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Population Density - China (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "chn_pop_den_2009.png"), chn_pop_den_2009, width = 10, height = 6)

# Russia - Population Total (2009)
rus_pop_2009 <- ggplot() +
  geom_tile(data = pop_pm_country_RUS,
            aes(x = lon, y = lat, fill = pop_tot_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Total",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Total Population - Russia (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "rus_pop_2009.png"), rus_pop_2009, width = 10, height = 6)

# Russia - Population Density (2009)
rus_pop_den_2009 <- ggplot() +
  geom_tile(data = pop_pm_country_RUS,
            aes(x = lon, y = lat, fill = pop_dens_2009)) +
  scale_fill_gradient(low = "yellow", high = "darkred",
                      name = "Population Density",
                      na.value = NA) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Population Density - Russia (2009)",
       x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

ggsave(here("images/pop", "rus_pop_den_2009.png"), rus_pop_den_2009, width = 10, height = 6)

print("All country-specific plots saved to images/pop folder")
