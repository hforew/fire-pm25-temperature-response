# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Regrid Zhao PM2.5, Merge with Existing pop_pm Dataset, Export ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Zhao source : 0.1 deg, centers on x.0
# Target grid : 0.5 deg, centers on x.25/x.75  (720 x 360 grid)
# Method      : area-weighted average via terra::resample(method = "average")
#
# New variables added to pop_pm_combined:
#   pm_2095_SSP245_Zhao  = SSP245 with fire
#   pm_2095_SSP585_Zhao  = SSP585 with fire
#   fpm_2095_SSP245_Zhao = SSP245 (with fire - without fire)
#   fpm_2095_SSP585_Zhao = SSP585 (with fire - without fire)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(here)                                        # relative paths
library(data.table)                                  # fast CSV I/O + merge
library(terra)                                       # raster ops & resample
library(ggplot2)                                     # plotting
library(viridis)                                     # color scales
library(patchwork)                                   # combine plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Paths ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pm_dir       <- here("input/Zhao_etal_2025/gridded_output")     # Zhao inputs
pop_pm_file  <- here("output", "pop_pm_combined_with_park2024.csv")
pop_zhao_file <- here("output", "pop_regrid_zhao_20252010.csv")
out_file     <- here("output", "pop_pm_combined_final.csv")
#fig_dir      <- here("images", "Zhao_checks")                       # save plots here
#dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Target Grid (0.5 deg) ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r_target <- rast(
  xmin = -180, xmax = 180, ymin = -90, ymax = 90,
  resolution = 0.5, crs = "EPSG:4326"
)                                                    # 720 x 360 = 259200 cells
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Helper: CSV -> regridded data.table on 0.5 deg grid ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
regrid_one <- function(csv_file, out_col) {
  dt <- fread(file.path(pm_dir, csv_file))           # lon, lat, pm25
  setnames(dt, tolower(names(dt)))
  
  # CSV -> 0.1 deg SpatRaster
  r_src <- rast(dt[, .(lon, lat, pm25)], type = "xyz", crs = "EPSG:4326")
  names(r_src) <- "val"
  
  # area-weighted resample to 0.5 deg
  r_out <- resample(r_src, r_target, method = "average", threads = TRUE)
  
  # SpatRaster -> data.table; rename by position
  out <- as.data.table(as.data.frame(r_out, xy = TRUE, na.rm = FALSE))
  setnames(out, 1:3, c("lon", "lat", out_col))
  return(out)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Regrid All 4 Scenarios ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wf_245 <- regrid_one("SSP245_with fire.csv",    "pm_2095_SSP245_Zhao")   # with fire
nf_245 <- regrid_one("SSP245_without fire.csv", "pm_nf_245_tmp")         # no fire
wf_585 <- regrid_one("SSP585_with fire.csv",    "pm_2095_SSP585_Zhao")   # with fire
nf_585 <- regrid_one("SSP585_without fire.csv", "pm_nf_585_tmp")         # no fire
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Build Merged Zhao Dataset (4 target variables) ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pm_2095_Zhao <- wf_245[nf_245, on = .(lon, lat)]
pm_2095_Zhao <- pm_2095_Zhao[wf_585, on = .(lon, lat)]
pm_2095_Zhao <- pm_2095_Zhao[nf_585, on = .(lon, lat)]

pm_2095_Zhao[, fpm_2095_SSP245_Zhao := pm_2095_SSP245_Zhao - pm_nf_245_tmp]
pm_2095_Zhao[, fpm_2095_SSP585_Zhao := pm_2095_SSP585_Zhao - pm_nf_585_tmp]
pm_2095_Zhao[, c("pm_nf_245_tmp", "pm_nf_585_tmp") := NULL]

rm(wf_245, nf_245, wf_585, nf_585)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global Sanity-Check Maps ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Two maps: one total PM (SSP5-8.5 with fire), one fire-attributable PM (SSP5-8.5)
# Use log10(x + 1) transform so extreme hotspots don't wash out the rest of the map

# Common map theme
map_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 8),
    axis.title       = element_blank(),
    legend.position  = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.key.height= unit(0.3, "cm"),
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey30")
  )

# Map 1: Total PM2.5 (with fire) under SSP5-8.5, 2095
p_pm <- ggplot(pm_2095_Zhao, aes(x = lon, y = lat, fill = pm_2095_SSP585_Zhao)) +
  geom_raster() +
  scale_fill_viridis_c(
    option = "inferno", trans = "log10",
    name = expression(PM[2.5]~"("*mu*"g/m"^3*")"),
    na.value = "grey95"
  ) +
  coord_quickmap(expand = FALSE) +
  labs(title    = "Total PM2.5 under SSP5-8.5 (with fire), 2095",
       subtitle = "Zhao et al. 2025, regridded 0.1° → 0.5°") +
  map_theme

# Map 2: Fire-attributable PM2.5 under SSP5-8.5, 2095
# fpm can be slightly negative in some cells due to nonlinear chemistry; floor at 0 for log
p_fpm <- ggplot(pm_2095_Zhao,
                aes(x = lon, y = lat, fill = pmax(fpm_2095_SSP585_Zhao, 0))) +
  geom_raster() +
  scale_fill_viridis_c(
    option = "magma", trans = "log10",
    name = expression("Fire "*PM[2.5]~"("*mu*"g/m"^3*")"),
    na.value = "grey95"
  ) +
  coord_quickmap(expand = FALSE) +
  labs(title    = "Fire-attributable PM2.5 under SSP5-8.5, 2095",
       subtitle = "With fire − Without fire (values < 0 clipped for log scale)") +
  map_theme

# Combine and save
p_combined <- p_pm / p_fpm                           # stack vertically
#ggsave(
#  filename = file.path(fig_dir, "zhao_ssp585_check.png"),
#  plot = p_combined,
#  width = 10, height = 10, dpi = 150
#)
#cat("Wrote figure:", file.path(fig_dir, "zhao_ssp585_check.png"), "\n")
print(p_combined)                                    # also display in RStudio
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Quick Numerical Sanity Check ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n--- Global summary stats ---\n")
cat(sprintf("pm_2095_SSP585_Zhao  : mean = %.3f  |  max = %.2f  |  NA = %d\n",
            mean(pm_2095_Zhao$pm_2095_SSP585_Zhao, na.rm = TRUE),
            max(pm_2095_Zhao$pm_2095_SSP585_Zhao,  na.rm = TRUE),
            sum(is.na(pm_2095_Zhao$pm_2095_SSP585_Zhao))))
cat(sprintf("fpm_2095_SSP585_Zhao : mean = %.3f  |  max = %.2f  |  min = %.3f  |  neg cells = %d\n",
            mean(pm_2095_Zhao$fpm_2095_SSP585_Zhao, na.rm = TRUE),
            max(pm_2095_Zhao$fpm_2095_SSP585_Zhao,  na.rm = TRUE),
            min(pm_2095_Zhao$fpm_2095_SSP585_Zhao,  na.rm = TRUE),
            sum(pm_2095_Zhao$fpm_2095_SSP585_Zhao < 0, na.rm = TRUE)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Load Existing pop_pm Dataset ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_pm <- fread(pop_pm_file)
cat("\nLoaded pop_pm_combined_with_park2024.csv:\n")
cat("  rows:", nrow(pop_pm), " cols:", ncol(pop_pm), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Load Zhao Regridded Population Dataset ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_zhao <- fread(pop_zhao_file)

cat("\nLoaded pop_regrid_zhao_20252010.csv:\n")
cat("  rows:", nrow(pop_zhao), " cols:", ncol(pop_zhao), "\n")

# If the population file uses x/y instead of lon/lat, rename them
if (all(c("x", "y") %in% names(pop_zhao))) {
  setnames(pop_zhao, c("x", "y"), c("lon", "lat"))
}

# Check that merge keys exist in both datasets
stopifnot(all(c("lon", "lat") %in% names(pop_zhao)))
stopifnot(all(c("lon", "lat") %in% names(pop_pm)))

# Rename the Zhao population column to pop_zhao_2010
if ("pop" %in% names(pop_zhao)) {
  setnames(pop_zhao, "pop", "pop_zhao_2010")
}

# Keep only lon, lat, and the Zhao 2010 population column
stopifnot("pop_zhao_2010" %in% names(pop_zhao))
pop_zhao <- pop_zhao[, .(lon, lat, pop_zhao_2010)]

# Remove duplicated coordinate pairs if any
pop_zhao <- unique(pop_zhao, by = c("lon", "lat"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge Zhao PM Data into pop_pm by (lon, lat) ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_pm_final <- pm_2095_Zhao[pop_pm, on = .(lon, lat)]

n_matched_pm <- pop_pm_final[!is.na(pm_2095_SSP245_Zhao), .N]
cat("\nPM merge: ", n_matched_pm, "/", nrow(pop_pm_final),
    " rows matched to Zhao PM grid (",
    round(100 * n_matched_pm / nrow(pop_pm_final), 2), "%)\n", sep = "")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Merge Zhao Regridded Population into Final Dataset ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_pm_final <- pop_zhao[pop_pm_final, on = .(lon, lat)]

n_matched_pop <- pop_pm_final[!is.na(pop_zhao_2010), .N]
cat("\nPopulation merge: added column = pop_zhao_2010\n")
cat("Population matched rows: ", n_matched_pop, "/", nrow(pop_pm_final),
    " (", round(100 * n_matched_pop / nrow(pop_pm_final), 2), "%)\n", sep = "")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Reorder Columns Before Writing ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop_cols <- grep("pop", names(pop_pm_final), value = TRUE)
pop_cols <- setdiff(pop_cols, c("lon", "lat"))

other_cols <- setdiff(names(pop_pm_final), c("lon", "lat", pop_cols))
setcolorder(pop_pm_final, c("lon", "lat", pop_cols, other_cols))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Write Final Output ###
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setorder(pop_pm_final, -lat, lon)
fwrite(pop_pm_final, out_file)

cat("\nWrote:", out_file, "\n")
cat("Final dataset: ", nrow(pop_pm_final), " rows x ", ncol(pop_pm_final), " cols\n", sep = "")

### THE END 