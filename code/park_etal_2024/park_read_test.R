

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ncdf4)   # read NetCDF files
library(here)    # build file paths relative to project root
library(ggplot2) # plotting
library(maps)    # country and state border data for ggplot borders()
library(scales)  # squish() for out-of-bounds color handling


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Two files from Park et al. 2024 GEOS-Chem CLASSIC runs:
#   - classic/obsclim: all PM2.5 sources including fire
#   - withoutfire: PM2.5 with fire emissions excluded
# Taking the difference gives PM2.5 attributable to fire alone.
#
# Both files:
#   Variable: PM25 (ug/m3)
#   Grid: 720 x 360 at 0.5deg resolution (global)
#   Vertical: lowest model level only (cdo select,levidx=1 in file history)
#   Time: single annual mean value (cdo yearavg, year 2016 selected)
#   Note: original GEOS-Chem output was at coarser resolution (~2deg x 2.5deg)
#   and was conservatively remapped to 0.5deg (cdo remapcon). Conservative
#   remapping means edge cells between two source cells get weighted average
#   values, while interior cells carry the original source value unchanged.
#   This produces visible blocks of identical values in the 0.5deg grid.

nc_path_classic_all <- here(
  "input/Park_etal_2024/GEOSChem_output/classic/obsclim",
  "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4"
)
nc_path_classic_nofire <- here(
  "input/Park_etal_2024/GEOSChem_output/withoutfire",
  "05x05_CEDS_2015_on_off_pm25_Surface_Re_yearavg.nc4"
)

# Open the all-sources file, extract coordinates and PM25, then close.
# ncvar_get() with singleton dimensions (lev=1, time=1) returns a 2D
# matrix [lon x lat] directly — no need to index the extra dimensions.
nc_classic      <- nc_open(nc_path_classic_all)
print(nc_classic)
lon             <- ncvar_get(nc_classic, "lon") # length 720
lat             <- ncvar_get(nc_classic, "lat") # length 360
pm_classic_all  <- ncvar_get(nc_classic, "PM25") # matrix [720 x 360]
nc_close(nc_classic)

# Verify grid spacing is uniform — confirmed 0.5deg in both directions.
# If spacing were irregular, geom_tile would misalign cells.
cat("lon spacing - unique diffs:\n")
print(unique(round(diff(lon), 6)))
cat("lat spacing - unique diffs:\n")
print(unique(round(diff(lat), 6)))

# Open the no-fire file and extract PM25.
nc_classic_nofire  <- nc_open(nc_path_classic_nofire)
print(nc_classic_nofire)
pm_classic_nofire  <- ncvar_get(nc_classic_nofire, "PM25")
nc_close(nc_classic_nofire)

# Build long-format data frames for ggplot.
# expand.grid() creates every lon-lat combination (lon varies fastest),
# matching the column-major storage order of the [lon x lat] matrix.
# as.vector() unrolls the matrix in the same column-major order,
# so pm values align correctly with their lon-lat coordinates.
df_classic_all    <- expand.grid(lon = lon, lat = lat)
df_classic_all$pm <- as.vector(pm_classic_all)

df_classic_nofire    <- expand.grid(lon = lon, lat = lat)
df_classic_nofire$pm <- as.vector(pm_classic_nofire)

# Fire PM2.5 = all sources minus no-fire run.
# This isolates the contribution of fire emissions to surface PM2.5.
df_classic_fire    <- df_classic_all
df_classic_fire$pm <- as.vector(pm_classic_all - pm_classic_nofire)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Color Scale #################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# White-to-dark-brown palette matching Park et al. 2024 Fig. paper style.
# trans = "log1p" applies a log(1+x) transform to the data before color
# mapping. This spreads color contrast across low values (0-30 ug/m3)
# where most of the geographic variation occurs, rather than compressing
# everything into the pale end of the scale.
pm25_colours <- c(
  "#FFFFFF", "#F5E6C8", "#E8C890", "#D49050",
  "#B06020", "#7A3010", "#3A0800"
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global Plot All PM  ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# geom_raster() is used for the global plot (no explicit cell size needed
# since the full global grid is regular and no boundary clipping occurs).
# coord_fixed() enforces equal axis scaling so 1 degree lon = 1 degree lat
# visually, giving correct geographic proportions.
ggplot(df_classic_all, aes(x = lon, y = lat, fill = pm)) +
  geom_raster() +
  borders("world", colour = "grey20", linewidth = 0.3) +
  scale_fill_gradientn(
    name     = expression(PM[2.5] ~ (mu * g ~ m^{-3})),
    colours  = pm25_colours,
    trans    = "log1p",
    na.value = "grey80"
  ) +
  coord_fixed(
    xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE
  ) +
  labs(
    title = "CLASSIC - PM2.5 (all sources, 2010s avg)",
    x = NULL, y = NULL
  ) +
  theme_minimal()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Global Plot Fire PM  ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(df_classic_fire, aes(x = lon, y = lat, fill = pm)) +
  geom_raster() +
  borders("world", colour = "grey20", linewidth = 0.3) +
  scale_fill_gradientn(
    name     = expression(PM[2.5] ~ (mu * g ~ m^{-3})),
    colours  = pm25_colours,
    trans    = "log1p",
    na.value = "grey80"
  ) +
  coord_fixed(
    xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE
  ) +
  labs(
    title = "CLASSIC - PM2.5 attributable to fire (2010s avg)",
    x = NULL, y = NULL
  ) +
  theme_minimal()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ USA Plot All PM  ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# usa_mask identifies grid cells within the CONUS bounding box.
# Used to compute color scale limits from USA data only, so the legend
# reflects the range of values actually visible in the plot rather than
# the global range (which is much higher, e.g. India, Sub-Saharan Africa).
usa_mask <- df_classic_all$lon >= -130 & df_classic_all$lon <= -60 &
  df_classic_all$lat >= 24 & df_classic_all$lat <= 50

# geom_tile(width=0.5, height=0.5) explicitly draws each cell at 0.5deg
# size. This is preferred over geom_raster() when cropping to a subregion
# because geom_raster() can produce artifacts at crop boundaries.
# borders("state") draws US state lines; borders("world", regions="usa")
# draws the thicker outer US national border on top.
# oob = squish clamps any values outside limits to the nearest limit color
# rather than rendering them as NA/grey.
ggplot(df_classic_all, aes(x = lon, y = lat, fill = pm)) +
  geom_tile(width = 0.5, height = 0.5) +
  borders("state", colour = "grey20", linewidth = 0.3) +
  borders("world", regions = "usa", colour = "grey10", linewidth = 0.5) +
  scale_fill_gradientn(
    name     = expression(PM[2.5] ~ (mu * g ~ m^{-3})),
    colours  = pm25_colours,
    trans    = "log1p",
    limits   = range(df_classic_all$pm[usa_mask], na.rm = TRUE),
    oob      = squish,
    na.value = "grey80"
  ) +
  coord_fixed(
    xlim = c(-130, -60), ylim = c(24, 50), expand = FALSE
  ) +
  labs(
    title = "CLASSIC - PM2.5 (all sources, 2010s avg) - USA",
    x = NULL, y = NULL
  ) +
  theme_minimal()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ USA Plot Fire PM  ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(df_classic_fire, aes(x = lon, y = lat, fill = pm)) +
  geom_tile(width = 0.5, height = 0.5) +
  borders("state", colour = "grey20", linewidth = 0.3) +
  borders("world", regions = "usa", colour = "grey10", linewidth = 0.5) +
  scale_fill_gradientn(
    name     = expression(PM[2.5] ~ (mu * g ~ m^{-3})),
    colours  = pm25_colours,
    trans    = "log1p",
    limits   = range(df_classic_fire$pm[usa_mask], na.rm = TRUE),
    oob      = squish,
    na.value = "grey80"
  ) +
  coord_fixed(
    xlim = c(-130, -60), ylim = c(24, 50), expand = FALSE
  ) +
  labs(
    title = "CLASSIC - PM2.5 attributable to fire (2010s avg) - USA",
    x = NULL, y = NULL
  ) +
  theme_minimal()


### THE END 
