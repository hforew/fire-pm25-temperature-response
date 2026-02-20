# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### check countries_gridded_0.1deg_v0.1.nc #############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ncdf4)
library(here)
library(ggplot2)

#setwd(here())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### read in countries_gridded_0.1deg_v0.1.nc #############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f <- here("input", "countries_gridded_0.1deg_v0.1.nc")

stopifnot(file.exists(f))
cat("NetCDF path:\n", f, "\n\n")

# open + read coords
nc  <- nc_open(f)
lat <- ncvar_get(nc, "lat")
lon <- ncvar_get(nc, "lon")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### check extent countries_gridded_0.1deg_v0.1.nc #############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
edge_by_extent <- (min(lat) <= -90 && max(lat) >= 90 &&
                     min(lon) <= -180 && max(lon) >= 180)

dlat <- unique(round(diff(lat), 10))
dlon <- unique(round(diff(lon), 10))

cat("=== Method 2: Extent check ===\n")
cat("lat range:", min(lat), "to", max(lat), "  (n =", length(lat), ")\n")
cat("lon range:", min(lon), "to", max(lon), "  (n =", length(lon), ")\n")
cat("unique dlat:", paste(dlat, collapse = ", "), "\n")
cat("unique dlon:", paste(dlon, collapse = ", "), "\n")
cat("Edges by extent? ", edge_by_extent, "\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### plot countries_gridded_0.1deg_v0.1.nc #############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lat <- ncvar_get(nc, "lat")
lon <- ncvar_get(nc, "lon")
nc_close(nc)

lat_sub <- lat[lat >= 33 & lat <= 42]
lon_sub <- lon[lon >= -125 & lon <= -114]
lat_sub <- lat_sub[seq(1, length(lat_sub), by = 2)]
lon_sub <- lon_sub[seq(1, length(lon_sub), by = 2)]

pts <- expand.grid(lon = lon_sub, lat = lat_sub)

# 0.1°
grid_lon <- seq(floor(min(lon_sub)*10)/10, ceiling(max(lon_sub)*10)/10, by = 0.1)
grid_lat <- seq(floor(min(lat_sub)*10)/10, ceiling(max(lat_sub)*10)/10, by = 0.1)

ggplot() +
  geom_vline(xintercept = grid_lon, linewidth = 0.3, linetype = "dotted") +
  geom_hline(yintercept = grid_lat, linewidth = 0.3, linetype = "dotted") +
  geom_point(data = pts, aes(x = lon, y = lat), size = 0.7) +
  coord_equal() +
  theme_bw() +
  labs(
    title = "Check: Are file coordinates on grid edges?",
    subtitle = "If points sit exactly on dotted grid-line intersections => EDGES",
    x = "lon", y = "lat"
  )