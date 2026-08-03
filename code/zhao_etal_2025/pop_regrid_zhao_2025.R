# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Zhao etal Population & Regrid 0.1 -> 0.5 deg ###############
# Source: 0.1 deg, centers x.0, edges x.05/x.95, extent [-180.05, 179.95] lon #
#                                                        [-90.05,   89.95] lat#
# Target: 0.5 deg, centers x.25/x.75, edges x.0/x.5, extent [-180,180]x[-90,90]#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Turn the Zhao et al. 2025 0.1° population .npz into a 0.5° long-format CSV,
#       conserving total population through the regrid.
#   1. Load the .npz via reticulate/numpy;
#      sum over the age dimension to a 2-D total-pop matrix (~3.82 B).
#   2. Flip lat S->N to terra's N->S, wrap as the SOURCE 0.1° raster (centers x.0),
#      build the TARGET 0.5° raster (centers x.25/x.75), and resample with
#      method = "sum" (area-weighted) so the global total is preserved (~0% diff).
#   3. Convert to long data.table (lon, lat, pop), run sanity checks — global before/after,
#      US total and 5 city totals using any-touch overlap (extract touches = TRUE) — and
#      write the CSV.
# Input  : input/Zhao_etal_2025/gridded_output_EM/pop_0.1x0.1_2010.npz
# Output : output/pm_joined/pop_regrid_zhao_20252010.csv
# Execution order:
#   files run before: pop_regrid_park_2024.R, and pm_regrid_park_2024.R
#   files run after: pm_regrid_zhao_2025.R / country_map_pop_pm_final.R ------> (Map gridded data to country codes)

rm(list = ls())


library(here)                                          # relative paths
library(reticulate)                                    # Python interop for .npz
library(terra)                                         # raster regrid
library(data.table)                                    # long-format output
library(rnaturalearth)                                 # country polygons
library(sf)                                            # sf <-> terra bridge

here::i_am("input/Zhao_etal_2025/gridded_output_EM/pop_0.1x0.1_2010.npz")

# ---- 1. Load .npz ----
np       <- import("numpy")                            # numpy module
npz_path <- here("input", "Zhao_etal_2025",
                 "gridded_output_EM", "pop_0.1x0.1_2010.npz")

npz      <- np$load(npz_path, allow_pickle = TRUE)     # load archive
pop_arr  <- npz[["var1"]]                              # 1800 x 3600 x 9

cat("pop_arr dim:", dim(pop_arr), "\n")                # 1800 3600 9
cat("pop_arr total:", sum(pop_arr) / 1e9, "B\n")       # ~3.82

# ---- 2. Collapse age dim: 3D -> 2D total-pop matrix ----
pop_2d <- apply(pop_arr, c(1, 2), sum)                 # 1800 (lat S->N) x 3600 (lon W->E)
cat("pop_2d dim:", dim(pop_2d), "\n")                  # 1800 3600
cat("pop_2d total:", sum(pop_2d) / 1e9, "B\n")         # same ~3.82

# ---- 3. Build SOURCE 0.1 deg raster (centers at x.0) ----
# Rows in pop_2d go S->N, but terra expects N->S, so flip the matrix
pop_2d_flip <- pop_2d[1800:1, ]                        # flip lat axis to N->S

r_zhao <- rast(pop_2d_flip,
               extent = ext(-180.05, 179.95,           # lon edges (center x.0)
                            -90.05,   89.95),          # lat edges (center x.0)
               crs    = "EPSG:4326")                   # WGS84
names(r_zhao) <- "pop"

cat("Zhao first 3 lon centers:", head(xFromCol(r_zhao, 1:3), 3), "\n")
# -180.0 -179.9 -179.8
cat("Zhao first 3 lat centers:", head(yFromRow(r_zhao, 1:3), 3), "\n")
#   89.9   89.8   89.7
cat("Zhao res:", res(r_zhao), "\n")                    # 0.1 0.1

# ---- 4. Build TARGET 0.5 deg raster (centers at x.25/x.75) ----
r_target <- rast(extent     = ext(-180, 180, -90, 90), # edges at x.0/x.5
                 resolution = 0.5,
                 crs        = "EPSG:4326")             # 360 rows x 720 cols

cat("Target first 3 lon centers:", head(xFromCol(r_target, 1:3), 3), "\n")
# -179.75 -179.25 -178.75
cat("Target first 3 lat centers:", head(yFromRow(r_target, 1:3), 3), "\n")
#   89.75   89.25   88.75

# ---- 5. Regrid with area-weighted SUM (preserves total pop) ----
# Source and target misaligned by 0.05 deg; method="sum" area-weights source
# cells into target cells so the global population total is conserved.
r_05 <- resample(r_zhao, r_target, method = "sum")     # 0.5 deg pop raster

# ---- 6. Sanity: total pop conservation check ----
tot_src <- global(r_zhao, "sum", na.rm = TRUE)[1, 1]
tot_tgt <- global(r_05,   "sum", na.rm = TRUE)[1, 1]
cat("Total pop before:", tot_src / 1e9, "B\n")
cat("Total pop after: ", tot_tgt / 1e9, "B\n")
cat("Rel diff:", (tot_tgt - tot_src) / tot_src * 100, "%\n")
# should be ~0

# ---- 7. Convert to long data.table (lat N->S outer, lon W->E inner) ----
dt_05 <- as.data.table(as.data.frame(r_05, xy = TRUE, na.rm = FALSE))
setnames(dt_05, c("x", "y", "pop"), c("lon", "lat", "pop"))
setcolorder(dt_05, c("lon", "lat", "pop"))
setorder(dt_05, -lat, lon)                             # screenshot order

cat("dt_05 rows:", nrow(dt_05), "\n")                  # 259,200 = 720 * 360
print(head(dt_05, 25))                                 # lat=89.75, lon=-179.75, ...

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Sanity Check: Global / US / City-level Pop Totals ###############
# 1. Global total (pre- vs post-regrid)                                       #
# 2. US total: ANY cell intersecting US polygon counts (touches=TRUE)         #
# 3. City totals: ANY cell intersecting each city's buffer (touches=TRUE)     #
#    Cities: NYC, LA, Denver, Anchorage, Honolulu                             #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---- Helper: sum pop over cells intersecting a polygon (ANY overlap) ----
sum_pop_touches <- function(rast, poly) {
  vals <- extract(rast, poly, touches = TRUE)[[2]]     # col 1 = ID, col 2 = pop
  sum(vals, na.rm = TRUE)
}

# ---- 1. Global totals ----
glob_src <- global(r_zhao, "sum", na.rm = TRUE)[1, 1]  # 0.1 deg
glob_tgt <- global(r_05,   "sum", na.rm = TRUE)[1, 1]  # 0.5 deg

# ---- 2. US totals (any-touch) ----
us_sf  <- ne_countries(country = "United States of America",
                       scale   = "medium", returnclass = "sf")
us_vec <- vect(us_sf)                                  # SpatVector

us_01  <- sum_pop_touches(r_zhao, us_vec)              # 0.1 deg, any-touch
us_05  <- sum_pop_touches(r_05,   us_vec)              # 0.5 deg, any-touch

# ---- 3. City totals (1km buffer + any-touch) ----
cities_pt <- data.table(
  name = c("NYC",     "LA",       "Denver",   "Anchorage", "Honolulu"),
  lon  = c(-74.0060, -118.2437, -104.9903,  -149.9003,   -157.8583),
  lat  = c( 40.7128,   34.0522,   39.7392,    61.2181,     21.3069)
)
cities_vec <- vect(cities_pt, geom = c("lon", "lat"), crs = "EPSG:4326")
cities_buf <- buffer(cities_vec, width = 1000)         # 1 km buffer

cities_pt[, pop_01 := sapply(seq_len(.N), function(i)
  sum_pop_touches(r_zhao, cities_buf[i]))]
cities_pt[, pop_05 := sapply(seq_len(.N), function(i)
  sum_pop_touches(r_05,   cities_buf[i]))]
cities_pt[, n_cells_01 := sapply(seq_len(.N), function(i)
  nrow(extract(r_zhao, cities_buf[i], touches = TRUE)))]
cities_pt[, n_cells_05 := sapply(seq_len(.N), function(i)
  nrow(extract(r_05,   cities_buf[i], touches = TRUE)))]

# ---- 4. Final summary cat (raw numbers, no units) ----
cat("\n=============== SANITY CHECK SUMMARY ===============\n")
cat(sprintf("[1] Global pop  | 0.1deg = %.0f\n",              glob_src))
cat(sprintf("                | 0.5deg = %.0f\n",              glob_tgt))
cat(sprintf("                | diff   = %.0f (%.6f%%)\n",
            glob_tgt - glob_src, (glob_tgt - glob_src) / glob_src * 100))

cat(sprintf("[2] US pop      | 0.1deg any-touch = %.0f\n",    us_01))
cat(sprintf("                | 0.5deg any-touch = %.0f\n",    us_05))

cat("[3] City pop (1km buffer, any-touch):\n")
cat(sprintf("     %-10s  %10s  %10s  %6s  %6s\n",
            "city", "pop_0.1", "pop_0.5", "n_0.1", "n_0.5"))
for (i in seq_len(nrow(cities_pt))) {
  cat(sprintf("     %-10s  %10.0f  %10.0f  %6d  %6d\n",
              cities_pt$name[i],
              cities_pt$pop_01[i],     cities_pt$pop_05[i],
              cities_pt$n_cells_01[i], cities_pt$n_cells_05[i]))
}
cat("=====================================================\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Export Regridded 0.5 deg Population to CSV ###############
# Columns: lon, lat, pop (long format, matches screenshot layout)             #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---- 1. Build output path via here() ----
out_dir  <- here("output")                             # Ford_BA_FPM25/output
out_file <- here("output", "pm_joined", "pop_regrid_zhao_20252010.csv") # full file path

if (!dir.exists(out_dir)) dir.create(out_dir,          # create dir if missing
                                     recursive = TRUE)

# ---- 2. Write CSV (data.table::fwrite for speed) ----
fwrite(dt_05, out_file)                                # 259,200 rows, 3 cols

# ---- 3. Confirm write ----
cat("\n=============== EXPORT ===============\n")
cat(sprintf("File:  %s\n",         out_file))
cat(sprintf("Rows:  %d\n",         nrow(dt_05)))
cat(sprintf("Cols:  %s\n",         paste(names(dt_05), collapse = ", ")))
cat(sprintf("Size:  %.2f MB\n",    file.info(out_file)$size / 1024^2))
cat("======================================\n")

##THE END