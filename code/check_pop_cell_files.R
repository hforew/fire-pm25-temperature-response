
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########### Confirm if PM data grid cell centers or edges ##########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove all objects from the environment
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages #####################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(ncdf4)
library(here)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##################Detect coordinate variable name (lat / lon)#################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pick_coord_name <- function(nc, what = c("lat","lon")) {
  what <- match.arg(what)
  dims <- names(nc$dim)
  vars <- names(nc$var)
  
  d <- dims[grep(what, dims, ignore.case = TRUE)]
  if (length(d) > 0) return(d[1])
  
  v <- vars[grep(what, vars, ignore.case = TRUE)]
  if (length(v) > 0) return(v[1])
  
  NA_character_
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##################Decide whether coordinates represent grid EDGE or CENTER####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
edge_or_center <- function(first_val, d, axis = c("lat","lon")) {
  axis <- match.arg(axis)
  if (is.na(d) || length(d) == 0) return("unknown")
  
  if (axis == "lat") {
    if (abs(first_val - (-90)) < 1e-6) return("edge")
    if (abs(first_val - (-90 + d/2)) < 1e-6) return("center")
  } else {
    if (abs(first_val - (-180)) < 1e-6) return("edge")
    if (abs(first_val - (-180 + d/2)) < 1e-6) return("center")
  }
  "unknown"
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##################Extract a vector aligned to 10 rows from a variable#########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_var_vec10 <- function(nc, varname, axis_pref = c("lon","lat","none"), n = 10) {
  axis_pref <- match.arg(axis_pref)
  
  v <- nc$var[[varname]]
  dnames <- sapply(v$dim, `[[`, "name")
  dlens  <- sapply(v$dim, `[[`, "len")
  
  # create slicing window
  build_slice <- function(axis_regex) {
    i <- which(grepl(axis_regex, dnames, ignore.case = TRUE))[1]
    if (is.na(i)) return(NULL)
    start <- rep(1, length(dnames))
    count <- rep(1, length(dnames))
    count[i] <- min(n, dlens[i])
    list(start = start, count = count, axis_i = i)
  }
  
  sl <- NULL
  if (axis_pref == "lon") sl <- build_slice("lon")
  if (is.null(sl) && axis_pref %in% c("lon","lat")) sl <- build_slice("lat")
  
  if (!is.null(sl)) {
    x <- ncvar_get(nc, varname, start = sl$start, count = sl$count)
    xv <- as.vector(x)
    
    out <- rep(NA_real_, n)
    out[1:min(n, length(xv))] <- xv[1:min(n, length(xv))]
    return(out)
  }
  
  # fallback: flatten variable
  x <- ncvar_get(nc, varname)
  xv <- as.vector(x)
  out <- rep(NA_real_, n)
  out[1:min(n, length(xv))] <- xv[1:min(n, length(xv))]
  out
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ################## print first 10 rows of ALL variables#######################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
check_nc_print10_allvars <- function(file, n = 10) {
  
  cat("\n====================================================\n")
  cat("FILE:", file, "\n")
  cat("====================================================\n")
  
  nc <- nc_open(file)
  on.exit(nc_close(nc), add = TRUE)
  
  latname <- pick_coord_name(nc, "lat")
  lonname <- pick_coord_name(nc, "lon")
  
  lat <- ncvar_get(nc, latname)
  lon <- ncvar_get(nc, lonname)
  
  dlat <- round(median(diff(lat[1:min(50, length(lat))])), 10)
  dlon <- round(median(diff(lon[1:min(50, length(lon))])), 10)
  
  cat("resolution: dlat =", dlat, " dlon =", dlon, "\n")
  cat("lat type:", edge_or_center(lat[1], dlat, "lat"), " | first lat =", lat[1], "\n")
  cat("lon type:", edge_or_center(lon[1], dlon, "lon"), " | first lon =", lon[1], "\n")
  
  # Decide slicing axis
  axis_pref <- "lon"
  if (length(lon) < 2 && length(lat) >= 2) axis_pref <- "lat"
  if (length(lon) < 2 && length(lat) < 2) axis_pref <- "none"
  
  # Build coordinate columns
  if (axis_pref == "lon") {
    df10 <- data.frame(
      lon = lon[1:min(n, length(lon))],
      lat = rep(lat[1], min(n, length(lon)))
    )
    if (nrow(df10) < n) df10[(nrow(df10)+1):n, ] <- NA
    
  } else if (axis_pref == "lat") {
    df10 <- data.frame(
      lon = rep(lon[1], min(n, length(lat))),
      lat = lat[1:min(n, length(lat))]
    )
    if (nrow(df10) < n) df10[(nrow(df10)+1):n, ] <- NA
    
  } else {
    df10 <- data.frame(
      lon = rep(NA_real_, n),
      lat = rep(NA_real_, n)
    )
  }
  
  # Identify data variables (exclude coordinates and bounds)
  all_vars <- names(nc$var)
  drop_pat <- "(^lat$|^latitude$|^lon$|^longitude$|time|bnds|bounds|nv|nbnd)"
  data_vars <- all_vars[!grepl(drop_pat, all_vars, ignore.case = TRUE)]
  
  cat("\nVariables (data vars):\n")
  cat(paste(data_vars, collapse = ", "), "\n")
  
  # Append each variable as a column
  for (vn in data_vars) {
    vec <- get_var_vec10(nc, vn, axis_pref = axis_pref, n = n)
    df10[[vn]] <- vec
  }
  
  cat("\nFirst 10 rows (lon/lat + ALL vars):\n")
  print(df10)
  
  invisible(df10)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ################## file path #################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

files <- c(
  here("input/population/clmforc.Li_2017_HYDEv3.2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2016_c170828.nc"),
  here("input/population/clmforc.Li_2018_SSP1_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc"),
  here("input/landmask_area/sfcarea4popgrid.nc"),
  here("input/population/gridcell_area_0.5deg.nc")
)

lapply(files, check_nc_print10_allvars)
