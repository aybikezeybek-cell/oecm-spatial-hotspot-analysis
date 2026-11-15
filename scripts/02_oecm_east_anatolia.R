# ==============================================================================
# OECM Spatial Gap & Kernel Density Hotspot Analysis: East Anatolia
# Author: Aybike ??elik
# Description: Filters WDPA terrestrial protected areas, retrieves GBIF species
#              occurrences, isolates records outside protected areas (OECM 
#              candidates), and generates 1 km resolution KDE hotspot rasters.
# ==============================================================================

library(sf)
library(tidyverse)
library(rgbif)
library(spatstat)
library(terra)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# 1. Extract and Merge WDPA Protected Area Polygons ----------------------------
zip_files <- list.files("data/raw", pattern = "\\.zip$", full.names = TRUE)
if (length(zip_files) > 0) {
  walk(zip_files, ~ unzip(.x, exdir = "data/raw"))
}

all_shp <- list.files(path = "data/raw", pattern = "(?i)\\.shp$", recursive = TRUE, full.names = TRUE)
poly_files <- all_shp[grep("polygon", all_shp, ignore.case = TRUE)]
if (length(poly_files) == 0) poly_files <- all_shp

tur_cleaned <- poly_files %>%
  map(~ read_sf(.x)) %>%
  bind_rows() %>%
  rename_all(tolower)

if ("status" %in% names(tur_cleaned)) {
  tur_cleaned <- tur_cleaned %>% filter(status != "Proposed")
}
if ("marine" %in% names(tur_cleaned)) {
  # 0 = terrestrial, 1 = coastal/partially marine
  tur_cleaned <- tur_cleaned %>% filter(as.character(marine) %in% c("0", "1"))
}

tur_cleaned <- tur_cleaned %>%
  st_make_valid() %>%
  st_transform(32637) # UTM Zone 37N

st_write(tur_cleaned, "data/processed/turkey_wdpa_cleaned.shp", delete_layer = TRUE)

# 2. Fetch & Clean GBIF Occurrence Records (East Anatolia Bounding Box) --------
gbif_res <- occ_data(
  decimalLatitude = "37.0,42.1",
  decimalLongitude = "38.5,44.8",
  hasCoordinate = TRUE,
  limit = 2000
)
gbif_data <- gbif_res$data

if (!is.null(gbif_data) && nrow(gbif_data) > 0) {
  gbif_simple <- gbif_data %>%
    select(where(~ !is.list(.))) %>%
    filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))
  
  if ("species" %in% names(gbif_simple)) {
    gbif_simple <- gbif_simple %>% filter(!is.na(species) & species != "")
  }
  
  main_cols <- intersect(
    c("gbifID", "species", "scientificName", "kingdom", "phylum", "class",
      "order", "family", "year", "decimalLongitude", "decimalLatitude"),
    names(gbif_simple)
  )
  
  gbif_sf <- gbif_simple %>%
    select(all_of(main_cols)) %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
    st_transform(32637)
  
  st_write(gbif_sf, "data/processed/east_anatolia_gbif_species.shp", delete_layer = TRUE)
}

# 3. Spatial Gap Analysis: Isolate Records Outside WDPA Boundaries -------------
wdpa_polygons <- read_sf("data/processed/turkey_wdpa_cleaned.shp")
gbif_points   <- read_sf("data/processed/east_anatolia_gbif_species.shp")

if (st_crs(wdpa_polygons) != st_crs(gbif_points)) {
  gbif_points <- st_transform(gbif_points, st_crs(wdpa_polygons))
}

points_inside_wdpa <- st_join(gbif_points, wdpa_polygons, join = st_intersects, left = FALSE)
oecm_candidate_points <- gbif_points %>%
  filter(!gbifID %in% points_inside_wdpa$gbifID)

st_write(oecm_candidate_points, "data/processed/oecm_potential_species_hotspots.shp", delete_layer = TRUE)

# 4. Kernel Density Estimation (KDE) Hotspot Raster (sigma = 10 km) ------------
coords <- st_coordinates(oecm_candidate_points)
bbox   <- st_bbox(oecm_candidate_points)
window <- owin(xrange = c(bbox["xmin"], bbox["xmax"]),
               yrange = c(bbox["ymin"], bbox["ymax"]))

points_ppp  <- ppp(x = coords[, 1], y = coords[, 2], window = window)
kde_density <- density.ppp(points_ppp, sigma = 10000, eps = 1000) # 1 km^2 resolution

kde_raster <- rast(kde_density)
crs(kde_raster) <- "EPSG:32637"

writeRaster(kde_raster, "data/processed/oecm_kde_hotspot.tif", overwrite = TRUE)
