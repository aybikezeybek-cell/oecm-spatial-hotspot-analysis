# ==============================================================================
# OECM Spatial Gap & Kernel Density Hotspot Analysis: National Scale (Turkey)
# Author: Aybike ??elik
# Description: National-scale spatial gap analysis filtering GBIF biodiversity
#              occurrences outside terrestrial WDPA polygons and generating a
#              15 km bandwidth KDE hotspot raster (EPSG:32637).
# ==============================================================================

library(rgbif)
library(sf)
library(tidyverse)
library(spatstat)
library(terra)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# 1. Fetch GBIF Occurrences Across Turkey (35.8?????42.2??N, 25.6?????44.8??E) ----------
gbif_res <- occ_data(
  decimalLatitude = "35.8,42.2",
  decimalLongitude = "25.6,44.8",
  hasCoordinate = TRUE,
  limit = 5000
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
  
  gbif_turkey_sf <- gbif_simple %>%
    select(all_of(main_cols)) %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
    st_transform(32637)
  
  st_write(gbif_turkey_sf, "data/processed/turkey_gbif_species.shp", delete_layer = TRUE)
  
  # 2. Spatial Gap Filtering Outside WDPA Polygons -----------------------------
  wdpa_polygons <- read_sf("data/processed/turkey_wdpa_cleaned.shp")
  if (st_crs(wdpa_polygons) != st_crs(gbif_turkey_sf)) {
    gbif_turkey_sf <- st_transform(gbif_turkey_sf, st_crs(wdpa_polygons))
  }
  
  inside_wdpa <- st_join(gbif_turkey_sf, wdpa_polygons, join = st_intersects, left = FALSE)
  oecm_points_turkey <- gbif_turkey_sf %>%
    filter(!gbifID %in% inside_wdpa$gbifID)
  
  st_write(oecm_points_turkey, "data/processed/oecm_potential_turkey_hotspots.shp", delete_layer = TRUE)
  
  # 3. National Kernel Density Estimation (sigma = 15 km, 1 km^2 grid) ---------
  coords <- st_coordinates(oecm_points_turkey)
  bbox   <- st_bbox(oecm_points_turkey)
  window <- owin(xrange = c(bbox["xmin"], bbox["xmax"]),
                 yrange = c(bbox["ymin"], bbox["ymax"]))
  
  points_ppp  <- ppp(x = coords[, 1], y = coords[, 2], window = window)
  kde_density <- density.ppp(points_ppp, sigma = 15000, eps = 1000)
  
  kde_raster <- rast(kde_density)
  crs(kde_raster) <- "EPSG:32637"
  
  writeRaster(kde_raster, "data/processed/oecm_kde_turkey_hotspot.tif", overwrite = TRUE)
}
