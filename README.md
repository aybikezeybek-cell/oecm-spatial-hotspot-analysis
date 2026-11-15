# Spatial Gap & Kernel Density Hotspot Analysis for OECM Identification in Turkey

**Reproducible GIS & Spatial Point-Pattern Pipeline (`R`, `sf`, `rgbif`, `spatstat`, `terra`)**

### Overview
Under the Kunming-Montreal Global Biodiversity Framework (Target 3) and **CBD Decision 14/8**, **Other Effective Area-based Conservation Measures (OECMs)** play a critical role in conserving biodiversity outside formally designated protected areas. 

This repository implements an automated spatial gap analysis workflow to identify candidate OECM biodiversity hotspots across **Turkey** (national scale) and **East Anatolia** (regional scale) by integrating World Database on Protected Areas (**WDPA**) boundaries with **GBIF** species occurrence records.

---

### Analytical Workflow

1. **WDPA Layer Standardization:**
   - Unzips and merges raw WDPA polygon shapefiles, excludes `Proposed` sites, retains terrestrial and coastal zones (`marine %in% c(0, 1)`), repairs invalid geometries (`st_make_valid()`), and projects layers to **UTM Zone 37N (`EPSG:32637`)**.
2. **GBIF Biodiversity Data Retrieval & Cleaning:**
   - Queries georeferenced occurrence records via `rgbif::occ_data()` for both **East Anatolia** (`37.0°–42.1°N, 38.5°–44.8°E`) and **Turkey** (`35.8°–42.2°N, 25.6°–44.8°E`), flattening nested list fields into ArcGIS/QGIS-compatible attribute tables.
3. **Spatial Gap Filtering (`st_intersects`):**
   - Performs a spatial join between GBIF species occurrences and terrestrial WDPA polygons to isolate biodiversity records occurring **outside** existing formal protected area networks.
4. **Kernel Density Estimation (KDE) Hotspot Mapping:**
   - Converts unprotected occurrence coordinates into `spatstat` planar point patterns (`ppp`) and computes continuous intensity surfaces at $1\text{ km}^2$ pixel resolution (`eps = 1000`):
     - **East Anatolia Regional Scale:** Gaussian smoothing bandwidth $\sigma = 10\text{ km}$ (`10,000 m`).
     - **National Scale (Turkey):** Gaussian smoothing bandwidth $\sigma = 15\text{ km}$ (`15,000 m`).
   - Exports ArcGIS Pro / QGIS-ready GeoTIFF rasters (`oecm_kde_hotspot.tif`, `oecm_kde_turkey_hotspot.tif`) and ESRI Shapefiles.

---

### Repository Structure
```text
├── README.md
├── .gitignore
└── scripts/
    ├── 01_oecm_turkey_national.R
    └── 02_oecm_east_anatolia.R
