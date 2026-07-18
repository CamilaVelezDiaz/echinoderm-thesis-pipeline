# ==============================================================================
# DISTRIBUTION MAPS
# ==============================================================================
# Map 1: all occurrence records, coloured by depth zone
# Map 2: species richness AND sampling density per grid cell
#
# These maps serve as the spatial/geographic component of this thesis in
# place of a formal Geographic Representativeness Index (GRih) - see
# Methods note on this scope decision. They illustrate geographic sampling
# patterns and biases descriptively, without requiring the additional
# analytical machinery (e.g. explaining WHY specific areas are under-
# sampled - accessibility, habitat, effort) that a formal index would
# demand.
#
# Input:  echino_wide.csv
# Output: map1_records_by_depth_zone.png
#         map2_species_richness_grid.png
#         map2_sampling_density_grid.png
#         grid_richness_density_summary.csv
# ==============================================================================

# =============================================================================
# SETUP
# =============================================================================

setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")

# =============================================================================
# Library
# =============================================================================
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(ggspatial)  # north arrow + scale bar

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE)

if (!"rank6" %in% names(echino_wide)) {
  rank_map <- c("Species"="Species","Subspecies"="Species","Genus"="Genus","Subgenus"="Genus",
                "Family"="Family","Order"="Order","Class"="Class","Phylum"="Phylum")
  echino_wide <- echino_wide %>% mutate(rank6 = recode(taxonomic_resolution_level, !!!rank_map))
}

zone_order <- c("Continental Shelf", "Upper Slope", "Lower Slope", "Abyssal")

# ------------------------------------------------------------------------------
# LOCATE COORDINATE COLUMNS (defensive - checks common candidate names
# rather than assuming one, since exact column names weren't confirmed
# for this specific file)
# ------------------------------------------------------------------------------

lat_candidates <- c("best_latitude", "decimalLatitude", "latitude", "lat")
lon_candidates <- c("best_longitude", "decimalLongitude", "longitude", "lon", "lng")

lat_col <- lat_candidates[lat_candidates %in% names(echino_wide)][1]
lon_col <- lon_candidates[lon_candidates %in% names(echino_wide)][1]

if (is.na(lat_col) || is.na(lon_col)) {
  stop("Could not find latitude/longitude columns automatically. Columns in ",
       "echino_wide.csv containing 'lat' or 'lon' (case-insensitive):\n",
       paste(grep("lat|lon", names(echino_wide), ignore.case = TRUE, value = TRUE), collapse = ", "),
       "\nEdit lat_col/lon_col below manually to match.")
}

message(sprintf("Using coordinate columns: %s (latitude), %s (longitude)", lat_col, lon_col))

echino_geo <- echino_wide %>%
  rename(lat = all_of(lat_col), lon = all_of(lon_col)) %>%
  filter(!is.na(lat), !is.na(lon))

message(sprintf("%d of %d records have usable coordinates (%.1f%%)",
                nrow(echino_geo), nrow(echino_wide), 100 * nrow(echino_geo) / nrow(echino_wide)))

# ------------------------------------------------------------------------------
# BASEMAP: Australia coastline via rnaturalearth
# ------------------------------------------------------------------------------
# Using scale = "large" (10m resolution, the finest rnaturalearth offers)
# throughout, rather than "medium" (50m). The coarser "medium" coastline was
# the main cause of the "dots on land" visual issue in earlier drafts of
# Map 1: it doesn't resolve small reef cays and islands (Lizard Island,
# Torres Strait islets, etc.) precisely, so genuinely-offshore/island-edge
# points could appear to sit on the mainland.
#
# Loaded once here and reused for both maps below, so the coastline is
# consistent across every figure. (Note: coordinate/land QC itself is now
# computed upstream in ECHINODERM_POST-PROCESSING.R Section 1d, not here -
# see the coord_qc_exclude_recommended check just below.)

australia <- ne_countries(scale = "large", country = "Australia", returnclass = "sf")

# Bounding box with a small buffer around the actual data extent
bbox <- echino_geo %>%
  summarise(xmin = min(lon, na.rm = TRUE) - 1, xmax = max(lon, na.rm = TRUE) + 1,
            ymin = min(lat, na.rm = TRUE) - 1, ymax = max(lat, na.rm = TRUE) + 1)

message(sprintf("\nMap extent: lon %.1f to %.1f, lat %.1f to %.1f",
                bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax))

# ------------------------------------------------------------------------------
# FIGURE SIZE: match the true aspect ratio of the map extent
# ------------------------------------------------------------------------------
# coord_sf() preserves true geographic aspect ratio (correcting for the fact
# that a degree of longitude covers less real distance than a degree of
# latitude, away from the equator). Previously ggsave used a fixed 9x8
# (aspect 1.125), which didn't match the map's true aspect - ggplot then
# padded the panel with blank space top/bottom to fit the true-aspect map
# into that frame, producing the large empty margins seen at output. Fixed
# here by deriving width/height from the actual bounding box and its
# central latitude, and sized up overall so the map fills more of the frame.

map_aspect <- (bbox$xmax - bbox$xmin) * cos(mean(c(bbox$ymin, bbox$ymax)) * pi / 180) /
  (bbox$ymax - bbox$ymin)
fig_height <- 9.5
fig_width  <- fig_height * map_aspect

message(sprintf("Figure size: %.1f x %.1f in (aspect %.2f)", fig_width, fig_height, map_aspect))

# ------------------------------------------------------------------------------
# REGIONAL CONTEXT: neighbouring countries + reference cities
# ------------------------------------------------------------------------------
# Some retained, in-scope records sit close to Papua New Guinea (Torres
# Strait/northern GBR - 818 records at lat > -10, individually verified in
# Section 1c as genuine QLD/Torres Strait collections, not PNG) and New
# Caledonia (2 MTQ CIDARIS records at lon ~166.4degE, retained as Coral Sea
# biogeographic province material). Without any neighbouring landmass shown,
# these can look like they're floating in open ocean for no reason. Adding
# PNG/New Caledonia/Solomon Islands/Vanuatu as background context (shown,
# not analysed) resolves that - the data selection itself is unchanged.
#
# Uses the whole Oceania continent at the same "large" (10m) resolution as
# the Australia layer, then drops Australia from it (already drawn above,
# separately, so it isn't duplicated).

world_context <- ne_countries(scale = "large", continent = "Oceania", returnclass = "sf") %>%
  filter(name != "Australia")

# Reference cities for spatial orientation, spanning the latitudinal extent
# of the study region (Torres Strait to southeast Queensland/northern NSW).
#
# Coordinates below are DELIBERATELY offset inland from each city's true
# location (by roughly 30-60km) rather than using exact city-centre
# coordinates. All six real cities sit right on the coast, in the densest
# part of the point cloud - plotting them at their literal location made
# the markers/labels unreadable, buried under echinoderm occurrence points.
# At this map's scale (spanning ~2,500km of coastline), a 30-60km inland
# nudge is visually imperceptible as a location error but makes every label
# legible. Each point below was checked against the coastline to confirm it
# lands on the mainland, not still offshore. Not intended as precise
# administrative coordinates - orientation markers only.
#
# label_lon/label_lat/hjust/vjust position each city's TEXT LABEL (as
# opposed to lon/lat, which position the diamond marker itself) - placed a
# little further inland again and right-justified so the text sits to the
# left of its marker, on land, rather than an automatic force-directed
# placement (ggrepel) occasionally pushing a label back out over the coast/
# ocean. 
reference_cities <- tribble(
  ~city,             ~lon,   ~lat,    ~label_lon, ~label_lat, ~hjust, ~vjust,
  "Cairns",           145.45, -16.95, 145.30,     -16.95,     1,      0.5,
  "Townsville",       146.45, -19.35, 146.30,     -19.35,     1,      0.5,
  "Mackay",           148.85, -21.20, 148.70,     -21.20,     1,      0.5,
  "Rockhampton",      150.05, -23.40, 149.90,     -23.40,     1,      0.5,
  "Brisbane",         152.55, -27.50, 152.40,     -27.50,     1,      0.5
)

# ------------------------------------------------------------------------------
# QC FLAG: now computed upstream, in ECHINODERM_POST-PROCESSING.R Section 1d
# ------------------------------------------------------------------------------
# coord_land_qc_flag / coord_qc_exclude_recommended are read directly from
# echino_wide.csv rather than recomputed here, so the maps and any other
# downstream analysis stay consistent with a single source of truth. See
# Section 1d for the full reasoning (near-island precision effects are kept;
# only high-uncertainty and far-inland records are flagged for exclusion).

if (!"coord_qc_exclude_recommended" %in% names(echino_geo)) {
  stop("coord_qc_exclude_recommended not found - run Section 1d of ",
       "ECHINODERM_POST-PROCESSING.R (coordinate/land QC flag) before this script.")
}

n_flagged <- sum(echino_geo$coord_qc_exclude_recommended, na.rm = TRUE)
cat(sprintf("\n%d of %d georeferenced records flagged coord_qc_exclude_recommended = TRUE (%.2f%%)\n",
            n_flagged, nrow(echino_geo), 100 * n_flagged / nrow(echino_geo)))
cat("(computed in ECHINODERM_POST-PROCESSING.R Section 1d - see coord_land_qc_flag\n")
cat("for the full category breakdown, e.g. table(echino_geo$coord_land_qc_flag))\n")


# ==============================================================================
# MAP 1: ALL RECORDS, COLOURED BY DEPTH ZONE
# ==============================================================================

echino_geo_zoned <- echino_geo %>%
  filter(depth_zone %in% zone_order) %>%
  mutate(depth_zone = factor(depth_zone, levels = zone_order)) %>%
  filter(!coord_qc_exclude_recommended)  # drop only the QC-flagged low-confidence points

n_no_zone <- nrow(echino_geo) - nrow(echino_geo %>% filter(depth_zone %in% zone_order))
if (n_no_zone > 0) {
  message(sprintf(
    "\nNote: %d georeferenced records have no depth zone assigned (no depth data) and are excluded from Map 1.",
    n_no_zone
  ))
}

map1 <- ggplot() +
  geom_sf(data = world_context, fill = "grey93", colour = "grey70", linewidth = 0.25) +
  geom_sf(data = australia, fill = "grey90", colour = "grey60", linewidth = 0.3) +
  geom_point(data = echino_geo_zoned, aes(x = lon, y = lat, colour = depth_zone),
             size = 0.6, alpha = 0.5) +
  geom_point(data = reference_cities, aes(x = lon, y = lat),
             shape = 18, size = 2.2, colour = "grey20") +
  geom_text(data = reference_cities,
            aes(x = label_lon, y = label_lat, label = city, hjust = hjust, vjust = vjust),
            size = 3, colour = "grey20", fontface = "italic") +
  annotate("text", x = 144.5, y = -23, label = "AUSTRALIA",
           size = 5, colour = "grey45", fontface = "bold") +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  scale_colour_manual(
    values = c("Continental Shelf" = "#2a78d6", "Upper Slope" = "#eda100",
               "Lower Slope" = "#e34948", "Abyssal" = "#4a3aa7"),
    limits = zone_order
  ) +
  labs(x = NULL, y = NULL, colour = "Depth zone",
       title = "Echinoderm occurrence records across the study region",
       subtitle = sprintf("n = %d georeferenced, depth-zoned records", nrow(echino_geo_zoned)),
       caption = "Map shown in geographic coordinates (WGS84); scale bar is approximate\nand most accurate at the map's central latitude.") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1.5, "cm"), width = unit(1.5, "cm"),
                         pad_x = unit(0.3, "in"), pad_y = unit(0.85, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", panel.grid = element_line(colour = "grey92"))

ggsave("map1_records_by_depth_zone.png", map1, width = fig_width, height = fig_height, dpi = 300)
cat("Figure saved: map1_records_by_depth_zone.png\n")


# ==============================================================================
# MAP 2: SPECIES RICHNESS AND SAMPLING DENSITY PER GRID CELL
# ==============================================================================
# Grid cell size: 0.5 x 0.5 degree (roughly 55km x 50km at this latitude) -
# a reasonable balance between spatial resolution and avoiding excessive
# empty cells given the data's clustering pattern. Adjust GRID_SIZE below
# if a coarser/finer view is preferred.

GRID_SIZE <- 0.5  # degrees

# Same QC exclusion as Map 1, applied here too so sampling density and
# species richness are computed from the same positionally-confident set
# of records (previously Map 2 used the full unfiltered echino_geo, which
# was inconsistent with Map 1 - see chat notes).
echino_geo_for_grid <- echino_geo %>% filter(!coord_qc_exclude_recommended)

cat(sprintf("\nMap 2 grid built from %d records (%d excluded by coord_qc_exclude_recommended)\n",
            nrow(echino_geo_for_grid), nrow(echino_geo) - nrow(echino_geo_for_grid)))

echino_sf <- st_as_sf(echino_geo_for_grid, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

grid <- st_make_grid(
  st_as_sfc(st_bbox(c(xmin = bbox$xmin, xmax = bbox$xmax, ymin = bbox$ymin, ymax = bbox$ymax), crs = 4326)),
  cellsize = GRID_SIZE, square = TRUE
) %>%
  st_sf(grid_id = seq_along(.), geometry = .)

# Spatial join: assign each record to its grid cell
echino_gridded <- st_join(echino_sf, grid, join = st_within)

# --- Sampling density: total records per cell ---
density_by_cell <- echino_gridded %>%
  st_drop_geometry() %>%
  filter(!is.na(grid_id)) %>%
  count(grid_id, name = "n_records")

# --- Species richness: distinct species per cell (species-level records only) ---
richness_by_cell <- echino_gridded %>%
  st_drop_geometry() %>%
  filter(!is.na(grid_id), rank6 == "Species", !is.na(accepted_name_final)) %>%
  group_by(grid_id) %>%
  summarise(n_species = n_distinct(accepted_name_final), .groups = "drop")

grid_summary <- grid %>%
  left_join(density_by_cell, by = "grid_id") %>%
  left_join(richness_by_cell, by = "grid_id") %>%
  mutate(n_records = coalesce(n_records, 0L), n_species = coalesce(n_species, 0L))

write_csv(st_drop_geometry(grid_summary), "grid_richness_density_summary.csv")

message(sprintf("\nGrid: %d cells total, %d with at least 1 record, %d with at least 1 species-level record",
                nrow(grid_summary), sum(grid_summary$n_records > 0), sum(grid_summary$n_species > 0)))

# --- Map 2a: sampling density (log scale, since a few cells are far denser
# than most - GBRSBD/CIDARIS survey concentration areas will otherwise
# wash out the rest of the map) ---
map2_density <- ggplot() +
  geom_sf(data = grid_summary %>% filter(n_records > 0),
          aes(fill = n_records), colour = NA) +
  geom_sf(data = world_context, fill = NA, colour = "grey60", linewidth = 0.25) +
  geom_sf(data = australia, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_point(data = reference_cities, aes(x = lon, y = lat),
             shape = 18, size = 2.2, colour = "grey20") +
  geom_text(data = reference_cities,
            aes(x = label_lon, y = label_lat, label = city, hjust = hjust, vjust = vjust),
            size = 3, colour = "grey20", fontface = "italic") +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  scale_fill_viridis_c(option = "plasma", trans = "log10",
                       name = "Records\n(log scale)", na.value = "grey95") +
  labs(x = NULL, y = NULL,
       title = "Sampling density",
       subtitle = sprintf("%.1f\u00b0 x %.1f\u00b0 grid cells", GRID_SIZE, GRID_SIZE),
       caption = "Map shown in geographic coordinates (WGS84); scale bar is approximate\nand most accurate at the map's central latitude.") +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1.4, "cm"), width = unit(1.4, "cm"),
                         pad_x = unit(0.3, "in"), pad_y = unit(0.8, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_line(colour = "grey92"))

ggsave("map2_sampling_density_grid.png", map2_density, width = fig_width, height = fig_height, dpi = 300)
cat("Figure saved: map2_sampling_density_grid.png\n")

# --- Map 2b: species richness ---
map2_richness <- ggplot() +
  geom_sf(data = grid_summary %>% filter(n_species > 0),
          aes(fill = n_species), colour = NA) +
  geom_sf(data = world_context, fill = NA, colour = "grey60", linewidth = 0.25) +
  geom_sf(data = australia, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_point(data = reference_cities, aes(x = lon, y = lat),
             shape = 18, size = 2.2, colour = "grey20") +
  geom_text(data = reference_cities,
            aes(x = label_lon, y = label_lat, label = city, hjust = hjust, vjust = vjust),
            size = 3, colour = "grey20", fontface = "italic") +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  scale_fill_viridis_c(option = "viridis", name = "Species\nrichness", na.value = "grey95") +
  labs(x = NULL, y = NULL,
       title = "Species richness",
       subtitle = sprintf("%.1f\u00b0 x %.1f\u00b0 grid cells, species-level records only", GRID_SIZE, GRID_SIZE),
       caption = "Map shown in geographic coordinates (WGS84); scale bar is approximate\nand most accurate at the map's central latitude.") +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1.4, "cm"), width = unit(1.4, "cm"),
                         pad_x = unit(0.3, "in"), pad_y = unit(0.8, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_line(colour = "grey92"))

ggsave("map2_species_richness_grid.png", map2_richness, width = fig_width, height = fig_height, dpi = 300)
cat("Figure saved: map2_species_richness_grid.png\n")

# ------------------------------------------------------------------------------
# DIAGNOSTIC: correlation between density and richness - a coarse but
# useful check on whether richness patterns are being driven mainly by
# sampling effort rather than genuine spatial diversity variation (worth
# noting explicitly if strongly correlated, as a limitation of this map)
# ------------------------------------------------------------------------------

check_data <- grid_summary %>% st_drop_geometry() %>% filter(n_records > 0)
cor_density_richness <- cor(check_data$n_records, check_data$n_species, method = "spearman")
cat(sprintf("\nSpearman correlation, sampling density vs species richness per cell: %.2f\n",
            cor_density_richness))
cat("(A strong positive correlation here would suggest richness patterns are\n")
cat("driven substantially by sampling effort rather than purely reflecting\n")
cat("genuine spatial diversity variation - worth noting as a limitation if so.)\n")

cat("\n=== Distribution maps complete ===\n")
cat("Outputs: map1_records_by_depth_zone.png, map2_species_richness_grid.png,\n")
cat("map2_sampling_density_grid.png, grid_richness_density_summary.csv\n")