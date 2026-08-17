# =============================================================================
# ECHINODERM POST-PROCESSING PIPELINE
# Camila Velez | JCU Marine Biology MSc | Phase 1
# Run this AFTER echinoderm_pipeline.R has produced echino_wide
#
# Note:
#   0. Case-duplicate record_key removal
#   1. Coordinate fixes
#     1a: Coordinate conflict resolution
#     1b: COORDINATE RECOVERY (idigbio geoPoint + locality gazetteer)
#     1c: Coordinate outlier fixes + geographic scope exclusions
#     1d: Coordinate / land QC flag
#   2. Country fixes
#   3. Year fixes
#   4. Basis of record review + standardisation
#   5. Depth fixes (best_min_depth / best_max_depth / depth_zone / depth_median)
#     5a-5e: GBIF artifact correction, zones, median/uncertainty, text fixes
#     5f: Depth recovery - event-sibling imputation
#     5g: Depth recovery - reef flat / intertidal text inference
#     5h: POST-AUDIT FIX - min/max swap correction (3,411 records)
#   6. Name normalisation + WoRMS resolution
#     6.0: Non-echinoderm contaminant removal
#     6a-6g: Normalisation, WoRMS passes, manual corrections, validation
#     6h: Full WoRMS validation of every accepted_name
#     6i: POST-AUDIT FIXES - accepted_name_final corrections (27 names total:
#     ochroleucus, intermedius, microplax capitalisation, 3 lowercase
#     subgenus fixes) + Hemiaster extinct-genus flag +
#     taxonomic_resolution_level recomputed from accepted_name_final
#     6j: echino_class resolution - direct from source class__* columns
#     (case-insensitive) + taxon->class lookup recovery for 144 records
#     with a finer-than-Phylum ID but no source class field
#     6k: Higher classification from WoRMS (kingdom, class, family, order)
#   7. Final save
#   8. Darwin Core Archive Export (occurrence core + eMoF extension)
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================

setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")

# =============================================================================
# Library
# =============================================================================

library(tidyverse)
library(worrms)
library(sf)
library(rnaturalearth)

# =============================================================================
# SOURCE PRIORITY (single source of truth)
# =============================================================================
# Used wherever a per-source conflict needs to be resolved by "trust the
# highest-priority source" rather than by majority vote or a different rule
# (e.g. echino_class direct resolution in Section 6j, best_basisOfRecord and
# best_eventDate in Sections 3/4). Defined ONCE here so every consumer stays
# in sync - previously this was redefined locally in Section 6j only, which
# meant AM_direct (the highest-trust source per the name-conflict criteria
# in Methods) was silently missing from echino_class resolution.
#
# CIDARIS/AM direct museum CMS exports are ranked above the broad aggregators
# (ALA/GBIF/OZCAM/iDigBio/OBIS), matching the source-trustworthiness
# criterion already used for name-conflict resolution (Section 6e/Appendix C).
# =============================================================================

SOURCE_PRIORITY <- c(
  "ALA_GBRSBD", "ALA_other", "ALA_echinodermata", "CSIRO_GBRSBD", "CSIRO_QM",
  "MTQ_CIDARIS", "CIDARIS_QMT", "AM_direct",
  "GBIF_QM", "GBIF_echinodermata", "OZCAM", "iDigBio", "OBIS"
)

# =============================================================================
# SECTION 0: CASE-DUPLICATE RECORD_KEY REMOVAL
#   When multiple records share the same record_key (case-insensitively), retain
#   the version integrated from the most sources (n_sources), breaking ties
#   alphabetically by record_key. This preserves the richest available data for
#   each specimen rather than arbitrarily keeping the first or last occurrence.
# =============================================================================

cat("── Section 0: Case-duplicate removal ──\n")

cat("Case-insensitive duplicates before fix:",
    sum(duplicated(str_to_lower(echino_wide$record_key))), "\n")

echino_wide <- echino_wide %>%
  mutate(record_key_lower = str_to_lower(record_key)) %>%
  group_by(record_key_lower) %>%
  arrange(desc(n_sources), record_key) %>%
  slice(1) %>%
  ungroup() %>%
  select(-record_key_lower)

cat("Records after case deduplication:", nrow(echino_wide), "\n")
cat("Case duplicates remaining:",
    sum(duplicated(str_to_lower(echino_wide$record_key))), "\n\n")


# =============================================================================
# SECTION 1: COORDINATE FIXES
#   OBIS is excluded from the majority-vote pool because OBIS systematically
#   re-projects coordinates from its contributing datasets, sometimes assigning
#   slightly different decimal precision than the original source. Including OBIS
#   in the vote would unfairly penalise correct coordinates that happen to differ
#   at the 3rd decimal place from an OBIS re-projection of the same point.
#   OBIS coordinates are retained as a fallback only when no other source has data.
# =============================================================================

cat("── Section 1a: Coordinate conflict resolution ──\n")

depth_min_cols  <- names(echino_wide)[str_detect(names(echino_wide), "^minimumDepthInMeters__")]
depth_max_cols  <- names(echino_wide)[str_detect(names(echino_wide), "^maximumDepthInMeters__")]
lat_cols        <- names(echino_wide)[str_detect(names(echino_wide), "^decimalLatitude__")]
lon_cols        <- names(echino_wide)[str_detect(names(echino_wide), "^decimalLongitude__")]
sciname_cols    <- names(echino_wide)[str_detect(names(echino_wide), "^scientificName__")]
year_cols_tmp   <- names(echino_wide)[str_detect(names(echino_wide), "^year__")]
family_cols_tmp <- names(echino_wide)[str_detect(names(echino_wide), "^family__")]

lat_mat_num <- echino_wide %>%
  select(all_of(lat_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

conflict_ranges <- apply(as.matrix(lat_mat_num), 1, function(row) {
  vals <- na.omit(row)
  if (length(vals) < 2) return(NA_real_)
  diff(range(vals))
})

cat("Coordinate conflict magnitude:\n")
tibble(lat_range = conflict_ranges) %>%
  filter(!is.na(lat_range)) %>%
  summarise(
    n_with_multiple_coords = n(),
    median_range_deg       = round(median(lat_range), 6),
    pct_under_0001         = round(mean(lat_range < 0.0001) * 100, 1),
    pct_under_001          = round(mean(lat_range < 0.001)  * 100, 1),
    pct_over_1deg          = round(mean(lat_range > 1)      * 100, 1)
  ) %>% print()

gbif_lat_col <- "decimalLatitude__GBIF_echinodermata"
gbif_lon_col <- "decimalLongitude__GBIF_echinodermata"
obis_lat_col <- "decimalLatitude__OBIS"
obis_lon_col <- "decimalLongitude__OBIS"

vote_lat_cols <- setdiff(lat_cols, obis_lat_col)
vote_lon_cols <- setdiff(lon_cols, obis_lon_col)

modal_coord <- function(row_vals) {
  vals <- round(na.omit(as.numeric(row_vals)), 3)
  if (length(vals) == 0) return(NA_real_)
  tbl <- sort(table(vals), decreasing = TRUE)
  as.numeric(names(tbl)[1])
}

resolve_coord <- function(lat_mat_vote, lon_mat_vote,
                          gbif_lat_vec, gbif_lon_vec,
                          lat_mat_all, lon_mat_all,
                          conflict_flag) {
  n <- nrow(lat_mat_vote)
  best_lat <- numeric(n)
  best_lon <- numeric(n)
  
  for (i in seq_len(n)) {
    if (!isTRUE(conflict_flag[i])) {
      best_lat[i] <- mean(as.numeric(lat_mat_all[i, ]), na.rm = TRUE)
      best_lon[i] <- mean(as.numeric(lon_mat_all[i, ]), na.rm = TRUE)
    } else {
      ml <- modal_coord(lat_mat_vote[i, ])
      mn <- modal_coord(lon_mat_vote[i, ])
      
      lat_votes <- round(na.omit(as.numeric(lat_mat_vote[i, ])), 3)
      lon_votes <- round(na.omit(as.numeric(lon_mat_vote[i, ])), 3)
      lat_n_agree <- sum(lat_votes == ml, na.rm = TRUE)
      lon_n_agree <- sum(lon_votes == mn, na.rm = TRUE)
      
      if (!is.na(ml) && lat_n_agree >= 2) {
        best_lat[i] <- ml
      } else if (!is.na(gbif_lat_vec[i])) {
        best_lat[i] <- gbif_lat_vec[i]
      } else {
        best_lat[i] <- mean(as.numeric(lat_mat_all[i, ]), na.rm = TRUE)
      }
      
      if (!is.na(mn) && lon_n_agree >= 2) {
        best_lon[i] <- mn
      } else if (!is.na(gbif_lon_vec[i])) {
        best_lon[i] <- gbif_lon_vec[i]
      } else {
        best_lon[i] <- mean(as.numeric(lon_mat_all[i, ]), na.rm = TRUE)
      }
    }
  }
  
  list(
    lat = if_else(is.nan(best_lat), NA_real_, best_lat),
    lon = if_else(is.nan(best_lon), NA_real_, best_lon)
  )
}

lon_mat_num <- echino_wide %>%
  select(all_of(lon_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

coord_conflict_flag <- conflict_ranges > 1 & !is.na(conflict_ranges)

lat_mat_vote_df <- echino_wide %>%
  select(all_of(vote_lat_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

lon_mat_vote_df <- echino_wide %>%
  select(all_of(vote_lon_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

gbif_lat_vec <- suppressWarnings(as.numeric(echino_wide[[gbif_lat_col]]))
gbif_lon_vec <- suppressWarnings(as.numeric(echino_wide[[gbif_lon_col]]))

resolved <- resolve_coord(
  lat_mat_vote = lat_mat_vote_df,
  lon_mat_vote = lon_mat_vote_df,
  gbif_lat_vec = gbif_lat_vec,
  gbif_lon_vec = gbif_lon_vec,
  lat_mat_all  = lat_mat_num,
  lon_mat_all  = lon_mat_num,
  conflict_flag = coord_conflict_flag
)

echino_wide <- echino_wide %>%
  mutate(
    coord_conflict  = coord_conflict_flag,
    coord_range_deg = round(conflict_ranges, 6),
    best_latitude   = resolved$lat,
    best_longitude  = resolved$lon
  )

cat("Records with meaningful coord conflict (>1°):",
    sum(echino_wide$coord_conflict, na.rm = TRUE), "\n")

echino_wide %>%
  filter(coord_conflict) %>%
  select(record_key, primary_source, sources_all,
         all_of(lat_cols), all_of(lon_cols), all_of(sciname_cols)) %>%
  write_csv("coord_conflicts_to_review.csv")

conflicted_rows <- which(coord_conflict_flag)
if (length(conflicted_rows) > 0) {
  n_majority <- sum(sapply(conflicted_rows, function(i) {
    ml <- modal_coord(lat_mat_vote_df[i, ])
    lat_votes <- round(na.omit(as.numeric(lat_mat_vote_df[i, ])), 3)
    !is.na(ml) && sum(lat_votes == ml, na.rm = TRUE) >= 2
  }))
  cat(sprintf(
    "  Resolved by majority vote (>=2 sources agree): %d (%.0f%%)\n",
    n_majority, n_majority / length(conflicted_rows) * 100
  ))
  cat(sprintf(
    "  Resolved by GBIF fallback: %d (%.0f%%)\n",
    length(conflicted_rows) - n_majority,
    (length(conflicted_rows) - n_majority) / length(conflicted_rows) * 100
  ))
}

# =============================================================================
# SECTION 1b: COORDINATE RECOVERY (idigbio geoPoint + locality gazetteer)
#  NOTE: recovered_coordinates.csv must be present in the working directory.
# This file was produced by a separate coordinate recovery pass and is not
# regenerated by this script. It contains 650 records with coordinates
# recovered from two sources:
#   (1) the idigbio:geoPoint JSON field in the source data, which carries
#       real recorded coordinates not read by the main pipeline (142 records);
#   (2) locality gazetteer geocoding of named-place text in locality__*
#       source columns not captured in verbatimLocality (508 records,
#       spot-checked against Geoscience Australia authoritative sources).
# See coord_recovery_method and coord_uncertainty_m columns for provenance
# and precision flags on each recovered record.
#
#   - Merges in 650 records' worth of recovered coordinates for records that
#     otherwise have no best_latitude/best_longitude:
#   - 142 records: idigbio:geoPoint field, a real recorded coordinate that
#     Section 1a's coordinate-resolution logic never reads (it only checks
#     decimalLatitude__iDigBio/decimalLongitude__iDigBio, missing this
#     separate JSON-structured field). VERIFIED: spot-checked 5 records
#     against the raw idigbio:geoPoint values in echino_wide - all matched
#     exactly. This is a genuine pipeline gap, not an estimate.
#   - 508 records: locality_gazetteer_geocoding, coordinates assigned from
#     named-place text in locality__* source columns that never made it
#     into the merged verbatimLocality field (e.g. "Lizard Island, Day Reef,
#     station 1", "Wheeler Reef, Townsville"). VERIFIED: spot-checked Lizard
#     Island (-14.677, 145.4675) and Heron Island (-23.4423, 151.9148)
#     against authoritative sources - both accurate to within ~150m-1km,
#     well inside their assigned 2-3km uncertainty radii. These are
#     ESTIMATES, not recorded coordinates - flagged via recovery_method and
#     carry coordinateUncertaintyInMeters_assigned so they can be filtered
#     or weighted separately from directly-recorded coordinates in any
#     spatial completeness or precision-sensitive analysis.
# =============================================================================

cat("── Section 1b: Coordinate recovery (geoPoint + gazetteer) ──\n\n")

recovered_coords <- read_csv("recovered_coordinates.csv", show_col_types = FALSE)

cat("Recovered coordinate records loaded:", nrow(recovered_coords), "\n")
cat("  - idigbio_geopoint_recovery (verified exact):",
    sum(recovered_coords$recovery_method == "idigbio_geopoint_recovery"), "\n")
cat("  - locality_gazetteer_geocoding (estimated, spot-checked):",
    sum(recovered_coords$recovery_method == "locality_gazetteer_geocoding"), "\n\n")

n_before_recovery <- sum(!is.na(echino_wide$best_latitude))

echino_wide <- echino_wide %>%
  left_join(
    recovered_coords %>%
      select(record_key, recovered_latitude, recovered_longitude,
             recovery_method, coordinateUncertaintyInMeters_assigned),
    by = "record_key"
  ) %>%
  mutate(
    coord_recovered_flag  = is.na(best_latitude) & !is.na(recovered_latitude),
    coord_recovery_method = if_else(coord_recovered_flag, recovery_method, NA_character_),
    coord_uncertainty_m   = if_else(coord_recovered_flag, coordinateUncertaintyInMeters_assigned, NA_real_),
    best_latitude         = if_else(coord_recovered_flag, recovered_latitude, best_latitude),
    best_longitude        = if_else(coord_recovered_flag, recovered_longitude, best_longitude)
  ) %>%
  select(-recovered_latitude, -recovered_longitude, -recovery_method,
         -coordinateUncertaintyInMeters_assigned)

n_after_recovery <- sum(!is.na(echino_wide$best_latitude))

cat("Records with best_latitude before recovery:", n_before_recovery, "\n")
cat("Records with best_latitude after recovery: ", n_after_recovery, "\n")
cat("Net coordinates gained:", n_after_recovery - n_before_recovery, "\n\n")

cat("coord_recovered_flag breakdown:\n")
echino_wide %>% count(coord_recovered_flag) %>% print()

cat("\ncoord_recovery_method breakdown (among recovered records):\n")
echino_wide %>%
  filter(coord_recovered_flag) %>%
  count(coord_recovery_method) %>%
  print()

# Reconciliation placed here so the reader sees the discrepancy
# (loaded: 142, applied: 48) before reading the explanation
idigbio_loaded  <- sum(recovered_coords$recovery_method == "idigbio_geopoint_recovery")
idigbio_applied <- sum(echino_wide$coord_recovery_method == "idigbio_geopoint_recovery",
                       na.rm = TRUE)

cat("\niDigBio geopoint reconciliation:\n")
cat("  Loaded from recovery file:                   ", idigbio_loaded, "\n")
cat("  Applied (record had no prior coordinates):   ", idigbio_applied, "\n")
cat("  Already had coordinates from another source: ",
    idigbio_loaded - idigbio_applied, "\n")
cat("  (coord_recovered_flag = FALSE for these —",
    "left_join found them but is.na(best_latitude) was already FALSE)\n\n")

cat("Remaining records still missing coordinates:",
    sum(is.na(echino_wide$best_latitude)), "\n\n")

# =============================================================================
# SECTION 1C: COORDINATE OUTLIER FIXES + GEOGRAPHIC SCOPE
# =============================================================================

cat("── Section 1c: Coordinate outlier fixes + geographic scope ──\n")

cat("Coordinate range check:\n")
cat("Lat < -29:  ", sum(echino_wide$best_latitude  < -29,  na.rm = TRUE), "\n")
cat("Lat > -10:  ", sum(echino_wide$best_latitude  > -10,  na.rm = TRUE), "\n")
cat("Lon < 130:  ", sum(echino_wide$best_longitude < 130,  na.rm = TRUE), "\n")
cat("Lon > 160:  ", sum(echino_wide$best_longitude > 160,  na.rm = TRUE), "\n")
cat("Lon 50-51:  ", sum(echino_wide$best_longitude > 50 &
                          echino_wide$best_longitude < 51,   na.rm = TRUE), "\n")

echino_wide %>%
  filter(!is.na(best_latitude) & !is.na(best_longitude)) %>%
  filter(
    best_latitude  < -29 | best_latitude  > -10 |
      best_longitude < 130 | best_longitude > 160
  ) %>%
  select(record_key, primary_source, best_latitude,
         best_longitude, country, verbatimLocality) %>%
  arrange(best_longitude) %>%
  print(n = Inf)

# Fix 1: Longitude 50° → 150° transcription error
echino_wide <- echino_wide %>%
  mutate(
    best_longitude = if_else(
      !is.na(best_longitude) & best_longitude > 50 & best_longitude < 51,
      best_longitude + 100,
      best_longitude
    )
  )
cat("Longitude 50° → 150° corrected\n")

# Fix 2: Remove NSW and WA outliers only
cat("Removing out-of-scope records (NSW + WA)\n")
echino_wide <- echino_wide %>%
  filter(!record_key %in% c(
    "MCZ:IZ:150955",
    "MCZ:IZ:150956",
    "URN:LSID:OZCAM.TAXONOMY.ORG.AU:QM:OTHERINVERTS:G233800"
  ))

cat("Coordinate range check:\n")
cat("Lat < -29:  ", sum(echino_wide$best_latitude  < -29,  na.rm = TRUE), "\n")
cat("Lat > -10:  ", sum(echino_wide$best_latitude  > -10,  na.rm = TRUE), "\n")
cat("Lon < 135:  ", sum(echino_wide$best_longitude < 135,  na.rm = TRUE), "\n")
cat("Lon > 167:  ", sum(echino_wide$best_longitude > 167,  na.rm = TRUE), "\n")
cat("Lon 50-51:  ", sum(echino_wide$best_longitude > 50 &
                          echino_wide$best_longitude < 51,   na.rm = TRUE), "\n")

# Diagnostic check BEFORE scope exclusions — depth columns removed,
# not yet resolved at this point in the script
echino_wide %>%
  filter(!is.na(best_longitude) & best_longitude < 135) %>%
  select(record_key, primary_source, best_latitude, best_longitude,
         verbatimLocality, country) %>%
  print()

echino_wide %>%
  filter(!is.na(best_latitude) & best_latitude > -10) %>%
  select(record_key, primary_source, best_latitude, best_longitude,
         verbatimLocality, country) %>%
  print()

# Fix 3: Geographic scope exclusions
cat("── Geographic scope exclusions ──\n")

scope_exclusions <- c(
  "CAT_G144076",                                                  # Darwin NT lon 131°E
  "09036793-E479-4A30-AD50-E3FE5E9FFFF5",                         # Inland (Russell river, between cairns and townsville) fungus MEL herbarium
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:AM:INVERTEBRATES - MARINE & OTHER:W.24503", # Inland (Yabbra State Forest NSW) Extinct Disaster genus (Echinodermata) 
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:AM:INVERTEBRATES - MARINE & OTHER:W.24504", # Inland (Mt Warning NSW) Extinct Disaster genus (Echinodermata)
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:AM:INVERTEBRATES - MARINE & OTHER:W.24548", # Inland (Yabbra State Forest NSW) Extinct Disaster genus (Echinodermata)
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:AM:INVERTEBRATES - MARINE & OTHER:W.24549", # Inland (Mt Warning NSW) Extinct Disaster genus (Echinodermata)
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:MAGNT:ECHINODERM:Q003029", # York Sound, WA - stateProvince incorrectly says Queensland; no coords to confirm, locality text consistent across 4 sources (ALA, GBIF, iDigBio, OZCAM)
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:MAGNT:ECHINODERM:Q003083", # York Sound, WA - same issue
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:MAGNT:ECHINODERM:Q003247",  # York Sound, WA - same issue
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:QM:OTHERINVERTS:G212779"        # Western Australia lon 118°E
)

echino_wide <- echino_wide %>%
  filter(!record_key %in% scope_exclusions)

cat("Records removed:", length(scope_exclusions), "\n")
cat("Records remaining:", nrow(echino_wide), "\n\n")

# Records outside the expected bounding box that are RETAINED after review:
#
# Lat < -29 (1 record): iDigBio specimen at -29.45°, 154°E. Locality text
#   places it off the coast of Moreton Bay, southeastern Queensland.
#   Moreton Bay falls within the Queensland marine jurisdiction and is
#   biogeographically relevant to this study's scope (southernmost extent
#   of the Great Barrier Reef Province sensu Briggs & Bowen 2012). Retained.
#
# Lat > -10 (818 records): Torres Strait and northern GBR specimens.
#   Queensland's marine jurisdiction extends to ~-9.1° (Bramble Cay,
#   northernmost point). All records in this range were verified against
#   their locality text and confirmed as genuine Queensland/Torres Strait
#   collections, not misassigned records from Papua New Guinea or Indonesia.
#   Retained.
#
# Lon > 160 (2 records after fixes): CAT_G143740 and 9727180D, both from
#   the MTQ CIDARIS collection with verbatimLocality = "New Caledonia" /
#   "off taku lighthouse near islet Amedée, Nouméa" (lon ~166.4°E).
#   Retained as Coral Sea biogeographic province material consistent with
#   the study scope — same decision applied to all other MTQ CIDARIS
#   New Caledonia records in this dataset (see scope_exclusions notes
#   for the York Sound contrast case).

cat("Final coordinate verification:\n")
cat("Lon < 135:   ", sum(echino_wide$best_longitude < 135, na.rm = TRUE), "\n")
cat("Lon > 167:   ", sum(echino_wide$best_longitude > 167, na.rm = TRUE), "\n")
cat("Lat < -29:   ", sum(echino_wide$best_latitude  < -29, na.rm = TRUE), "(Moreton Bay - confirmed QLD)\n")
cat("Lat > -10:   ", sum(echino_wide$best_latitude  > -10, na.rm = TRUE),
    "(Torres Strait / northern GBR — confirmed legitimate)\n")

# =============================================================================
# SECTION 1D: COORDINATE / LAND QC FLAG
# =============================================================================
# Adds a data-quality flag rather than removing records (consistent with
# is_straddler, crosses_zone_boundary, obs_quality_flag, depth_imputed_flag
# elsewhere in this pipeline). Distinguishes two things that Section 1c's
# bounding-box check did not, and could not, distinguish - both fall well
# within the accepted lat -29 to -10 / lon 135-167 study bounding box:
#
#   1. Genuine small-island/reef-station localities (Lizard Island Research
#      Station, Murray Island/Mer, etc.) where the recorded coordinate sits
#      on the island itself rather than the exact wet substrate nearby.
#      Echinoderms legitimately occur here - NOT an error, kept as-is.
#   2. Records where the coordinate is actually unreliable: self-reported
#      coordinateUncertaintyInMeters >= 10 km, or a coordinate that lands
#      >5 km inland with no such precision effect (e.g. a citizen-science
#      record defaulting to a home address, or a locality/coordinate
#      mismatch between two different real places).
#
# IMPORTANT: which records fall "on land" is itself sensitive to basemap
# resolution choice (documented separately - the single-country boundary
# from ne_countries() undercounts small reef/island features relative to
# Natural Earth's full physical land layer). For that reason the flag below
# is NOT based on "on land" alone - it combines distance-to-coastline with
# self-reported coordinateUncertaintyInMeters, both basemap-resolution-
# independent criteria, so the flag itself is stable regardless of which
# coastline polygon a given analysis happens to use downstream.
# =============================================================================

cat("── Section 1d: Coordinate / land QC flag ──\n")

australia_hires <- ne_countries(scale = "large", country = "Australia", returnclass = "sf")
coastline_boundary <- st_boundary(st_union(australia_hires))

echino_coords_sf <- echino_wide %>%
  filter(!is.na(best_latitude), !is.na(best_longitude)) %>%
  st_as_sf(coords = c("best_longitude", "best_latitude"), crs = 4326, remove = FALSE)

on_land <- st_within(echino_coords_sf, australia_hires, sparse = FALSE)[, 1]
cat(sprintf("%d of %d georeferenced records fall within the coastline polygon (%.2f%%)\n",
            sum(on_land), nrow(echino_coords_sf), 100 * mean(on_land)))

on_land_sf <- echino_coords_sf[on_land, ]
on_land_sf$dist_to_coast_m <- as.numeric(st_distance(on_land_sf, coastline_boundary))

qc_lookup <- on_land_sf %>%
  st_drop_geometry() %>%
  mutate(
    coord_land_qc_flag = case_when(
      dist_to_coast_m <= 2000 ~ "on_land_near_coast_or_island",
      !is.na(coordinateUncertaintyInMeters) & coordinateUncertaintyInMeters >= 10000 ~ "on_land_high_uncertainty",
      dist_to_coast_m > 5000 ~ "on_land_far_inland_review",
      TRUE ~ "on_land_moderate_distance_review"
    )
  ) %>%
  select(record_key, coord_land_qc_flag, dist_to_coast_m)

echino_wide <- echino_wide %>%
  select(-any_of(c("coord_land_qc_flag", "dist_to_coast_m",
                   "coord_qc_exclude_recommended"))) %>%
  left_join(qc_lookup, by = "record_key")

echino_wide <- echino_wide %>%
  mutate(
    coord_land_qc_flag = if_else(is.na(coord_land_qc_flag),
                                 "not_on_land", coord_land_qc_flag),
    # Convenience boolean for anything that wants a simple exclude switch -
    # only the two genuinely-suspect categories are TRUE. Near-coast/island
    # and not-on-land records are both usable as-is.
    coord_qc_exclude_recommended = coord_land_qc_flag %in%
      c("on_land_high_uncertainty", "on_land_far_inland_review")
  )

cat("\ncoord_land_qc_flag distribution:\n")
print(echino_wide %>% count(coord_land_qc_flag, sort = TRUE))

cat(sprintf("\n%d records flagged coord_qc_exclude_recommended = TRUE (%.2f%% of all records)\n",
            sum(echino_wide$coord_qc_exclude_recommended),
            100 * mean(echino_wide$coord_qc_exclude_recommended)))
cat("(These are NOT removed here - flagged only. Exclude at the point of use,\n")
cat("e.g. in distribution-map scripts, with filter(!coord_qc_exclude_recommended).\n\n")


# =============================================================================
# SECTION 2: COUNTRY FIXES
# =============================================================================

cat("── Section 2: Country fixes ──\n")

echino_wide %>%
  count(country, sort = TRUE) %>%
  print(n = Inf)

# Non-Australian countries
echino_wide %>%
  filter(!country %in% c("Australia", "australia",
                         "Commonwealth of Australia")) %>%
  filter(!is.na(best_latitude) & !is.na(best_longitude)) %>%
  group_by(country) %>%
  summarise(
    n         = n(),
    lat_min   = round(min(best_latitude), 2),
    lat_max   = round(max(best_latitude), 2),
    lon_min   = round(min(best_longitude), 2),
    lon_max   = round(max(best_longitude), 2),
    .groups   = "drop"
  ) %>%
  print()

# Papua New Guinea country
# PNG investigation: OBIS assigns country = "Papua New Guinea" to records in the
# Torres Strait border zone where the Australian/PNG boundary is ambiguous.
# The following blocks establish whether these records are genuinely PNG territory
# (Milne Bay: confirmed) or plausibly Australian (NA stateProvince: retained as
# PNG on precautionary grounds - no locality data to confirm either way, and
# coordinates fall within the known OBIS border-zone assignment region).

echino_wide %>%
  filter(country == "Papua New Guinea") %>%
  summarise(
    n           = n(),
    lat_min     = round(min(best_latitude, na.rm = TRUE), 3),
    lat_max     = round(max(best_latitude, na.rm = TRUE), 3),
    lon_min     = round(min(best_longitude, na.rm = TRUE), 3),
    lon_max     = round(max(best_longitude, na.rm = TRUE), 3)
  ) %>% print()

echino_wide %>%
  filter(country == "Papua New Guinea") %>%
  count(stateProvince, sort = TRUE) %>%
  print()

echino_wide %>%
  filter(country == "Papua New Guinea") %>%
  count(primary_source, sort = TRUE) %>%
  print()

echino_wide %>%
  filter(country == "Papua New Guinea" & is.na(stateProvince)) %>%
  summarise(
    n       = n(),
    lat_min = round(min(best_latitude, na.rm = TRUE), 3),
    lat_max = round(max(best_latitude, na.rm = TRUE), 3),
    lon_min = round(min(best_longitude, na.rm = TRUE), 3),
    lon_max = round(max(best_longitude, na.rm = TRUE), 3)
  ) %>% print()

# Depth check for PNG 
cat("Raw depth coverage for these PNG NA-stateProvince records:\n")
echino_wide %>%
  filter(country == "Papua New Guinea" & is.na(stateProvince)) %>%
  mutate(
    has_any_min_raw = rowSums(!is.na(
      select(., all_of(depth_min_cols)) %>%
        mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
    )) > 0,
    has_any_max_raw = rowSums(!is.na(
      select(., all_of(depth_max_cols)) %>%
        mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
    )) > 0
  ) %>%
  summarise(
    n = n(),
    with_any_min_raw = sum(has_any_min_raw),
    with_any_max_raw = sum(has_any_max_raw)
  ) %>%
  print()

echino_wide %>%
  filter(country == "Papua New Guinea" & is.na(stateProvince)) %>%
  select(record_key, best_latitude, best_longitude,
         verbatimLocality) %>%
  head(20) %>% print()

echino_wide %>%
  filter(is.na(country)) %>%
  summarise(
    n         = n(),
    lat_min   = round(min(best_latitude, na.rm = TRUE), 2),
    lat_max   = round(max(best_latitude, na.rm = TRUE), 2),
    lon_min   = round(min(best_longitude, na.rm = TRUE), 2),
    lon_max   = round(max(best_longitude, na.rm = TRUE), 2),
    no_coords = sum(is.na(best_latitude))
  ) %>%
  print()

# =============================================================================
# Summary of decisions:
# French Polynesia (7)   → fix to Australia (coords in Coral Sea, mislabelled)
# UNITED STATES (5)      → fix to Australia (collector nationality, not location)
# PNG - Milne Bay (133)  → keep as Papua New Guinea (confirmed PNG territory)
# PNG - NA stateProvince (884) → keep as Papua New Guinea (OBIS-assigned,
#                                border zone, no locality data to confirm either way)
# australia (20,456)     → standardise to Australia
# Commonwealth of Australia (31) → standardise to Australia
# NA (4,698)             → assign Australia (all within bbox, confirmed above)
# =============================================================================

echino_wide <- echino_wide %>%
  mutate(
    country = case_when(
      country %in% c("French Polynesia", "UNITED STATES") ~ "Australia",
      country %in% c("australia", "Commonwealth of Australia") ~ "Australia",
      is.na(country) ~ "Australia",
      TRUE ~ country
    )
  )

cat("Final country distribution:\n")
echino_wide %>% count(country, sort = TRUE) %>% print()
cat("Records remaining:", nrow(echino_wide), "\n\n")


# =============================================================================
# SECTION 3: YEAR FIXES
# =============================================================================

cat("── Section 3: Year fixes (corrected for new record_keys) ──\n")

year_cols <- names(echino_wide)[str_detect(names(echino_wide), "^year__")]

# DIAGNOSTIC

cat("── Diagnostic: scanning for implausible years across all sources ──\n")

year_long_check <- echino_wide %>%
  select(record_key, primary_source, sources_all, all_of(year_cols)) %>%
  mutate(across(all_of(year_cols), as.character)) %>%
  pivot_longer(
    cols = all_of(year_cols),
    names_to = "year_source",
    values_to = "year_value"
  ) %>%
  mutate(year_value_num = suppressWarnings(as.integer(year_value))) %>%
  filter(!is.na(year_value_num))

cat("Years below 1700 (likely data entry error):\n")
year_long_check %>%
  filter(year_value_num < 1700) %>%
  select(record_key, primary_source, sources_all, year_source, year_value_num) %>%
  arrange(year_value_num) %>%
  print(n = Inf)

cat("\nYears between 1700-1860 (check for placeholder/range artifacts):\n")
year_long_check %>%
  filter(year_value_num >= 1700 & year_value_num <= 1860) %>%
  count(year_source, year_value_num, sort = TRUE) %>%
  print(n = Inf)

cat("\nYears after current year (impossible future dates):\n")
year_long_check %>%
  filter(year_value_num > 2026) %>%
  select(record_key, primary_source, sources_all, year_source, year_value_num) %>%
  print(n = Inf)

cat("\nOverall year range across all sources (pre-fix):",
    min(year_long_check$year_value_num), "–",
    max(year_long_check$year_value_num), "\n\n")

# Fix 1: Four iDigBio records carry impossible collection years (154, 954, 983,
# 1472) confirmed as data-entry artefacts — likely catalogue numbers or specimen
# counts accidentally recorded in the year field during digitisation. These are
# targeted by VALUE rather than record_key because the record_keys for these
# specimens changed after the upstream catalogNumber case-insensitivity fix
# (Section 0) collapsed duplicate entries, making the original record_key list
# obsolete. Nulling the iDigBio year leaves the year derivable from other source
# columns if available; if not, the record is correctly treated as undated.

bad_years <- c(154, 954, 983, 1472)

echino_wide <- echino_wide %>%
  mutate(
    `year__iDigBio` = case_when(
      suppressWarnings(as.numeric(`year__iDigBio`)) %in% bad_years ~ NA_real_,
      TRUE ~ suppressWarnings(as.numeric(`year__iDigBio`))
    )
  )

# Verify
all_years <- echino_wide %>%
  select(all_of(year_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.integer(.x)))) %>%
  unlist() %>%
  na.omit()

cat("Year range after fix:", min(all_years), "–", max(all_years), "\n")
cat("Years before 1800:", sum(all_years < 1800), "\n\n")

# Fix 2: 40 MCZ (Museum of Comparative Zoology, Harvard) records accessed via
# iDigBio carry years of 1800, 1830, or 1850 — these are collection date range
# start points assigned during historical batch digitisation, not actual
# collection years. The MCZ digitised 19th-century collections in decade-range
# batches; 1800/1830/1850 are range boundaries, not real dates. These records
# are flagged as is_historical_collection = TRUE before nulling the year so
# their historical provenance is preserved in the dataset even though the
# specific year cannot be recovered.

placeholder_year_keys <- echino_wide %>%
  filter(
    !is.na(suppressWarnings(as.numeric(`year__iDigBio`))) &
      suppressWarnings(as.numeric(`year__iDigBio`)) >= 1800 &
      suppressWarnings(as.numeric(`year__iDigBio`)) <= 1850
  ) %>%
  pull(record_key)

cat("Historical collection placeholders found:", length(placeholder_year_keys), "\n")

echino_wide <- echino_wide %>%
  mutate(
    is_historical_collection = record_key %in% placeholder_year_keys,
    `year__iDigBio` = if_else(
      record_key %in% placeholder_year_keys,
      NA_real_,                  # ← was NA_character_
      `year__iDigBio`
    )
  )

cat("Historical collection flag:", sum(echino_wide$is_historical_collection), "records\n")

# Final year verification
all_years <- echino_wide %>%
  select(all_of(year_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.integer(.x)))) %>%
  unlist() %>%
  na.omit()

cat("Year range after all fixes:", min(all_years), "–", max(all_years), "\n")
cat("Years before 1800:", sum(all_years < 1800), "\n")
cat("Years 1800-1850:  ", sum(all_years >= 1800 & all_years <= 1850), "\n\n")

# ADDITION: RESOLVE best_year (majority-vote consensus column)
# =============================================================================
# Mirrors best_latitude/best_longitude (Section 1a) and best_min_depth/
# best_max_depth (Section 5) - a single, canonical consensus column built
# ONCE here, rather than being re-derived inconsistently by every
# downstream analysis script that happens to need a year field.
#
# Uses majority-vote resolution (matching the modal_coord() logic already
# used for coordinates in Section 1a), not a simple first-non-missing
# coalesce - 4 records were found to have genuine cross-source year
# disagreement (0.01%), each with a clear majority value and exactly one
# dissenting source (iDigBio in 3 of 4 cases, MTQ_CIDARIS in the 4th).
# A first-non-missing approach would have happened to pick the same
# values here by coincidence of column ordering, but majority-vote is the
# principled, order-independent choice.
# =============================================================================

cat("\n── Section 3 addition: resolving best_year ──\n\n")

year_mat <- echino_wide %>%
  select(all_of(year_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

year_conflicts <- apply(as.matrix(year_mat), 1, function(row) {
  vals <- unique(na.omit(row))
  length(vals) > 1
})

cat(sprintf("Records where year sources disagree: %d (%.2f%%)\n",
            sum(year_conflicts), 100 * mean(year_conflicts)))

if (sum(year_conflicts) > 0) {
  conflict_detail <- echino_wide %>%
    select(record_key, all_of(year_cols)) %>%
    filter(year_conflicts)
  write_csv(conflict_detail, "year_source_conflicts.csv")
  cat("Conflicting records written to year_source_conflicts.csv for reference.\n")
}

modal_year <- function(row_vals) {
  vals <- na.omit(row_vals)
  if (length(vals) == 0) return(NA_real_)
  tbl <- sort(table(vals), decreasing = TRUE)
  as.numeric(names(tbl)[1])
}

echino_wide <- echino_wide %>%
  mutate(best_year = apply(as.matrix(year_mat), 1, modal_year))

n_with_year <- sum(!is.na(echino_wide$best_year))
cat(sprintf("Records with best_year resolved: %d of %d (%.1f%%)\n",
            n_with_year, nrow(echino_wide), 100 * n_with_year / nrow(echino_wide)))
cat(sprintf("Year range: %d - %d\n",
            min(echino_wide$best_year, na.rm = TRUE),
            max(echino_wide$best_year, na.rm = TRUE)))

# =============================================================================
# ADDITION: RESOLVE best_eventDate (source-priority consensus column)
# =============================================================================
# eventDate strings vary in format across sources (e.g. "2009-11-11 12:00:00"
# vs bare years vs date ranges), so majority-vote-on-exact-string (as used
# for best_year, above) would rarely agree even where sources mean the same
# date. Instead this takes the highest-priority source's eventDate per
# record, using the single shared SOURCE_PRIORITY defined at the top of the
# script - the same list primary_source assignment was already based on
# (Appendix C), so this is really just "read the field the record's
# priority source actually populated" rather than a separate resolution
# rule of its own.
#
# best_year remains authoritative for year-level analysis; best_eventDate is
# supplementary full-date precision where available, and will likely need
# format cleanup before being published as its own column - check the
# printed sample of distinct values below for mixed formats.
# =============================================================================

cat("\n── Section 3 addition: resolving best_eventDate ──\n\n")

ed_cols <- names(echino_wide)[str_detect(names(echino_wide), "^eventDate__")]

ordered_ed_cols <- intersect(paste0("eventDate__", SOURCE_PRIORITY), ed_cols)
ordered_ed_cols <- c(ordered_ed_cols, setdiff(ed_cols, ordered_ed_cols))

first_non_missing <- function(row_vals) {
  # row_vals arrives already in source-priority order (see ordered_ed_cols)
  vals <- row_vals[!is.na(row_vals) & str_trim(coalesce(row_vals, "")) != ""]
  if (length(vals) == 0) return(NA_character_)
  vals[1]
}

ed_mat <- echino_wide %>% select(all_of(ordered_ed_cols)) %>% as.matrix()

echino_wide <- echino_wide %>%
  mutate(best_eventDate = apply(ed_mat, 1, first_non_missing))

# Time components are frequently placeholder values (e.g. exact noon, or
# 00:01/00:02) rather than real collection times, and aren't needed for this
# study -- keep only the leading YYYY-MM-DD date portion. Values that don't
# match a full ISO date (e.g. a bare year, or a date range) are left as-is
# rather than forced through this pattern.
echino_wide <- echino_wide %>%
  mutate(best_eventDate = if_else(
    str_detect(coalesce(best_eventDate, ""), "^\\d{4}-\\d{2}-\\d{2}"),
    str_extract(best_eventDate, "^\\d{4}-\\d{2}-\\d{2}"),
    best_eventDate
  ))

n_with_ed <- sum(!is.na(echino_wide$best_eventDate))
cat(sprintf("Records with best_eventDate resolved: %d of %d (%.1f%%)\n",
            n_with_ed, nrow(echino_wide), 100 * n_with_ed / nrow(echino_wide)))
cat("Sample of distinct formats present (check for cleanup before publishing):\n")
print(head(unique(echino_wide$best_eventDate), 15))


# =============================================================================
# SECTION 4: BASIS OF RECORD REVIEW + STANDARDISATION
# =============================================================================

cat("── Section 4: Basis of record ──\n")

# Standardise basisOfRecord values to Darwin Core camelCase across all sources
# Applied before any filtering so the distribution printout shows clean values
# and the downstream str_detect filters (which already use (?i) are unaffected)

bor_std <- function(x) {
  case_when(
    str_detect(coalesce(x, ""), "(?i)^preserved.specimen$|^PRESERVED_SPECIMEN$|^preservedspecimen$") ~ "PreservedSpecimen",
    str_detect(coalesce(x, ""), "(?i)^MATERIAL_SAMPLE;PRESERVED_SPECIMEN$")                          ~ "PreservedSpecimen",
    str_detect(coalesce(x, ""), "(?i)^human.observation$|^HUMAN_OBSERVATION$")                       ~ "HumanObservation",
    str_detect(coalesce(x, ""), "(?i)^observation$")                                                 ~ "HumanObservation",
    str_detect(coalesce(x, ""), "(?i)^MATERIAL_SAMPLE$")                                             ~ "MaterialSample",
    str_detect(coalesce(x, ""), "(?i)^living.specimen$|^LIVING_SPECIMEN$")                           ~ "LivingSpecimen",
    str_detect(coalesce(x, ""), "(?i)^fossil.specimen$|^fossilspecimen$")                            ~ "FossilSpecimen",
    str_detect(coalesce(x, ""), "(?i)^machine.observation")                                          ~ "MachineObservation",
    str_detect(coalesce(x, ""), "(?i)^occurrence")                                                   ~ "HumanObservation",
    is.na(x) | x == ""  ~ NA_character_,
    TRUE                ~ x
  )
}

bor_cols    <- names(echino_wide)[str_detect(names(echino_wide), "^basisOfRecord__")]

echino_wide <- echino_wide %>%
  mutate(across(all_of(bor_cols), bor_std))

cat("basisOfRecords standardised to Darwin Core camelCase across all sources.\n")

bor_idigbio <- names(echino_wide)[str_detect(names(echino_wide), "^basisOfRecord__iDigBio")]

cat("Current basisOfRecord distribution:\n")
echino_wide %>%
  mutate(bor = coalesce(!!!lapply(bor_cols, function(col) echino_wide[[col]]))) %>%
  count(bor, sort = TRUE) %>%
  print()

# Check non-PerservedSpecimens records

cat("\nFossil records:\n")
cat(sum(str_detect(str_to_lower(coalesce(
  echino_wide[[bor_idigbio]], "")), "fossil"), na.rm = TRUE), "\n")

cat("\nMachine observation records:\n")
cat(sum(apply(echino_wide %>% select(all_of(bor_cols)), 1,
              function(r) any(str_detect(coalesce(r, ""), "(?i)machine"))),
        na.rm = TRUE), "\n")

n_before <- nrow(echino_wide)

echino_wide <- echino_wide %>%
  filter(!str_detect(str_to_lower(coalesce(
    .data[[bor_idigbio]], "")), "fossil")) %>%
  filter(!apply(select(., all_of(bor_cols)), 1,
                function(r) any(str_detect(coalesce(r, ""), "(?i)machine"))))

cat("Fossil + machine observation records removed:",
    n_before - nrow(echino_wide), "\n")
cat("Records remaining:", nrow(echino_wide), "\n\n")

cat("\nLiving specimen records by source and dataset:\n")
echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""), "(?i)living")))) %>%
  count(primary_source, datasetName, sort = TRUE) %>%
  print(n = Inf)

# Check depth information
echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""), "(?i)living")))) %>%
  select(record_key, primary_source, datasetName, verbatimLocality) %>%
  head(20) %>%
  print()

echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""), "(?i)living")))) %>%
  mutate(is_gbr = str_detect(coalesce(record_key, ""), "^GBR")) %>%
  count(is_gbr) %>%
  print()

n_before <- nrow(echino_wide)

echino_wide <- echino_wide %>%
  filter(!apply(select(., all_of(bor_cols)), 1,
                function(r) any(str_detect(coalesce(r, ""), "(?i)living"))))

cat("Living specimens removed:", n_before - nrow(echino_wide), "\n")
cat("Records remaining:", nrow(echino_wide), "\n\n")

cat("\nMaterial sample sources:\n")
echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""), "(?i)material")))) %>%
  count(primary_source, sort = TRUE) %>%
  print()

# Remove all material samples except 3 confirmed physical specimens
# KEEP (3): ECHOZ097-09, ECHOZ134-09, ECHOZ145-09
#   — Museums Victoria whole specimens used for iBOL barcoding
#   — PRESERVED_SPECIMEN confirmed in GBIF, OBIS, OZCAM
#
# REMOVE (1,651): all other material samples, and  KX duplicate records
# (PreservedSpecimen, GBIF duplicates of QM G22502)

n_before <- nrow(echino_wide)

echino_wide <- echino_wide %>%
  filter(
    !apply(select(., all_of(bor_cols)), 1,
           function(r) any(str_detect(coalesce(r, ""), "(?i)material"))) |
      apply(select(., all_of(bor_cols)), 1,
            function(r) any(str_detect(coalesce(r, ""), "(?i)preserved")))
  ) %>%
  filter(!record_key %in% c("KX844574", "KX844583", "KX844598"))

cat("Material samples removed (basisOfRecord = MaterialSample):",
    n_before - nrow(echino_wide) - 3, "\n")
cat(" - OBIS:", 1508, "| ALA_echinodermata:", 143, "\n")
cat("KX duplicate records removed (PreservedSpecimen, GBIF duplicates of QM G22502): 3\n")
cat("Total removed this step:", n_before - nrow(echino_wide), "\n")

cat("Records remaining:", nrow(echino_wide), "\n\n")

cat("ECHOZ records retained:\n")
echino_wide %>%
  filter(record_key %in% c("ECHOZ097-09", "ECHOZ134-09", "ECHOZ145-09")) %>%
  select(record_key, primary_source) %>%
  print()

cat("BCUH records remaining (should be 0):",
    sum(str_detect(echino_wide$record_key, "^BCUH")), "\n")
cat("CUKES records remaining (should be 0):",
    sum(str_detect(echino_wide$record_key, "^CUKES")), "\n")
cat("KX records remaining (should be 0):",
    sum(echino_wide$record_key %in%
          c("KX844574", "KX844583", "KX844598")), "\n")

# Human observation records

cat("Human observation records by source:\n")
echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""),"(?i)human|(?i)observation")))) %>%
  count(primary_source, sort = TRUE) %>%
  print()

echino_wide <- echino_wide %>%
  mutate(
    is_gbrsbd_human_obs = apply(select(., all_of(bor_cols)), 1,
                                function(r) any(str_detect(coalesce(r, ""),
                                                           "(?i)human|(?i)observation"))) &
      str_detect(coalesce(sources_all, ""), "GBRSBD")
  )

cat("GBRSBD human observations flagged:",
    sum(echino_wide$is_gbrsbd_human_obs), "\n")

cat("── ALA_echinodermata human observations ──\n")
echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""),
                                          "(?i)human|(?i)observation")))) %>%
  filter(primary_source == "ALA_echinodermata") %>%
  count(dataResourceName, sort = TRUE) %>%
  print(n = 20)

cat("\n── OBIS human observations ──\n")
echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""),
                                          "(?i)human|(?i)observation")))) %>%
  filter(primary_source == "OBIS") %>%
  count(dataResourceName, sort = TRUE) %>%
  print(n = 20)

# Raw, unresolved depth coverage check for OBIS human observations —
# used here only as a diagnostic to decide whether these 3,504 records
# carry potentially valuable depth info before Section 5 resolves it properly

cat("Raw (unresolved) depth coverage — OBIS human observations:\n")
echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""),
                                          "(?i)human|(?i)observation")))) %>%
  filter(primary_source == "OBIS") %>%
  mutate(
    has_any_min_raw = rowSums(!is.na(
      select(., all_of(depth_min_cols)) %>%
        mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
    )) > 0,
    has_any_max_raw = rowSums(!is.na(
      select(., all_of(depth_max_cols)) %>%
        mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
    )) > 0
  ) %>%
  summarise(
    n = n(),
    with_any_min_raw = sum(has_any_min_raw),
    with_any_max_raw = sum(has_any_max_raw),
    with_coords      = sum(!is.na(best_latitude))
  ) %>%
  print()

echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""),
                                          "(?i)human|(?i)observation"))) &
           primary_source == "OBIS" &
           is.na(institutionCode) &
           is.na(collectionCode)) %>%
  mutate(
    has_any_min_raw = rowSums(!is.na(
      select(., all_of(depth_min_cols)) %>%
        mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
    )) > 0
  ) %>%
  summarise(
    n = n(),
    with_any_min_raw = sum(has_any_min_raw),
    with_coords      = sum(!is.na(best_latitude))
  ) %>%
  print()

echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""),
                                          "(?i)human|(?i)observation"))) &
           primary_source == "OBIS" &
           is.na(institutionCode) &
           is.na(collectionCode)) %>%
  mutate(name = coalesce(
    `scientificName__ALA_GBRSBD`,
    `scientificName__ALA_other`,
    `scientificName__ALA_echinodermata`,
    `scientificName__GBIF_echinodermata`,
    `scientificName__OBIS`
  )) %>%
  select(record_key, name, verbatimLocality) %>%
  head(20) %>%
  print()

echino_wide %>%
  filter(apply(select(., all_of(bor_cols)), 1,
               function(r) any(str_detect(coalesce(r, ""),
                                          "(?i)human|(?i)observation"))) &
           primary_source == "OBIS" &
           is.na(institutionCode) &
           is.na(collectionCode)) %>%
  count(datasetID, sort = TRUE) %>%
  print(n = 20)

# Remove IMOS National Reference Station zooplankton records: these are net-tow
# zooplankton samples from the Integrated Marine Observing System (IMOS) that
# record echinoderm larvae (primarily holothurian and echinoid plutei) caught in
# plankton nets. They represent larval stages, not adult benthic specimens, and
# their depth fields record net tow depth rather than the benthic habitat depth
# relevant to this study's depth-gradient analysis.

n_before <- nrow(echino_wide)

echino_wide <- echino_wide %>%
  filter(
    !str_detect(coalesce(dataResourceName, ""),
                "IMOS National Reference Station.*Zooplankton") &
      !str_detect(coalesce(collectionCode, ""),
                  "ANMN_NRS_BGC_ZOOPLANKTON")
  )

cat("IMOS zooplankton removed:", n_before - nrow(echino_wide), "\n")
cat("Records remaining:", nrow(echino_wide), "\n\n")

# Three-tier observational quality framework:
#   Flag 1 (specimen-based): physical vouchers in institutional collections
#          (basisOfRecord = PreservedSpecimen) and GBRSBD survey records
#          (Great Barrier Reef Seabed Biodiversity project - expert-identified
#          trawl/dredge specimens with institutional backing).
#   Flag 2 (formal survey): structured scientific surveys without full voucher
#          deposition - AIMS long-term monitoring, CReefs, PROCFish, SPRFMO VME
#          surveys, and similar peer-reviewed dataset sources.
#   Flag 3 (citizen science): opportunistic observations without expert
#          verification - iNaturalist, BowerBird, and equivalent platforms;
#          also OBIS records with no institutionCode or collectionCode, which
#          cannot be attributed to a verifiable collection or survey.
# Assignment is hierarchical: Flag 1 takes precedence over 2, which takes
# precedence over 3. Records not matching any specific criterion default to
# Flag 1 (the majority are PreservedSpecimen).

bor_cols <- names(echino_wide)[str_detect(names(echino_wide), "^basisOfRecord__")]

echino_wide <- echino_wide %>%
  mutate(
    obs_quality_flag = case_when(
      
      apply(select(., all_of(bor_cols)), 1,
            function(r) any(str_detect(coalesce(r, ""),
                                       "(?i)preserved")))                     ~ 1L,
      
      str_detect(coalesce(sources_all, ""), "GBRSBD") ~ 1L,
      
      str_detect(coalesce(dataResourceName, ""),
                 paste(c(
                   "Lizard Island Research Station",
                   "East Torres Strait Seacucumber Survey",
                   "Australian Institute of Marine Science.*Lizard",
                   "Fish and invertebrate communities",
                   "AIMS Long-term Monitoring",
                   "Crown-of-thorns starfish movement",
                   "CSIRO.*Cruises.*Bycatch",
                   "CReefs.*Heron Island",
                   "Vulnerable marine ecosystems.*South Pacific",
                   "WEST CAPE YORK MARINE PARK"
                 ), collapse = "|"))                      ~ 2L,
      
      primary_source == "OBIS" &
        coalesce(collectionCode, "") %in% c(
          "Ophiuroidea",
          "PROCFish/C",
          "LIRS",
          "CMAR_GOC",
          "Comfish_2000_2002",
          "CReefs - Heron Island 2009",
          "CReefs - Heron Island 2008",
          "CReefs - Lizard Island 2009",
          "TS_CRC_Seacucumber_2005",
          "SPRFMO VME"
        )                                        ~ 2L,
      
      primary_source == "OBIS" &
        coalesce(institutionCode, "") ==
        "Institute for Marine and Antarctic Studies (IMAS)" ~ 2L,
      
      str_detect(coalesce(dataResourceName, ""),
                 paste(c(
                   "iNaturalist",
                   "Earth Guardians",
                   "Gaia Guide",
                   "ALA species sightings",
                   "Rocky Shore Citizen Science",
                   "BowerBird",
                   "ClimateWatch",
                   "Diveboard",
                   "Entangled Wildlife"
                 ), collapse = "|"))                      ~ 3L,
      
      primary_source == "OBIS" &
        is.na(institutionCode) &
        is.na(collectionCode)                    ~ 3L,
      
      apply(select(., all_of(bor_cols)), 1,
            function(r) any(str_detect(coalesce(r, ""),
                                       "(?i)human|(?i)observation")))         ~ 3L,
      
      TRUE                                       ~ 1L
    )
  )

cat("Data quality flag distribution:\n")
echino_wide %>%
  count(obs_quality_flag) %>%
  mutate(
    label = case_when(
      obs_quality_flag == 1L ~ "Flag 1 — Specimen-based + GBRSBD",
      obs_quality_flag == 2L ~ "Flag 2 — Formal scientific survey",
      obs_quality_flag == 3L ~ "Flag 3 — Citizen science",
      TRUE                   ~ "Unassigned"
    ),
    pct = round(n / sum(n) * 100, 1)
  ) %>%
  arrange(obs_quality_flag) %>%
  print()

cat("\nGBRSBD flag check:\n")
echino_wide %>%
  filter(str_detect(coalesce(sources_all, ""), "GBRSBD")) %>%
  count(obs_quality_flag) %>%
  print()

cat("\niNaturalist flag check:\n")
echino_wide %>%
  filter(str_detect(coalesce(dataResourceName, ""), "iNaturalist")) %>%
  count(obs_quality_flag) %>%
  print()

cat("\nRecords remaining:", nrow(echino_wide), "\n")

# =============================================================================
# ADDITION: RESOLVE best_basisOfRecord (majority-vote consensus column)
# =============================================================================
# Mirrors best_year (Section 3 addition). Built AFTER bor_std() standardisation
# and AFTER filtering, so it reflects only the values actually retained.
# Uses the single shared SOURCE_PRIORITY defined at the top of the script.
# =============================================================================

cat("\n── Section 4 addition: resolving best_basisOfRecord ──\n\n")

ordered_bor_cols <- intersect(paste0("basisOfRecord__", SOURCE_PRIORITY), bor_cols)
ordered_bor_cols <- c(ordered_bor_cols, setdiff(bor_cols, ordered_bor_cols))

bor_conflicts <- echino_wide %>%
  select(all_of(bor_cols)) %>%
  apply(1, function(row) length(unique(na.omit(row))) > 1)

cat(sprintf("Records where basisOfRecord sources disagree: %d (%.2f%%)\n",
            sum(bor_conflicts), 100 * mean(bor_conflicts)))

if (sum(bor_conflicts) > 0) {
  write_csv(
    echino_wide %>% select(record_key, primary_source, all_of(bor_cols)) %>% filter(bor_conflicts),
    "basisOfRecord_source_conflicts.csv"
  )
  cat("Conflicting records written to basisOfRecord_source_conflicts.csv for reference.\n")
}

modal_priority <- function(row_vals) {
  # row_vals arrives already in source-priority order (see ordered_bor_cols)
  vals <- na.omit(row_vals)
  if (length(vals) == 0) return(NA_character_)
  tbl <- sort(table(vals), decreasing = TRUE)
  if (sum(tbl == tbl[1]) == 1) return(names(tbl)[1])
  tied <- names(tbl)[tbl == tbl[1]]
  row_vals[row_vals %in% tied][1]   # first (highest-priority) tied value
}

bor_mat <- echino_wide %>% select(all_of(ordered_bor_cols)) %>% as.matrix()

echino_wide <- echino_wide %>%
  mutate(best_basisOfRecord = apply(bor_mat, 1, modal_priority))

n_with_bor <- sum(!is.na(echino_wide$best_basisOfRecord))
cat(sprintf("Records with best_basisOfRecord resolved: %d of %d (%.1f%%)\n",
            n_with_bor, nrow(echino_wide), 100 * n_with_bor / nrow(echino_wide)))
print(table(echino_wide$best_basisOfRecord))

# =============================================================================
# SECTION 5: DEPTH FIXES (COMPLETE, FINAL VERSION)
# Camila Velez | JCU MMB | Phase 1
#
# Steps:
#   5a. GBIF midpoint-flattening exclusion 
#   5b. Min/max swap correction: systematic min/max swap correction (7 + 3,411 records)
#       Found via independent cell-by-cell audit of the output CSV: ~96% of
#       cases trace to ALA_other (QM/OZCAM) source records where
#       minimumDepthInMeters/maximumDepthInMeters are sometimes transcribed
#       with the larger value under "minimum". Confirmed genuine, mixed
#       source-data inconsistency (not a uniform mislabel - 66% of
#       ALA_other's dual-depth records are correctly ordered already, so a
#       blanket swap would have been wrong). Differences are mostly small
#       (median 0.5m) and do not cross any depth_zone boundary, but make
#       depth_range_m/depth_uncertainty negative for affected records.
#   5c. Adopt ecologically-grounded depth zones
#   5d. depth_median / depth_uncertainty / is_straddler
#   5e. Negative depth fix + verbatimDepth text-only fixes
#   5f. Final recompute of all depth-derived columns
#   5g. DEPTH RECOVERY STEP 1: event-sibling imputation
#       Records sharing a dated collection event (same non-blank eventDate +
#       same coordinates to 4dp) with a sibling that HAS depth are assigned
#       that sibling's depth, restricted to events where ALL known depths
#       among siblings agree exactly (no conflicting readings). Blank-date
#       records are excluded from grouping entirely, since "same coordinates"
#       alone does not reliably indicate "same sampling event" in this
#       dataset (confirmed: one site held 416 records across 2 sources with
#       genuinely varying depths and no shared date - a false aggregation).
#   5h. DEPTH RECOVERY STEP 2: reef flat / intertidal text inference
#       Records whose verbatimLocality describes the reef flat / intertidal
#       zone are assigned depth = 0m, consistent with the existing
#       verbatimDepth == "Intertidal" logic in 5d. Does NOT match text
#       describing horizontal distance from shore (e.g. "150 metres off
#       Kirra Beach"), which is unrelated to depth.
#
# All records recovered via 5f/5g are flagged in depth_imputed_flag = TRUE,
# distinguishing them from directly measured/reported depth.
#
# =============================================================================

cat("── Section 5: Depth resolution ──\n\n")

# -----------------------------------------------------------------------------
# DISCOVERY SUMMARY (5a):
# All depth_source_conflict records investigated individually.
#
#      GBIF systematic midpoint-flattening:
#      GBIF_echinodermata reports min==max at a single value that equals
#      the exact midpoint of the genuine range reported by other sources
#      (OZCAM, MTQ_CIDARIS, ALA_other, ALA_echinodermata).
#      RESOLUTION: best_min_depth/best_max_depth computed EXCLUDING GBIF
#      whenever any other source has depth data for that record. GBIF is
#      used as fallback only when it is the SOLE source with depth info.
# -----------------------------------------------------------------------------
cat("── Section 5a: GBIF midpoint-flattening exclusion ──\n\n")

depth_min_cols <- names(echino_wide)[str_detect(names(echino_wide), "^minimumDepthInMeters__")]
depth_max_cols <- names(echino_wide)[str_detect(names(echino_wide), "^maximumDepthInMeters__")]

gbif_min_col <- "minimumDepthInMeters__GBIF_echinodermata"
gbif_max_col <- "maximumDepthInMeters__GBIF_echinodermata"
non_gbif_min_cols <- setdiff(depth_min_cols, gbif_min_col)
non_gbif_max_cols <- setdiff(depth_max_cols, gbif_max_col)

non_gbif_min_mat <- echino_wide %>%
  select(all_of(non_gbif_min_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

non_gbif_max_mat <- echino_wide %>%
  select(all_of(non_gbif_max_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

has_non_gbif_min <- rowSums(!is.na(non_gbif_min_mat)) > 0
has_non_gbif_max <- rowSums(!is.na(non_gbif_max_mat)) > 0

best_min_depth_nongbif <- apply(as.matrix(non_gbif_min_mat), 1, function(r) {
  vals <- na.omit(r); if (length(vals) == 0) return(NA_real_); min(vals)
})
best_max_depth_nongbif <- apply(as.matrix(non_gbif_max_mat), 1, function(r) {
  vals <- na.omit(r); if (length(vals) == 0) return(NA_real_); max(vals)
})

gbif_min_vec <- suppressWarnings(as.numeric(echino_wide[[gbif_min_col]]))
gbif_max_vec <- suppressWarnings(as.numeric(echino_wide[[gbif_max_col]]))

echino_wide <- echino_wide %>%
  mutate(
    best_min_depth = if_else(has_non_gbif_min, best_min_depth_nongbif, gbif_min_vec),
    best_max_depth = if_else(has_non_gbif_max, best_max_depth_nongbif, gbif_max_vec),
    depth_quality_note = NA_character_
  )

cat("GBIF-excluded depth recomputed.\n")
cat("Records with best_min_depth:", sum(!is.na(echino_wide$best_min_depth)), "\n")
cat("Records with best_max_depth:", sum(!is.na(echino_wide$best_max_depth)), "\n\n")

# -----------------------------------------------------------------------------
# DISCOVERY SUMMARY (5b):
#   1. Source-level min/max swap (7 records - G138869, G139141, G139190,
#      G139191, G139213, G139518, G139677):
#      OZCAM/original museum record has minimumDepth > maximumDepth.
#      Confirmed as a genuine swap (not GBIF artifact) via GBIF's
#      flattened value matching the midpoint of the swapped values, and
#      geographic validation (coords sit at the outer edge of the
#      Townsville continental shelf, where a 16.9-429.6m trawl range is
#      geographically realistic).
#      RESOLUTION: min/max values swapped back to correct order.
#   2. Found via independent cell-by-cell audit of the output CSV: ~96% of
#      cases trace to ALA_other (QM/OZCAM) source records where
#      minimumDepthInMeters/maximumDepthInMeters are sometimes transcribed
#      with the larger value under "minimum". Confirmed genuine, mixed
#      source-data inconsistency (not a uniform mislabel - 66% of
#      ALA_other's dual-depth records are correctly ordered already, so a
#      blanket swap would have been wrong). Differences are mostly small
#      (median 0.5m) and do not cross any depth_zone boundary, but make
#      depth_range_m/depth_uncertainty negative for affected records.
# -----------------------------------------------------------------------------

# --- Fix 1: 7 swapped min/max records ---
swapped_minmax_keys <- c(
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:QM:OTHERINVERTS:G138869",
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:QM:OTHERINVERTS:G139141",
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:QM:OTHERINVERTS:G139190",
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:QM:OTHERINVERTS:G139191",
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:QM:OTHERINVERTS:G139213",
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:QM:OTHERINVERTS:G139518",
  "URN:LSID:OZCAM.TAXONOMY.ORG.AU:QM:OTHERINVERTS:G139677"
)

echino_wide <- echino_wide %>%
  mutate(
    best_min_depth = if_else(record_key %in% swapped_minmax_keys, 16.9, best_min_depth),
    best_max_depth = if_else(record_key %in% swapped_minmax_keys, 429.6, best_max_depth),
    depth_quality_note = if_else(
      record_key %in% swapped_minmax_keys,
      "Valores invertidos corregidos tras validacion geografica en el talud continental",
      depth_quality_note
    )
  )

cat("Swapped min/max records fixed:", length(swapped_minmax_keys), "\n\n")

# --- Recompute depth_range_m, has_depth ---
echino_wide <- echino_wide %>%
  mutate(
    depth_range_m = if_else(
      !is.na(best_min_depth) & !is.na(best_max_depth),
      best_max_depth - best_min_depth,
      NA_real_
    ),
    has_depth = !is.na(best_min_depth) | !is.na(best_max_depth)
  )

# --- Fix 2: rest of depth min/max swap correction ---

both_depth <- !is.na(echino_wide$best_min_depth) & !is.na(echino_wide$best_max_depth)
n_swapped <- sum(echino_wide$best_min_depth[both_depth] > echino_wide$best_max_depth[both_depth])
cat("Records with best_min_depth > best_max_depth before fix:", n_swapped, "\n")

echino_wide <- echino_wide %>%
  mutate(
    .min_tmp = pmin(best_min_depth, best_max_depth, na.rm = FALSE),
    .max_tmp = pmax(best_min_depth, best_max_depth, na.rm = FALSE),
    best_min_depth = if_else(!is.na(best_min_depth) & !is.na(best_max_depth), .min_tmp, best_min_depth),
    best_max_depth = if_else(!is.na(best_min_depth) & !is.na(best_max_depth), .max_tmp, best_max_depth)
  ) %>%
  select(-.min_tmp, -.max_tmp)

both_depth_after <- !is.na(echino_wide$best_min_depth) & !is.na(echino_wide$best_max_depth)
n_swapped_after <- sum(echino_wide$best_min_depth[both_depth_after] > echino_wide$best_max_depth[both_depth_after])
cat("Records with best_min_depth > best_max_depth after fix:", n_swapped_after, "\n\n")

cat("depth_range_m and has_depth recomputed.\n")
cat("Records with depth:", sum(echino_wide$has_depth), "\n\n")

# -----------------------------------------------------------------------------
# 5c. ECOLOGICALLY-GROUNDED DEPTH ZONES
#   Zone boundaries follow the depth-zonation scheme described for Indo-Pacific
#   echinoderms: the ~200m calcification limit marks the shelf-slope transition
#   (Gage & Tyler, 1991); the 1,000m boundary approximates the
#   upper/lower-slope faunal transition where deposit-feeding holothurians
#   become dominant; and 4,000m separates the slope from abyssal fauna
#   where pressure-specialised taxa dominate. These thresholds are consistent
#   with those used in Australian deep-sea biodiversity assessments (Byrne &
#   O'Hara, 2017).
#
#   REVISED (post Phase-1 review): Mid Slope (1000-3000m) and Deep Slope
#   (3000-4000m) are merged into a single "Lower Slope" zone (1000-4000m).
#   Standard bathyal-zone definitions place the upper/lower bathyal split
#   anywhere from 800m to 3,000m depending on the source, with the whole
#   bathyal commonly spanning 200-4,000m (e.g. Gage & Tyler 1991; Rex &
#   Etter 2010) — so this merge stays within accepted usage rather than
#   introducing an arbitrary split. It also fixes a sample-size problem:
#   the original Deep Slope bin held only 17 records, too few to support
#   any statistical statement about that zone alone.
#
#     Continental Shelf   :    0 - 200 m
#     Upper Slope         :  200 - 1000 m
#     Lower Slope         : 1000 - 4000 m
#     Abyssal             :  > 4000 m
# -----------------------------------------------------------------------------

cat("── Adopting ecologically-grounded depth zones ──\n\n")

echino_wide <- echino_wide %>%
  mutate(
    min_zone = case_when(
      is.na(best_min_depth)  ~ NA_character_,
      best_min_depth <= 200  ~ "Continental Shelf",
      best_min_depth <= 1000 ~ "Upper Slope",
      best_min_depth <= 4000 ~ "Lower Slope",
      TRUE                   ~ "Abyssal"
    ),
    max_zone = case_when(
      is.na(best_max_depth)  ~ NA_character_,
      best_max_depth <= 200  ~ "Continental Shelf",
      best_max_depth <= 1000 ~ "Upper Slope",
      best_max_depth <= 4000 ~ "Lower Slope",
      TRUE                   ~ "Abyssal"
    ),
    depth_zone = case_when(
      is.na(min_zone) & is.na(max_zone) ~ "No depth data",
      is.na(min_zone)                   ~ max_zone,
      is.na(max_zone)                   ~ min_zone,
      min_zone == max_zone              ~ min_zone,
      TRUE ~ paste0("Spans ", min_zone, " to ", max_zone)
    )
  ) %>%
  select(-min_zone, -max_zone)

cat("Depth zone distribution (pre-recovery):\n")
echino_wide %>%
  count(depth_zone, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = Inf)


# -----------------------------------------------------------------------------
# 5d. depth_median (primary continuous variable), depth_uncertainty, is_straddler
# -----------------------------------------------------------------------------

echino_wide <- echino_wide %>%
  mutate(
    depth_median = if_else(
      !is.na(best_min_depth) & !is.na(best_max_depth),
      (best_min_depth + best_max_depth) / 2,
      coalesce(best_min_depth, best_max_depth)
    ),
    depth_uncertainty = depth_range_m,
    is_straddler = !is.na(depth_uncertainty) & depth_uncertainty > 100
  )

cat("\nRecords with depth_median:", sum(!is.na(echino_wide$depth_median)), "\n")
cat("Records flagged as straddlers (uncertainty > 100m):",
    sum(echino_wide$is_straddler, na.rm = TRUE), "\n\n")


# -----------------------------------------------------------------------------
# 5e. NEGATIVE DEPTH FIX + VERBATIMDEPTH TEXT-ONLY FIXES
# -----------------------------------------------------------------------------

cat("── Negative depth check ──\n")
cat("Negative depths present:", sum(echino_wide$best_min_depth < 0, na.rm = TRUE), "\n")

# Fix: URN:CATALOG:CAS:IZ:101407 - ALA shows -0.3m (digitisation error),
# GBIF correctly shows +0.305m (inferred) for the same record
echino_wide <- echino_wide %>%
  mutate(
    best_min_depth = if_else(
      record_key == "URN:CATALOG:CAS:IZ:101407",
      0.305,
      best_min_depth
    ),
    best_max_depth = if_else(
      record_key == "URN:CATALOG:CAS:IZ:101407" & is.na(best_max_depth),
      0.305,
      best_max_depth
    )
  )

cat("Negative depths remaining:", sum(echino_wide$best_min_depth < 0, na.rm = TRUE), "\n\n")

cat("── VerbatimDepth text-only records ──\n")
echino_wide %>%
  filter(!is.na(verbatimDepth) & is.na(best_min_depth) & is.na(best_max_depth)) %>%
  count(verbatimDepth, sort = TRUE) %>%
  print(n = 30)

echino_wide <- echino_wide %>%
  mutate(
    best_min_depth = case_when(
      str_detect(coalesce(verbatimDepth, ""), "(?i)^intertidal$") & is.na(best_min_depth) ~ 0,
      str_detect(coalesce(verbatimDepth, ""), "^0m$") & is.na(best_min_depth) ~ 0,
      str_detect(coalesce(verbatimDepth, ""), "^10-25 m$") & is.na(best_min_depth) ~ 10,
      TRUE ~ best_min_depth
    ),
    best_max_depth = case_when(
      str_detect(coalesce(verbatimDepth, ""), "(?i)^intertidal$") & is.na(best_max_depth) ~ 0,
      str_detect(coalesce(verbatimDepth, ""), "^0m$") & is.na(best_max_depth) ~ 0,
      str_detect(coalesce(verbatimDepth, ""), "^10-25 m$") & is.na(best_max_depth) ~ 25,
      TRUE ~ best_max_depth
    )
  )

cat("Records with depth after verbatimDepth fixes:",
    sum(!is.na(echino_wide$best_min_depth) | !is.na(echino_wide$best_max_depth)), "\n\n")


# -----------------------------------------------------------------------------
# 5f. FINAL RECOMPUTE of all depth-derived columns (pre-recovery state)
# -----------------------------------------------------------------------------

echino_wide <- echino_wide %>%
  mutate(
    depth_range_m = if_else(!is.na(best_min_depth) & !is.na(best_max_depth),
                            best_max_depth - best_min_depth, NA_real_),
    has_depth = !is.na(best_min_depth) | !is.na(best_max_depth)
  ) %>%
  mutate(
    min_zone = case_when(
      is.na(best_min_depth)  ~ NA_character_,
      best_min_depth <= 200  ~ "Continental Shelf",
      best_min_depth <= 1000 ~ "Upper Slope",
      best_min_depth <= 4000 ~ "Lower Slope",
      TRUE                   ~ "Abyssal"
    ),
    max_zone = case_when(
      is.na(best_max_depth)  ~ NA_character_,
      best_max_depth <= 200  ~ "Continental Shelf",
      best_max_depth <= 1000 ~ "Upper Slope",
      best_max_depth <= 4000 ~ "Lower Slope",
      TRUE                   ~ "Abyssal"
    ),
    depth_zone = case_when(
      is.na(min_zone) & is.na(max_zone) ~ "No depth data",
      is.na(min_zone)                   ~ max_zone,
      is.na(max_zone)                   ~ min_zone,
      min_zone == max_zone              ~ min_zone,
      TRUE ~ paste0("Spans ", min_zone, " to ", max_zone)
    )
  ) %>%
  select(-min_zone, -max_zone) %>%
  mutate(
    depth_median = if_else(
      !is.na(best_min_depth) & !is.na(best_max_depth),
      (best_min_depth + best_max_depth) / 2,
      coalesce(best_min_depth, best_max_depth)
    ),
    depth_uncertainty = depth_range_m,
    is_straddler = !is.na(depth_uncertainty) & depth_uncertainty > 100,
    depth_imputed_flag = FALSE  # initialise before recovery steps below
  )

cat("── Pre-recovery depth state ──\n")
cat("Records with depth (has_depth):", sum(echino_wide$has_depth), "\n")
cat("Pct with depth:", round(sum(echino_wide$has_depth) / nrow(echino_wide) * 100, 1), "%\n\n")


# =============================================================================
# 5g. DEPTH RECOVERY STEP 1: EVENT-SIBLING IMPUTATION
# =============================================================================
# Records sharing a dated collection event (same non-blank eventDate + same
# coordinates to 4 decimal places) with a sibling record that HAS depth are
# assigned that sibling's depth - restricted to events where ALL known
# depths among siblings agree exactly. Blank-date records are excluded
# entirely from this grouping: "same coordinates" alone does not reliably
# indicate "same sampling event" in this dataset. This was confirmed by
# inspecting one site (-14.6, 145.0 rounded) holding 416 records from 2
# different sources, all with blank eventDate, spanning genuinely different
# real depths - a false aggregation of unrelated visits to the same site,
# not a single sampling event.
#
# Of 140 dated, multi-record, depth-recoverable events, 110 had fully
# agreeing depths among siblings (safe to impute) and 30 had genuinely
# conflicting depths among siblings (excluded - these represent real
# multi-depth sampling, e.g. dredge transects or sequential reef survey
# passes at different depth bands on the same date, confirmed by direct
# inspection of two example conflicting events).
# =============================================================================

cat("── Section 5f: Depth recovery via event-sibling imputation ──\n\n")

event_date_col <- names(echino_wide)[str_detect(names(echino_wide), "^eventDate__")]

echino_wide_events <- echino_wide %>%
  mutate(
    event_date_clean = coalesce(!!!syms(event_date_col)),
    lat_round = round(best_latitude, 4),
    lon_round = round(best_longitude, 4)
  )

# Identify all dated, multi-record events with at least one depth-bearing
# and one depth-missing sibling
depth_agreement_check <- echino_wide_events %>%
  filter(!is.na(event_date_clean) & event_date_clean != "" &
           !is.na(lat_round) & !is.na(lon_round)) %>%
  group_by(event_date_clean, lat_round, lon_round) %>%
  filter(n() > 1) %>%
  filter(sum(has_depth) > 0 & sum(!has_depth) > 0) %>%
  summarise(
    n_records = n(),
    n_distinct_min_depths = n_distinct(best_min_depth, na.rm = TRUE),
    n_distinct_max_depths = n_distinct(best_max_depth, na.rm = TRUE),
    .groups = "drop"
  )

cat("Recoverable events (dated, multi-record, mixed depth coverage):",
    nrow(depth_agreement_check), "\n")

safe_events <- depth_agreement_check %>%
  filter(n_distinct_min_depths <= 1 & n_distinct_max_depths <= 1)

conflicting_events <- depth_agreement_check %>%
  filter(n_distinct_min_depths > 1 | n_distinct_max_depths > 1)

cat("Safe to impute (siblings agree):", nrow(safe_events), "events\n")
cat("Excluded - conflicting depths among siblings:", nrow(conflicting_events), "events\n\n")

# Build event-level depth lookup from the safe events only
event_depth_lookup <- echino_wide_events %>%
  filter(!is.na(event_date_clean) & event_date_clean != "" &
           !is.na(lat_round) & !is.na(lon_round)) %>%
  semi_join(safe_events, by = c("event_date_clean", "lat_round", "lon_round")) %>%
  filter(has_depth) %>%
  group_by(event_date_clean, lat_round, lon_round) %>%
  summarise(
    imputed_min_depth = first(na.omit(best_min_depth)),
    imputed_max_depth = first(na.omit(best_max_depth)),
    .groups = "drop"
  )

cat("Event-level depth lookup built:", nrow(event_depth_lookup), "events\n")

# Identify the exact record_keys eligible for imputation (avoids relying on
# left_join floating-point coordinate matching, which can silently fail to
# match identical-looking rounded values due to binary float representation)
records_to_impute <- echino_wide_events %>%
  filter(!is.na(event_date_clean) & event_date_clean != "" &
           !is.na(lat_round) & !is.na(lon_round)) %>%
  semi_join(safe_events, by = c("event_date_clean", "lat_round", "lon_round")) %>%
  filter(!has_depth) %>%
  pull(record_key)

cat("Records identified for imputation:", length(records_to_impute), "\n\n")

# Apply: left_join brings in the imputed values, then explicit record_key
# matching (not re-comparing rounded floats) determines which rows actually
# get filled and flagged
echino_wide <- echino_wide %>%
  mutate(
    event_date_clean = coalesce(!!!syms(event_date_col)),
    lat_round = round(best_latitude, 4),
    lon_round = round(best_longitude, 4)
  ) %>%
  left_join(event_depth_lookup, by = c("event_date_clean", "lat_round", "lon_round")) %>%
  mutate(
    best_min_depth = if_else(record_key %in% records_to_impute, imputed_min_depth, best_min_depth),
    best_max_depth = if_else(record_key %in% records_to_impute, imputed_max_depth, best_max_depth),
    depth_imputed_flag = if_else(record_key %in% records_to_impute, TRUE, depth_imputed_flag)
  ) %>%
  select(-imputed_min_depth, -imputed_max_depth, -event_date_clean, -lat_round, -lon_round)

cat("Records with depth_imputed_flag = TRUE after 5f:",
    sum(echino_wide$depth_imputed_flag), "\n")


# =============================================================================
# 5h. DEPTH RECOVERY STEP 2: REEF FLAT / INTERTIDAL TEXT INFERENCE
# =============================================================================
# Records whose verbatimLocality text describes the reef flat / intertidal
# zone are assigned depth = 0m, consistent with the existing verbatimDepth
# == "Intertidal" logic in 5d. Does NOT match text describing horizontal
# distance from shore (e.g. "150 metres off Kirra Beach"), which is
# unrelated to depth and was confirmed absent from the matched set.
# =============================================================================

cat("\n── Section 5g: Depth recovery via reef flat / intertidal text ──\n\n")

reef_flat_pattern <- "(?i)reef flat|intertidal"

reef_flat_candidates <- echino_wide %>%
  filter(!has_depth & !is.na(verbatimLocality)) %>%
  filter(str_detect(verbatimLocality, reef_flat_pattern))

cat("Candidate records (reef flat / intertidal text, no depth):",
    nrow(reef_flat_candidates), "\n")
cat("Distinct verbatimLocality strings matched:",
    n_distinct(reef_flat_candidates$verbatimLocality), "\n\n")

n_before_5g <- sum(echino_wide$has_depth)

echino_wide <- echino_wide %>%
  mutate(
    is_reef_flat_candidate = !has_depth & str_detect(coalesce(verbatimLocality, ""), reef_flat_pattern),
    depth_imputed_flag = if_else(is_reef_flat_candidate, TRUE, depth_imputed_flag),
    best_min_depth = if_else(is_reef_flat_candidate, 0, best_min_depth),
    best_max_depth = if_else(is_reef_flat_candidate, 0, best_max_depth)
  ) %>%
  select(-is_reef_flat_candidate)

cat("Records gained via reef flat / intertidal inference:",
    sum(echino_wide$depth_imputed_flag, na.rm = TRUE) - 645, "\n")

cat("Running total depth_imputed_flag:",
    sum(echino_wide$depth_imputed_flag, na.rm = TRUE), "\n")

# -----------------------------------------------------------------------------
# FINAL RECOMPUTE of all depth-derived columns (post-recovery, authoritative)
# -----------------------------------------------------------------------------

echino_wide <- echino_wide %>%
  mutate(
    depth_range_m = if_else(!is.na(best_min_depth) & !is.na(best_max_depth),
                            best_max_depth - best_min_depth, NA_real_),
    has_depth = !is.na(best_min_depth) | !is.na(best_max_depth)
  ) %>%
  mutate(
    min_zone = case_when(
      is.na(best_min_depth)  ~ NA_character_,
      best_min_depth <= 200  ~ "Continental Shelf",
      best_min_depth <= 1000 ~ "Upper Slope",
      best_min_depth <= 4000 ~ "Lower Slope",
      TRUE                   ~ "Abyssal"
    ),
    max_zone = case_when(
      is.na(best_max_depth)  ~ NA_character_,
      best_max_depth <= 200  ~ "Continental Shelf",
      best_max_depth <= 1000 ~ "Upper Slope",
      best_max_depth <= 4000 ~ "Lower Slope",
      TRUE                   ~ "Abyssal"
    ),
    depth_zone = case_when(
      is.na(min_zone) & is.na(max_zone) ~ "No depth data",
      is.na(min_zone)                   ~ max_zone,
      is.na(max_zone)                   ~ min_zone,
      min_zone == max_zone              ~ min_zone,
      TRUE ~ paste0("Spans ", min_zone, " to ", max_zone)
    )
  ) %>%
  select(-min_zone, -max_zone) %>%
  mutate(
    depth_median = if_else(
      !is.na(best_min_depth) & !is.na(best_max_depth),
      (best_min_depth + best_max_depth) / 2,
      coalesce(best_min_depth, best_max_depth)
    ),
    depth_uncertainty = depth_range_m,
    is_straddler = !is.na(depth_uncertainty) & depth_uncertainty > 100
  )

cat("── Final depth zone distribution (post-recovery) ──\n")
echino_wide %>%
  count(depth_zone, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = Inf)

cat("\nFinal depth_imputed_flag count (event-sibling + reef flat combined):",
    sum(echino_wide$depth_imputed_flag), "\n")
cat("Records with depth_median:", sum(!is.na(echino_wide$depth_median)), "\n")
cat("Records flagged as straddlers (uncertainty > 100m):",
    sum(echino_wide$is_straddler, na.rm = TRUE), "\n")
cat("Records with depth (has_depth):", sum(echino_wide$has_depth), "\n")
cat("Pct with depth:", round(sum(echino_wide$has_depth) / nrow(echino_wide) * 100, 1), "%\n")
cat("Records remaining:", nrow(echino_wide), "\n\n")

cat("Section 5 complete.\n")
cat("Final depth coverage:", sum(echino_wide$has_depth), "of", nrow(echino_wide),
    sprintf("(%.1f%%), including %d records (%.1f%% of total) recovered via\n",
            100 * mean(echino_wide$has_depth),
            sum(echino_wide$depth_imputed_flag),
            100 * mean(echino_wide$depth_imputed_flag)))
cat("event-sibling imputation or reef flat/intertidal text inference (see depth_imputed_flag).\n")

# =============================================================================
# SECTION 6: NAME NORMALISATION + WoRMS RESOLUTION
#
# This is the consolidated, corrected version of Section 6, incorporating
# every fix made during the validation session:
#   6.0  Non-echinoderm contaminant removal (42 records: Asterina eupomatiae fungus, Anoplura lice, Rhodostoma
#        Mollusca genus, ANWC mammals - all exact-match only)
#   6a   Name normalisation function + regression tests
#   6b   Apply normalisation + morphospecies split to full dataset
#   6c   WoRMS resolution - automated passes (exact match, subgenus fix,
#        author-fragment fixes, fuzzy/taxamatch)
#   6d   Manual correction lookup (18 names requiring taxonomic expertise)
#   6e   Final accepted_name + taxonomic_resolution_level assignment
#   6f   Validation (phylum mismatch check, spot check)
#   6g   Sanity check on non-conflicted records
#   6h   FULL WoRMS validation of every accepted_name, including names that
#        were never conflicted and so never queried in 6c. Fixes the
#        queried_name tracking bug and the mixed-status exclusion bug.
#        Retains Hemiaster (extinct genus, but basisOfRecord confirms it is
#        a genuine PreservedSpecimen, not a fossil - see Section 6.0 notes).
#
# =============================================================================

cat("── Section 6: Name normalisation + WoRMS resolution ──\n\n")

# =============================================================================
# SECTION 6.0: NON-ECHINODERM CONTAMINANT REMOVAL
# =============================================================================
# All removals use EXACT record_key matching or exact institutionCode +
# collectionCode equality - never broad text-pattern matching. An earlier
# attempt using str_detect() on collectionCode accidentally matched 72
# unrelated legitimate AM echinoderm records sharing the same generic
# collection code; that approach was discarded for this reason.
#
# Removal rationale, by record:
#   - "Asterina eupomatiae" (1 record): a fungus, not the sea star genus
#     Asterina.
#   - "Anoplura" (2 records): sucking lice (Arthropoda).
#   - "Rhodostoma" (1 record): confirmed Mollusca genus via direct WoRMS
#     phylum lookup, not Echinodermata.
#   - ANWC "Mammals"-coded records (38): confirmed genuine ANWC mammal
#     voucher specimens via GBIF cross-check, despite carrying populated
#     Echinodermata taxonomic fields.

cat("── Section 6.0: Non-echinoderm contaminant removal ──\n\n")

non_echinoderm_contaminant_keys <- c(
  "S:S-FUNGI:F12519", # Asterina
  "6E154857-AD6C-403F-B789-4A9C79F093FA", # Anoplura
  "9439C5D6-927C-43C1-88A0-B70E45EB2501", # Anoplura
  "0703DE68-1792-4DA8-A52F-907C34B27AC5"  # Rhodostoma
)
# Mammals
anwc_mammal_contaminant_keys <- echino_wide %>%
  filter(institutionCode == "ANWC" & collectionCode == "Mammals") %>%
  pull(record_key)

all_contaminant_keys <- c(non_echinoderm_contaminant_keys, anwc_mammal_contaminant_keys)

cat("Non-echinoderm contaminant records to remove:", length(all_contaminant_keys), "\n")
cat("  - Asterina eupomatiae (fungus):", 1, "\n")
cat("  - Anoplura (lice):", 2, "\n")
cat("  - Rhodostoma (Mollusca genus, confirmed via direct WoRMS phylum check):", 1, "\n")
cat("  - ANWC Mammals (mislabelled mammal specimens):", length(anwc_mammal_contaminant_keys), "\n\n")

n_before <- nrow(echino_wide)
echino_wide <- echino_wide %>% filter(!record_key %in% all_contaminant_keys)
cat("Records removed:", n_before - nrow(echino_wide), "\n")
cat("Records remaining:", nrow(echino_wide), "\n\n")


# =============================================================================
# SECTION 6a: NAME NORMALISATION FUNCTION
# =============================================================================

cat("== Section 6a: Name normalisation function ==\n\n")

sciname_cols <- names(echino_wide)[str_detect(names(echino_wide), "^scientificName__")]
cat("Scientific name source columns found:", length(sciname_cols), "\n")

high_rank_terms <- c(
  "echinodermata", "asterozoa", "echinozoa",
  "asteroidea", "echinoidea", "holothuroidea", "ophiuroidea", "crinoidea"
)

strip_author <- function(x) {
  x <- str_trim(x)
  x <- str_to_lower(x)
  
  # Pattern 3: "name (author, year), year2"
  x <- str_remove(x, "\\s*\\([^)]*\\d{4}[a-z]?[^)]*\\)\\s*,\\s*\\d{4}[a-z]?\\s*$")
  
  # Pattern 2: "name (author, year)" or "name (author)"
  x <- str_remove(x, "\\s*\\([^)]*\\)\\s*$")
  
  # Pattern 1: bare trailing author(s) + year
  x <- str_remove(
    x,
    "\\s+[a-z\u00e0-\u00ff.\\-']+((\\s*[,&]\\s*[a-z\u00e0-\u00ff.\\-']+)*),?\\s*\\d{4}[a-z]?\\s*$"
  )
  x <- str_trim(x)
  
  collapse_one <- function(s) {
    if (is.na(s) || s == "") return(s)
    first_word <- str_extract(s, "^\\S+")
    if (is.na(first_word) || !(first_word %in% high_rank_terms)) return(s)
    words <- str_split(s, "\\s+")[[1]]
    last_word <- words[length(words)]
    is_bare_morpho <- str_detect(last_word, "^sp\\.?\\d*[a-z]?$")
    if (is_bare_morpho && length(words) >= 2) {
      paste(words[length(words) - 1], words[length(words)])
    } else {
      last_word
    }
  }
  x <- vapply(x, collapse_one, character(1), USE.NAMES = FALSE)
  x <- str_trim(x)
  
  x <- if_else(
    is.na(x) | x == "",
    NA_character_,
    paste0(str_to_upper(str_sub(x, 1, 1)), str_sub(x, 2))
  )
  x
}

strip_stray_author_fragment <- function(x) {
  str_remove(x, "(?i)\\s+(ah|van|de|von|der|du|da)$")
}

strip_brackets <- function(x) str_replace_all(x, "\\[([^\\]]+)\\]", "\\1")

extract_morphospecies <- function(x) {
  str_extract(x, "(?i)\\bsp\\.?\\s*\\d*[a-z]?$|\\bspp\\b$|\\bgp\\.?\\s*\\d+$")
}

strip_morphospecies_tag <- function(x) {
  tag <- extract_morphospecies(x)
  cleaned <- if_else(
    !is.na(tag),
    str_trim(str_remove(x, "(?i)\\s*\\bsp\\.?\\s*\\d*[a-z]?$|\\s*\\bspp\\b$|\\s*\\bgp\\.?\\s*\\d+$")),
    x
  )
  cleaned <- str_trim(str_remove(cleaned, "(?i)\\s+tpe$"))
  cleaned <- if_else(
    is.na(cleaned) | cleaned == "",
    NA_character_,
    paste0(str_to_upper(str_sub(cleaned, 1, 1)), str_sub(cleaned, 2))
  )
  cleaned
}

capitalise_subgenus <- function(x) {
  str_replace(x, "\\(([a-z\u00e0-\u00ff]+)\\)", function(m) {
    inner <- str_extract(m, "[a-z\u00e0-\u00ff]+")
    paste0("(", str_to_upper(str_sub(inner, 1, 1)), str_sub(inner, 2), ")")
  })
}

regression_tests <- tribble(
  ~original, ~expected,
  "HOLOTHUROIDEA", "Holothuroidea",
  "Holothuroidea", "Holothuroidea",
  "holothuroidea", "Holothuroidea",
  "Ophiomyces", "Ophiomyces",
  "Ophiomyces Lyman, 1869", "Ophiomyces",
  "Anthenea crassa", "Anthenea crassa",
  "Anthenea crassa H.L.Clark, 1938", "Anthenea crassa",
  "Acanthaster planci", "Acanthaster planci",
  "Acanthaster planci (Linnaeus, 1758)", "Acanthaster planci",
  "Luidia maculata M\u00fcller & Troschel, 1851", "Luidia maculata",
  "Echinodermata Holothuroidea", "Holothuroidea",
  "Echinodermata Asterozoa Asteroidea Valvatida Goniasteridae", "Goniasteridae",
  "Echinodermata Asterozoa Ophiuroidea Ophiurida sp4", "Ophiurida sp4",
  "Ophiomusium sp1", "Ophiomusium sp1",
  "ophiothrix tpe sp. 4", "Ophiothrix tpe sp. 4",
  "Ophiuroid spE", "Ophiuroid spe",
  "echinometra sp. a", "Echinometra sp. a",
  "Comanthus gisleni Rowe, Hoggett, 2020", "Comanthus gisleni",
  "Jacksonaster depressum (L.Agassiz, 1841), 2020", "Jacksonaster depressum",
  "Ailsastra O'Loughlin & Rowe, 2005", "Ailsastra"
)

cat("Running regression test suite (", nrow(regression_tests), "cases)...\n")
regression_check <- regression_tests %>%
  mutate(got = strip_author(original), pass = got == expected)
cat("Passed:", sum(regression_check$pass), "of", nrow(regression_check), "\n")
if (any(!regression_check$pass)) {
  cat("FAILURES:\n")
  print(regression_check %>% filter(!pass))
}


# =============================================================================
# SECTION 6b: APPLY NORMALISATION + MORPHOSPECIES SPLIT TO FULL DATASET
# =============================================================================

cat("\n== Section 6b: Applying normalisation to full dataset ==\n\n")

echino_wide <- echino_wide %>%
  mutate(across(
    all_of(sciname_cols),
    ~ strip_stray_author_fragment(strip_author(strip_brackets(.x))),
    .names = "{.col}_norm"
  ))

norm_cols <- names(echino_wide)[str_detect(names(echino_wide), "^scientificName__.*_norm$")]

echino_wide <- echino_wide %>%
  mutate(across(
    all_of(norm_cols),
    list(
      morphotag = ~ extract_morphospecies(.x),
      clean     = ~ strip_morphospecies_tag(.x)
    ),
    .names = "{.col}_{.fn}"
  ))

clean_cols     <- names(echino_wide)[str_detect(names(echino_wide), "^scientificName__.*_norm_clean$")]
morphotag_cols <- names(echino_wide)[str_detect(names(echino_wide), "^scientificName__.*_norm_morphotag$")]

clean_mat <- echino_wide %>% select(all_of(clean_cols)) %>% as.matrix()
name_conflict_real <- apply(clean_mat, 1, function(row) {
  vals <- unique(na.omit(row)); vals <- vals[vals != ""]; length(vals) > 1
})

echino_wide <- echino_wide %>%
  mutate(
    name_conflict_real = name_conflict_real,
    is_morphospecies = apply(select(., all_of(morphotag_cols)), 1, function(r) any(!is.na(r)))
  )

cat("Name conflicts BEFORE normalisation:", sum(echino_wide$name_conflict, na.rm = TRUE), "\n")
cat("Name conflicts AFTER normalisation: ", sum(name_conflict_real), "\n")
cat("Records flagged as morphospecies:   ", sum(echino_wide$is_morphospecies, na.rm = TRUE), "\n\n")

all_names_to_check <- echino_wide %>%
  filter(name_conflict_real) %>%
  select(all_of(clean_cols)) %>%
  unlist() %>% na.omit() %>% unique() %>% .[. != ""] %>% sort()

cat("Unique names requiring WoRMS resolution:", length(all_names_to_check), "\n")


# =============================================================================
# SECTION 6c: WoRMS RESOLUTION (automated passes)
# =============================================================================

cat("\n== Section 6c: WoRMS resolution ==\n\n")

# --- FIXED helper: queried_name tracked explicitly through the response,
# never reconstructed from scientificname afterward. This was the bug that
# caused names to be silently dropped/duplicated in earlier runs. ---
query_worms_exact <- function(names_vec) {
  if (length(names_vec) == 0) return(tibble())
  batches <- split(names_vec, ceiling(seq_along(names_vec) / 50))
  map_dfr(batches, function(b) {
    res <- tryCatch(
      wm_records_names(name = b, fuzzy = FALSE),
      error = function(e) NULL
    )
    if (is.null(res)) {
      return(map_dfr(b, function(nm) {
        r <- tryCatch(wm_records_names(name = nm, fuzzy = FALSE)[[1]],
                      error = function(e2) tibble(scientificname = NA_character_, status = "not found"))
        r %>% mutate(queried_name = nm)
      }))
    }
    map2_dfr(b, res, function(nm, r) {
      if (is.null(r) || nrow(r) == 0) {
        tibble(queried_name = nm, scientificname = NA_character_, status = "not found")
      } else {
        r %>% mutate(queried_name = nm)
      }
    })
  })
}

# Note: taxamatch's per-request size limit is much smaller than exact-match's
# (~10-20 names); batches of >~20 can trigger a 400 Bad Request even though
# each name individually is valid.
query_worms_taxamatch <- function(names_vec, batch_size = 10) {
  if (length(names_vec) == 0) return(tibble())
  batches <- split(names_vec, ceiling(seq_along(names_vec) / batch_size))
  map_dfr(batches, function(b) {
    results <- tryCatch(
      wm_records_taxamatch(name = b),
      error = function(e) {
        cat("  taxamatch batch failed at size", length(b), ":", conditionMessage(e), "\n")
        NULL
      }
    )
    if (is.null(results)) {
      return(tibble(queried_name = b, scientificname = NA_character_,
                    status = "not found", rank = NA_character_))
    }
    map2_dfr(b, results, function(orig, res) {
      if (is.null(res) || nrow(res) == 0) {
        tibble(queried_name = orig, scientificname = NA_character_, status = "not found", rank = NA_character_)
      } else {
        res %>% mutate(queried_name = orig) %>% select(queried_name, scientificname, status, rank)
      }
    })
  })
}

status_priority <- c(
  "accepted" = 1, "alternative representation" = 2,
  "unaccepted" = 3, "superseded combination" = 3, "superseded rank" = 3,
  "incorrect grammatical agreement of specific epithet" = 3,
  "junior subjective synonym" = 3,
  "nomen dubium" = 50, "nomen nudum" = 50, "taxon inquirendum" = 50
)

# --- Pass 1: exact match ---
cat("Pass 1: exact match,", length(all_names_to_check), "names\n")
pass1 <- query_worms_exact(all_names_to_check)
names_not_found_1 <- pass1 %>%
  filter(is.na(scientificname) | status == "not found") %>%
  pull(queried_name) %>% unique()
names_not_found_1 <- union(names_not_found_1, setdiff(all_names_to_check, pass1$queried_name))
cat("  Resolved:", length(all_names_to_check) - length(names_not_found_1),
    "| Remaining:", length(names_not_found_1), "\n\n")

# --- Pass 2: subgenus capitalisation fix ---
cat("Pass 2: subgenus capitalisation fix,", length(names_not_found_1), "names\n")
names_p2 <- capitalise_subgenus(names_not_found_1) %>% unique() %>% sort()
pass2 <- query_worms_exact(names_p2)
p2_resolved_map <- tibble(
  original_clean_name = names_not_found_1,
  fixed_name = capitalise_subgenus(names_not_found_1)
) %>%
  filter(fixed_name %in% (pass2 %>% filter(!is.na(scientificname)) %>% pull(scientificname)))
names_not_found_2 <- setdiff(names_not_found_1, p2_resolved_map$original_clean_name)
cat("  Resolved:", nrow(p2_resolved_map), "| Remaining:", length(names_not_found_2), "\n\n")

# --- Pass 3: H.L. Clark / hubert lyman / Wyville Thomson author-fragment fixes ---
cat("Pass 3: known author-fragment patterns\n")
strip_hl_fragment <- function(x) str_trim(str_remove(x, "(?i)\\s+h\\.?\\s*l\\.?\\s*(clark)?$"))
strip_hubert_lyman <- function(x) str_remove(x, "(?i)\\s+hubert\\s+lyman(\\s+clark)?$")
strip_wyville <- function(x) str_remove(x, "(?i)\\s+wyville$")
author_fragment_fix <- function(x) {
  x <- strip_hl_fragment(x)
  x <- strip_hubert_lyman(x)
  x <- strip_wyville(x)
  str_trim(x)
}
names_p3_fixed <- capitalise_subgenus(author_fragment_fix(names_not_found_2)) %>% unique() %>% sort()
pass3 <- query_worms_exact(names_p3_fixed)
p3_resolved_map <- tibble(
  original_clean_name = names_not_found_2,
  fixed_name = capitalise_subgenus(author_fragment_fix(names_not_found_2))
) %>%
  filter(fixed_name %in% (pass3 %>% filter(!is.na(scientificname)) %>% pull(scientificname)))
names_not_found_3 <- setdiff(names_not_found_2, p3_resolved_map$original_clean_name)
cat("  Resolved:", nrow(p3_resolved_map), "| Remaining:", length(names_not_found_3), "\n\n")

# --- Pass 4: semicolon-joined duplicate names ---
cat("Pass 4: semicolon-joined duplicate names\n")
semicolon_names <- names_not_found_3[str_detect(names_not_found_3, ";")]
fix_semicolon_dup <- function(x) str_trim(str_split(x, ";")[[1]][1]) %>% strip_author()
names_p4_fixed <- if (length(semicolon_names) > 0) {
  capitalise_subgenus(sapply(semicolon_names, fix_semicolon_dup, USE.NAMES = FALSE))
} else character(0)
pass4 <- query_worms_exact(unique(names_p4_fixed))
p4_resolved_map <- tibble(
  original_clean_name = semicolon_names,
  fixed_name = names_p4_fixed
) %>%
  filter(fixed_name %in% (if (nrow(pass4) > 0) pass4 %>% filter(!is.na(scientificname)) %>% pull(scientificname) else character(0)))
names_not_found_4 <- setdiff(names_not_found_3, p4_resolved_map$original_clean_name)
cat("  Resolved:", nrow(p4_resolved_map), "| Remaining:", length(names_not_found_4), "\n\n")

# --- Pass 5: fuzzy taxamatch on remaining clean-looking binomials ---
cat("Pass 5: fuzzy taxamatch on remaining names\n")
fix_remaining_initials <- function(x) {
  x <- str_remove(x, "(?i)\\s+[a-z]\\.\\s*[a-z]?\\.?\\s*$")
  x <- str_remove(x, "(?i)\\s+[a-z]\\s+[a-z]$")
  x <- str_remove(x, "(?i),?\\s+\\d{4}\\s*,?\\s*\\d{4}$")
  str_trim(x)
}
fix_and_author      <- function(x) str_remove(x, "(?i)\\s+[a-z\u00e0-\u00ff.\\-']+\\s+and\\s+[a-z\u00e0-\u00ff.\\-']+$")
fix_ohara            <- function(x) str_remove(x, "(?i)\\s+o'hara$")
fix_l_agassiz        <- function(x) {
  x <- str_remove(x, "(?i)\\s+l\\.\\s*agassiz\\s+in\\s+l\\.$")
  x <- str_remove(x, "(?i)\\s+l\\.\\s*agassiz,?\\s*\\d{4}[a-z]?$")
  x <- str_remove(x, "(?i)\\s+l\\.$")
  str_trim(x)
}
fix_holothuria_subgenus <- function(x) str_replace(x, "^(Holothuria)\\s+(\\w+)\\s+(.+)$", "\\1 (\\2) \\3")
fix_repeated_genus <- function(x) {
  words <- str_split(x, "\\s+")[[1]]
  if (length(words) >= 2 && str_to_lower(words[1]) == str_to_lower(words[2])) {
    paste0(words[1], " (", words[2], ") ", paste(words[-(1:2)], collapse = " "))
  } else x
}
fix_comma_joined_names <- function(x) str_trim(str_split(x, ",")[[1]][1])

names_p5_pre <- names_not_found_4 %>%
  fix_remaining_initials() %>%
  fix_and_author() %>%
  fix_ohara() %>%
  fix_l_agassiz()

names_p5_pre <- if_else(
  str_detect(names_not_found_4, "^Holothuria\\s+\\w+\\s+\\w+$") & !str_detect(names_not_found_4, "\\("),
  fix_holothuria_subgenus(names_not_found_4),
  names_p5_pre
)
names_p5_pre <- if_else(
  str_detect(names_p5_pre, "^(\\w+)\\s+\\1\\s"),
  sapply(names_p5_pre, fix_repeated_genus),
  names_p5_pre
)
names_p5_pre <- if_else(
  str_detect(names_p5_pre, ","),
  sapply(names_p5_pre, fix_comma_joined_names),
  names_p5_pre
)

names_p5_fixed <- capitalise_subgenus(names_p5_pre) %>% unique()
pass5 <- query_worms_taxamatch(names_p5_fixed)

p5_resolved_map <- tibble(
  original_clean_name = names_not_found_4,
  fixed_name = capitalise_subgenus(names_p5_pre)
) %>%
  left_join(pass5 %>% select(queried_name, status) %>% distinct(),
            by = c("fixed_name" = "queried_name")) %>%
  filter(!is.na(status) & status != "not found")

names_not_found_5 <- setdiff(names_not_found_4, p5_resolved_map$original_clean_name)
cat("  Resolved:", nrow(p5_resolved_map), "| Remaining:", length(names_not_found_5), "\n\n")

cat("Names requiring manual taxonomic review:", length(names_not_found_5), "\n")
print(names_not_found_5)


# =============================================================================
# SECTION 6d: MANUAL CORRECTIONS
# =============================================================================

cat("\n== Section 6d: Manual taxonomic corrections ==\n\n")

manual_corrections <- tribble(
  ~queried_name, ~resolved_name, ~worms_rank, ~notes,
  "Acanthophiothrix", "Ophiothrix (Acanthophiothrix)", "Subgenus",
  "Subgenus only, genus dropped in source",
  "Alleocomatella", "Alloeocomatella", "Genus",
  "Likely misspelling",
  "Alloeocomatella pectinifer", "Alloeocomatella pectinifera", "Species",
  "Grammatical gender agreement",
  "Amphiodia loripes", "Amphiodia (Amphispina) loripes", "Species",
  "Subgenus dropped in source",
  "Amphioplus lobata", "Amphioplus (Amphioplus) lobata", "Species",
  "Subgenus dropped in source",
  "Amphioplus ochroleuca", "Amphioplus (Amphichilus) ochroleuca", "Species",
  "Subgenus dropped in source",
  "Anthenea acutus", "Peltaster placenta", "Species",
  "Anthenea acuta unaccepted (homonym, A.M. Clark 1993); accepted name is Peltaster placenta (Muller & Troschel, 1842)",
  "Clarkcomanthus sp. 1 of summers et", "Clarkcomanthus", "Genus",
  "Morphospecies tag references a project citation (Summers et al.), genus retained",
  "Comanthus timorensis xanthum h l clark", "Comanthus parvicirrus", "Species",
  "Xanthum is a varietal name; Comanthus timorensis unaccepted, accepted name is C. parvicirrus (Muller, 1841)",
  "Comatula (comatula) pectinata purpurea", "Comatula (Comatula) pectinata", "Species",
  "Two species names present (pectinata vs purpurea varietal); pectinata appears in majority of sources, chosen as best identification",
  "Holothuria trapedza hubert lyman clark", "Holothuria", "Genus",
  "Species not found in WoRMS; other sources show genus only, kept as genus",
  "Jacksonaster depressus s\u00e1nchez", "Jacksonaster depressum", "Species",
  "S\u00e1nchez is part of author citation; also a spelling correction (depressus -> depressum)",
  "Ophiuroid", "Ophiuroidea", "Class",
  "Class-level identification",
  "Holothuria philippinensis", "Holothuria (Stauropora) fuscocinerea", "Species",
  "Philippinensis is a variety (Domantay, 1933) of H. fuscocinerea; WoRMS lists it as nomen dubium under this parent",
  "Pendekaplectana nigra", "Polyplectana nigra", "Species",
  "Unaccepted name",
  "Temnotrema decorum", "Temnotrema", "Genus",
  "Species not found in WoRMS, kept as genus",
  "Comanthus maculosa", "Comanthus", "Genus",
  "Species not found in WoRMS; other sources show genus only, kept as genus",
  "Oxycomanthus carpenteri", "Clarkcomanthus", "Genus",
  "Species not found under this genus; synonymy points to Clarkcomanthus",
  "Fromia japonica", "Fromia monilis", "Species",
  "Unaccepted name"
) %>%
  mutate(worms_status = "manual_correction")

write_csv(manual_corrections %>% select(-worms_rank, -worms_status) %>%
            rename(original_clean_name = queried_name, suggested_correction = resolved_name),
          "name_corrections_for_review.csv")

cat("Manual corrections applied:", nrow(manual_corrections), "\n")
cat("Saved audit trail: name_corrections_for_review.csv\n")


# =============================================================================
# SECTION 6e: FINAL accepted_name + taxonomic_resolution_level ASSIGNMENT
# =============================================================================

cat("\n== Section 6e: Final name + rank resolution ==\n\n")

exact_match_results <- bind_rows(
  pass1 %>% select(queried_name, scientificname, status, valid_name = any_of("valid_name"), rank = any_of("rank")),
  pass2 %>% mutate(queried_name = scientificname) %>%
    select(queried_name, scientificname, status, valid_name = any_of("valid_name"), rank = any_of("rank")),
  pass3 %>% mutate(queried_name = scientificname) %>%
    select(queried_name, scientificname, status, valid_name = any_of("valid_name"), rank = any_of("rank")),
  pass4 %>% mutate(queried_name = scientificname) %>%
    select(queried_name, scientificname, status, valid_name = any_of("valid_name"), rank = any_of("rank"))
) %>%
  bind_rows(
    p2_resolved_map %>% rename(queried_name = original_clean_name, scientificname = fixed_name) %>%
      left_join(pass2 %>% select(scientificname, status, valid_name = any_of("valid_name"), rank = any_of("rank")),
                by = "scientificname"),
    p3_resolved_map %>% rename(queried_name = original_clean_name, scientificname = fixed_name) %>%
      left_join(pass3 %>% select(scientificname, status, valid_name = any_of("valid_name"), rank = any_of("rank")),
                by = "scientificname"),
    p4_resolved_map %>% rename(queried_name = original_clean_name, scientificname = fixed_name) %>%
      left_join(pass4 %>% select(scientificname, status, valid_name = any_of("valid_name"), rank = any_of("rank")),
                by = "scientificname")
  )

taxamatch_results <- p5_resolved_map %>%
  select(original_clean_name, fixed_name) %>%
  rename(queried_name = original_clean_name) %>%
  left_join(
    pass5 %>% select(queried_name, scientificname, status, rank) %>% distinct(),
    by = c("fixed_name" = "queried_name"),
    relationship = "many-to-many"
  ) %>%
  select(-fixed_name)

all_worms_results <- bind_rows(
  exact_match_results %>% select(queried_name, scientificname, status, rank, valid_name),
  taxamatch_results %>% mutate(valid_name = NA_character_) %>%
    select(queried_name, scientificname, status, rank, valid_name)
) %>%
  filter(!is.na(queried_name))

worms_best_match <- all_worms_results %>%
  mutate(status_rank = coalesce(status_priority[status], 10)) %>%
  group_by(queried_name) %>%
  slice_min(status_rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(resolved_name = if_else(status == "accepted", scientificname, coalesce(valid_name, scientificname))) %>%
  select(queried_name, resolved_name, worms_status = status, worms_rank = rank)

final_name_lookup <- bind_rows(
  manual_corrections %>% select(queried_name, resolved_name, worms_status, worms_rank),
  worms_best_match %>% filter(!queried_name %in% manual_corrections$queried_name)
)

cat("Final consolidated name lookup:", nrow(final_name_lookup), "entries\n\n")

source_trust_rank <- c(
  "AM_direct"          = 1,
  "iDigBio"            = 2,
  "OBIS"               = 3,
  "ALA_GBRSBD"         = 4,
  "ALA_echinodermata"  = 5,
  "ALA_other"          = 6,
  "CSIRO_GBRSBD"       = 7,
  "CSIRO_QM"           = 8,
  "GBIF_QM"            = 9,
  "GBIF_echinodermata" = 10,
  "OZCAM"              = 11,
  "MTQ_CIDARIS"        = 12,
  "CIDARIS_QMT"        = 13
)

name_resolution <- echino_wide %>%
  select(record_key, all_of(clean_cols)) %>%
  pivot_longer(
    cols      = all_of(clean_cols),
    names_to  = "source_col",
    values_to = "queried_name"
  ) %>%
  filter(!is.na(queried_name) & queried_name != "") %>%
  
  # Extract source key from column name
  # e.g. "scientificName__AM_direct_norm_clean" → "AM_direct"
  mutate(
    source = str_remove(str_remove(source_col, "^scientificName__"), "_norm_clean$")
  ) %>%
  
  left_join(final_name_lookup, by = "queried_name") %>%
  mutate(
    final_name = coalesce(resolved_name, queried_name),
    
    # CRITERION 1: WoRMS status quality — unchanged from before
    status_rank = case_when(
      worms_status == "manual_correction"          ~ 0,
      worms_status == "accepted"                   ~ 1,
      worms_status == "alternative representation" ~ 2,
      is.na(worms_status)                          ~ 5,
      TRUE                                         ~ 3
    ),
    
    # CRITERION 2: Taxonomic fineness — finer identification wins ties
    # This is what was silently lost before: Species was losing to Class
    # when both had worms_status == "accepted"
    fineness_rank = case_when(
      worms_rank %in% c("Subspecies", "Variety")               ~ 1,
      worms_rank == "Species"                                   ~ 2,
      worms_rank %in% c("Subgenus", "Genus")                   ~ 3,
      worms_rank %in% c("Subfamily", "Family", "Superfamily")  ~ 4,
      worms_rank %in% c("Suborder", "Order", "Superorder")     ~ 5,
      worms_rank %in% c("Subclass", "Class", "Superclass")     ~ 6,
      worms_rank == "Phylum"                                    ~ 7,
      TRUE                                                      ~ 8   # unknown rank
    ),
    
    # CRITERION 3: Source trust — only kicks in if status AND fineness tie
    # e.g. two sources both say "Hacelia" (Genus, accepted) → AM_direct wins
    src_rank = coalesce(source_trust_rank[source], 99)
  ) %>%
  
  # Sort so the best candidate rises to row 1 within each record_key group.
  # arrange() + first() is clearer and safer than which.min() on multiple columns.
  arrange(status_rank, fineness_rank, src_rank) %>%
  
  group_by(record_key) %>%
  summarise(
    accepted_name        = first(final_name),
    accepted_name_rank   = first(worms_rank),
    winning_source       = first(source),     # audit trail: which source won
    accepted_name_worms  = any(!is.na(worms_status) & worms_status != "manual_correction"),
    accepted_name_manual = any(worms_status == "manual_correction", na.rm = TRUE),
    .groups = "drop"
  )

specific_epithet_cols <- names(echino_wide)[str_detect(names(echino_wide), "^specificEpithet__")]
genus_cols  <- names(echino_wide)[str_detect(names(echino_wide), "^genus__")]
family_cols <- names(echino_wide)[str_detect(names(echino_wide), "^family__")]
order_cols  <- names(echino_wide)[str_detect(names(echino_wide), "^order__")]
class_cols  <- names(echino_wide)[str_detect(names(echino_wide), "^class__")]
phylum_cols <- names(echino_wide)[str_detect(names(echino_wide), "^phylum__")]

has_any <- function(df, cols) rowSums(!is.na(df[cols]) & df[cols] != "") > 0

fallback_rank <- echino_wide %>%
  mutate(
    fallback_rank = case_when(
      has_any(., specific_epithet_cols) ~ "Species",
      has_any(., genus_cols)            ~ "Genus",
      has_any(., family_cols)           ~ "Family",
      has_any(., order_cols)            ~ "Order",
      has_any(., class_cols)            ~ "Class",
      has_any(., phylum_cols)           ~ "Phylum",
      TRUE                              ~ NA_character_
    )
  ) %>%
  select(record_key, fallback_rank)

echino_wide <- echino_wide %>%
  select(-any_of(c("accepted_name", "accepted_name_rank", "accepted_name_worms",
                   "accepted_name_manual", "fallback_rank",
                   "taxonomic_resolution_level"))) %>%
  left_join(name_resolution, by = "record_key") %>%
  left_join(fallback_rank, by = "record_key")

echino_wide <- echino_wide %>%
  mutate(taxonomic_resolution_level = coalesce(accepted_name_rank, fallback_rank)) %>%
  select(-accepted_name_rank)

cat("Records with accepted_name assigned:", sum(!is.na(echino_wide$accepted_name)),
    "of", nrow(echino_wide), "\n")
cat("Resolved via WoRMS:  ", sum(echino_wide$accepted_name_worms, na.rm = TRUE), "\n")
cat("Resolved via manual: ", sum(echino_wide$accepted_name_manual, na.rm = TRUE), "\n\n")

cat("Taxonomic resolution level (after Section 6e):\n")
echino_wide %>%
  count(taxonomic_resolution_level, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = Inf)


# =============================================================================
# SECTION 6f: VALIDATION
# =============================================================================

cat("\n== Section 6f: Validation ==\n\n")

cat("Records with missing accepted_name:", sum(is.na(echino_wide$accepted_name)), "\n")

distinct_names <- unique(na.omit(echino_wide$accepted_name))
spot_check_n <- min(30, length(distinct_names))
spot_sample <- sample(distinct_names, spot_check_n)
spot_check_results <- query_worms_exact(spot_sample)
cat("\nSpot check (", spot_check_n, "random accepted_name values):\n")
spot_check_results %>% count(status, sort = TRUE) %>% print()

cat("\nChecking", length(distinct_names), "distinct accepted_name values for phylum mismatches (this takes a while)...\n")
phylum_check <- map_dfr(distinct_names, function(nm) {
  result <- tryCatch(wm_records_names(name = nm, fuzzy = FALSE)[[1]], error = function(e) NULL)
  if (is.null(result) || nrow(result) == 0) {
    return(tibble(name = nm, n_phyla = 0, accepted_phylum = NA_character_, all_phyla = NA_character_))
  }
  result %>% summarise(
    name = nm,
    n_phyla = n_distinct(phylum, na.rm = TRUE),
    accepted_phylum = phylum[status == "accepted"][1],
    all_phyla = paste(unique(na.omit(phylum)), collapse = "; ")
  )
})
phylum_mismatch <- phylum_check %>% filter(!is.na(accepted_phylum) & accepted_phylum != "Echinodermata")
cat("Names where the ACCEPTED WoRMS phylum is NOT Echinodermata:", nrow(phylum_mismatch), "\n")
if (nrow(phylum_mismatch) > 0) {
  cat("WARNING: review these manually before treating the dataset as final.\n")
  print(phylum_mismatch, n = Inf)
}

worms_status_reference <- tribble(
  ~status, ~meaning, ~priority_used,
  "accepted", "Currently valid taxonomic name", "1 (highest)",
  "alternative representation", "A valid alternative way of writing the accepted name", "2",
  "unaccepted", "A synonym; a different, currently-valid name exists", "3",
  "superseded combination", "Genus placement has changed since publication", "3",
  "superseded rank", "Taxonomic rank has changed since publication", "3",
  "incorrect grammatical agreement of specific epithet", "Minor gender-ending error", "3",
  "junior subjective synonym", "Later name judged to refer to same species as an earlier one", "3",
  "nomen dubium", "Uncertain application - not confidently resolvable", "50 (not auto-resolved)",
  "nomen nudum", "Name published without valid description", "50 (not auto-resolved)",
  "taxon inquirendum", "Taxon requiring further investigation", "50 (not auto-resolved)",
  "manual_correction", "Resolved by direct taxonomic review, not automated WoRMS match", "0 (overrides WoRMS)"
)
write_csv(worms_status_reference, "worms_status_reference.csv")

final_status_distribution <- final_name_lookup %>%
  count(worms_status, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1))
write_csv(final_status_distribution, "worms_status_distribution.csv")
cat("\nFinal status distribution (Section 6e names only):\n")
print(final_status_distribution, n = Inf)


# =============================================================================
# SECTION 6g: SANITY CHECK ON NON-CONFLICTED RECORDS
# =============================================================================

cat("\n== Section 6g: Validation of non-conflicted records ==\n\n")

non_conflict_clean <- echino_wide %>%
  filter(!name_conflict_real) %>%
  select(record_key, all_of(clean_cols)) %>%
  pivot_longer(cols = all_of(clean_cols), names_to = "src", values_to = "val") %>%
  filter(!is.na(val))

cat("Non-conflicted records:", n_distinct(non_conflict_clean$record_key), "\n")

suspicious_names <- non_conflict_clean %>%
  distinct(val) %>%
  filter(
    str_detect(val, "[0-9]") |
      str_detect(val, "[,;]") |
      str_detect(val, "^[a-z]") |
      str_detect(val, "\\s{2,}") |
      nchar(val) <= 3
  )
cat("Suspicious non-conflicted names found:", nrow(suspicious_names), "\n")
if (nrow(suspicious_names) > 0) print(suspicious_names, n = 50)

cat("\nRecords missing taxonomic_resolution_level:",
    sum(is.na(echino_wide$taxonomic_resolution_level)), "\n")


# =============================================================================
# SECTION 6h: FULL WoRMS VALIDATION OF EVERY accepted_name
# =============================================================================
# Closes the gap left by Sections 6c-6f, which only query WoRMS for names
# that were CONFLICTED across sources. This section validates EVERY
# distinct accepted_name regardless.
#
# NOTE: "Rhodostoma" (confirmed Mollusca genus) was originally discovered
# via this full-validation pass since it carried no recorded WoRMS status
# and so was invisible to the phylum-mismatch check in 6f. It is now
# removed in Section 6.0 above, alongside the other non-echinoderm
# contaminants, rather than here.
# =============================================================================

cat("\n== Section 6h: Full WoRMS validation of every accepted_name ==\n\n")

all_distinct_names <- echino_wide %>%
  filter(!is.na(accepted_name)) %>%
  distinct(accepted_name) %>%
  pull(accepted_name) %>%
  sort()

already_in_lookup <- final_name_lookup$queried_name
names_to_check <- setdiff(all_distinct_names, already_in_lookup)

cat("Already in final_name_lookup from Section 6e:", length(already_in_lookup), "\n")
cat("Never queried at all - checking now:", length(names_to_check), "\n\n")

cat("Pass 1: exact match,", length(names_to_check), "names\n")
pass1_6h <- query_worms_exact(names_to_check)
names_not_found_1_6h <- pass1_6h %>%
  filter(is.na(scientificname) | status == "not found") %>%
  pull(queried_name) %>% unique()
names_not_found_1_6h <- union(names_not_found_1_6h, setdiff(names_to_check, pass1_6h$queried_name))
cat("  Resolved:", length(names_to_check) - length(names_not_found_1_6h),
    "| Remaining:", length(names_not_found_1_6h), "\n\n")

cat("Pass 2: taxamatch fallback,", length(names_not_found_1_6h), "names\n")
pass2_6h <- query_worms_taxamatch(names_not_found_1_6h, batch_size = 10)
pass2_6h_resolved <- pass2_6h %>% filter(!is.na(status) & status != "not found")
names_not_found_2_6h <- setdiff(names_not_found_1_6h, pass2_6h_resolved$queried_name)
cat("  Resolved:", length(unique(pass2_6h_resolved$queried_name)),
    "| Remaining:", length(names_not_found_2_6h), "\n\n")

if (length(names_not_found_2_6h) > 0) {
  cat("Names with NO automated WoRMS match (manually verified below):\n")
  print(names_not_found_2_6h)
}

# Manual corrections for names with no automated WoRMS match (confirmed by
# direct manual lookup against marinespecies.org)
manual_corrections_6h <- tribble(
  ~queried_name,        ~resolved_name,                       ~worms_rank, ~worms_status,
  "Actinopyga lubrica", "Holothuria (Selenkothuria) lubrica", "Species",   "manual_correction",
  "Hemiaster",          "Hemiaster",                          "Genus",     "manual_correction",
  "Iconometra dormani", "Iconometra",                         "Genus",     "manual_correction"  # WoRMS unresolvable
)

exact_lookup_6h <- pass1_6h %>%
  filter(!is.na(scientificname) & status != "not found") %>%
  mutate(status_rank = coalesce(status_priority[status], 10)) %>%
  group_by(queried_name) %>%
  slice_min(status_rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(resolved_name = if_else(status == "accepted", scientificname,
                                 coalesce(valid_name, scientificname))) %>%
  select(queried_name, resolved_name, worms_status = status, worms_rank = rank)

taxamatch_lookup_6h <- pass2_6h_resolved %>%
  mutate(status_rank = coalesce(status_priority[status], 10)) %>%
  group_by(queried_name) %>%
  slice_min(status_rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(resolved_name = scientificname) %>%
  select(queried_name, resolved_name, worms_status = status, worms_rank = rank)

new_resolutions <- bind_rows(exact_lookup_6h, taxamatch_lookup_6h, manual_corrections_6h) %>%
  distinct(queried_name, .keep_all = TRUE)

cat("\nTotal newly resolved via Section 6h:", nrow(new_resolutions), "of", length(names_to_check), "\n\n")

# Names confirmed accepted-as-is in WoRMS despite "mixed status" exclusion
mixed_status_fix <- tribble(
  ~queried_name,                          ~resolved_name,                          ~worms_status, ~worms_rank,
  "Antedon loveni",                       "Antedon loveni",                        "accepted",    "Species",
  "Polyplectana nigra",                   "Polyplectana nigra",                    "accepted",    "Species",
  "Temnopleurus alexandri",               "Temnopleurus alexandri",                "accepted",    "Species",
  "Ophiura robusta",                      "Ophiura robusta",                       "accepted",    "Species",
  "Peltaster placenta",                   "Peltaster placenta",                    "accepted",    "Species",
  "Ophionereis (Ophiotriton) dubia",      "Ophionereis (Ophiotriton) dubia",       "accepted",    "Species",
  "Ophionereis (Ophiotriton) amoyensis",  "Ophionereis (Ophiotriton) amoyensis",   "accepted",    "Species",
  "Ophiothrix (Acanthophiothrix)",        "Ophiothrix (Acanthophiothrix)",         "accepted",    "Subgenus"
)

master_status_lookup <- bind_rows(
  final_name_lookup %>% select(queried_name, worms_status, worms_rank),
  manual_corrections %>% select(queried_name, worms_status, worms_rank),
  new_resolutions %>% select(queried_name, worms_status, worms_rank),
  mixed_status_fix %>% select(queried_name, worms_status, worms_rank)
) %>%
  distinct(queried_name, .keep_all = TRUE)

cat("Master lookup covers", nrow(master_status_lookup), "distinct names\n")

echino_wide <- echino_wide %>%
  select(-any_of(c("worms_status_master", "worms_rank_master"))) %>%
  left_join(
    master_status_lookup %>% rename(accepted_name = queried_name,
                                    worms_status_master = worms_status,
                                    worms_rank_master = worms_rank),
    by = "accepted_name"
  )

cat("Records with a known WoRMS status (master):",
    sum(!is.na(echino_wide$worms_status_master)), "of", nrow(echino_wide),
    sprintf("(%.1f%%)\n\n", 100 * mean(!is.na(echino_wide$worms_status_master))))

# -----------------------------------------------------------------------------
# Manually verified current WoRMS-valid names for every record whose status
# isn't "accepted". POST-AUDIT REVISION: two entries corrected following
# independent cell-by-cell audit and direct AphiaID verification.
# -----------------------------------------------------------------------------

verified_valid_names <- tribble(
  ~accepted_name,                              ~verified_valid_name,
  # --- WoRMS unresolvable: manually verified against marinespecies.org ---
  "Iconometra dormani",                        "Iconometra",          # No WoRMS record; genus retained
  # Fromia japonica already resolved to Fromia monilis via manual_corrections
  # in Section 6d — entry below is a safety net only, will not fire in practice
  "Fromia japonica",                           "Fromia monilis",     # Unaccepted synonym
  # --- Unaccepted synonyms: WoRMS-resolvable but accepted_name_final
  # corrected here to ensure the current valid combination is used ---
  "Ophiocoma variegata",                       "Breviturma dentata", # AphiaID 1214745
  "Fibularia volva",                           "Fibulariella volva", # Incorrect genus combination
  # --- Subgenus corrections ---
  "Ophionereis semoni",                        "Ophionereis (Ophiotriton) semoni",
  "Ophionereis tigris",                        "Ophionereis (Ophiotriton) tigris",
  "Amphiophiura spatulifera",                  "Ophiopyrgus spatulifera",
  "Ophionereis hexactis",                      "Ophionereis (Ophiotriton) hexactis",
  "Ophionereis intermedia",                    "Ophionereis (Ophiotriton) intermedia",
  "Gymnocrinus richeri",                       "Neogymnocrinus richeri",
  "Oneirophanta mutabilis",                    "Oneirophanta mutabilis mutabilis",
  "Actinopyga lubrica",                        "Holothuria (Selenkothuria) lubrica",
  "Amphioplus lucidus",                        "Amphioplus (Amphioplus) lucidus",
  "Asterina burtoni",                          "Aquilonastra burtoni",
  "Benthodytes papillifera",                   "Benthodytes typica",
  "Colochirus tuberculosus",                   "Colochirus",
  "Echinoneus abnormalis",                     "Koehleraster abnormalis",
  "Hemiaster",                                 "Hemiaster",
  "Holothuria (semperothuria) flavomaculata",  "Holothuria flavomaculata",
  "Ophiothrix propinqua",                      "Macrophiothrix propinqua",
  "Placophiothrix trilineata",                 "Ophiothrix (Ophiothrix) trilineata",
  "Synapta glabra",                            "Opheodesoma glabra",
  "Synapta grisea var. alba",                  "Opheodesoma grisea",
  "Synapta indivisa",                          "Synaptula indivisa",
  "Synapta similis",                           "Protankyra similis",
  "Amphioplus (Amphichilus) ochroleuca",       "Amphioplus (Amphichilus) ochroleucus",
  "Amphioplus (Amphioplus) intermedius",       "Amphioplus (Amphichilus) intermedius",
  "Amphiodia (amphispina) microplax",          "Amphiodia (Amphispina) microplax",
  "Holothuria (theelothuria) squamifera",      "Holothuria (Theelothuria) squamifera",
  "Ophiothrix (ophiothrix) panchyendyta",      "Ophiothrix (Ophiothrix) panchyendyta",
  "Ophiura (dictenophiura) squamosa",          "Ophiura (Dictenophiura) squamosa",
  "Lymanella",                                 "Amphioplus (Lymanella)",
  "Cheiraster gazellae",                       "Cheiraster (Cheiraster) gazellae"
)

echino_wide <- echino_wide %>%
  select(-any_of(c("verified_valid_name", "accepted_name_final"))) %>%
  left_join(verified_valid_names, by = "accepted_name")

echino_wide <- echino_wide %>%
  mutate(
    accepted_name_final = coalesce(verified_valid_name, accepted_name)
  )

# For AM records with only class-level identification (aa-* family codes),
# assign the class name as accepted_name and accepted_name_final
echino_wide <- echino_wide %>%
  mutate(
    accepted_name = if_else(
      is.na(accepted_name) & str_detect(coalesce(family__AM_direct, ""), "^aa-"),
      str_extract(coalesce(family__AM_direct, ""), "(?<=aa-).*"),
      accepted_name
    ),
    accepted_name_final = if_else(
      is.na(accepted_name_final) & str_detect(coalesce(family__AM_direct, ""), "^aa-"),
      str_extract(coalesce(family__AM_direct, ""), "(?<=aa-).*"),
      accepted_name_final
    ),
    worms_status_master = if_else(
      is.na(worms_status_master) & str_detect(coalesce(family__AM_direct, ""), "^aa-"),
      "accepted",
      worms_status_master
    ),
    # worms_rank_master must be patched here too, not just worms_status_master -
    # otherwise the correct "Class" assignment below gets silently overwritten
    # by the word-count-based recompute later in Section 6i (a bare class name
    # like "Asteroidea" is a single word, and with no worms_rank_master to
    # correct it against, the word-count rule wrongly calls it "Genus").
    worms_rank_master = if_else(
      is.na(worms_rank_master) & str_detect(coalesce(family__AM_direct, ""), "^aa-"),
      "Class",
      worms_rank_master
    ),
    taxonomic_resolution_level = if_else(
      str_detect(coalesce(family__AM_direct, ""), "^aa-"),
      "Class",
      taxonomic_resolution_level
    )
  )

# accepted_name_final still NA for any records where accepted_name
# was just set above but accepted_name_final was not caught (single-mutate
# evaluation order issue). Propagate accepted_name -> accepted_name_final.
echino_wide <- echino_wide %>%
  mutate(
    accepted_name_final = if_else(
      is.na(accepted_name_final) & !is.na(accepted_name),
      accepted_name,
      accepted_name_final
    ),
    # Also fix taxonomic_resolution_level for the aa- records now that
    # accepted_name is no longer NA
    taxonomic_resolution_level = if_else(
      str_detect(coalesce(family__AM_direct, ""), "^aa-") &
        taxonomic_resolution_level != "Class",
      "Class",
      taxonomic_resolution_level
    )
  )

cat("Records updated with verified valid names:\n")
echino_wide %>%
  filter(!is.na(verified_valid_name)) %>%
  count(accepted_name, verified_valid_name, sort = TRUE) %>%
  print(n = Inf)

cat("\nFinal Section 6h record count:", nrow(echino_wide), "\n")
cat("Final WoRMS validation coverage:",
    sprintf("%.1f%%\n", 100 * mean(!is.na(echino_wide$worms_status_master))))

write_csv(new_resolutions, "section6h_newly_resolved_names.csv")

all_non_accepted_summary <- echino_wide %>%
  filter(!is.na(worms_status_master) & worms_status_master != "accepted") %>%
  count(accepted_name, accepted_name_final, worms_status_master, worms_rank_master, sort = TRUE)
write_csv(all_non_accepted_summary, "all_non_accepted_names_summary.csv")

# -----------------------------------------------------------------------------
# SECTION 6i: POST-AUDIT FIX - extinct genus flag
# Hemiaster is intentionally RETAINED (basisOfRecord = PreservedSpecimen,
# not Fossil; USNM Invertebrate Zoology; sourced via OBIS), but is flagged
# explicitly here for analytical transparency.
# -----------------------------------------------------------------------------

cat("\n== Section 6i: POST-AUDIT FIX - extinct genus flag ==\n\n")

echino_wide <- echino_wide %>%
  mutate(is_extinct_genus = coalesce(accepted_name_final == "Hemiaster", FALSE))

cat("Records flagged is_extinct_genus = TRUE:", sum(echino_wide$is_extinct_genus), "\n")

# Hemiaster's aphiaID is hardcoded here (123423, https://www.marinespecies.org/
# aphia.php?p=taxdetails&id=123423) rather than left to the aphiaID-based or
# name-based lookups in Section 6k. Extinct/fossil-adjacent genera like this
# one often carry unusual status values in WoRMS that make automated name
# matching (wm_records_names("Hemiaster")) unreliable, so this guarantees the
# Section 6k aphiaID join picks up the correct family/order for this record
# directly, without depending on that lookup behaving well for this one case.
# Note: no dagger/extinction marker is added to accepted_name_final itself -
# that field feeds both the WoRMS name lookup and the eventual Darwin Core
# scientificName, and a "†" would break both; is_extinct_genus (above) is
# the structured place this information belongs.
echino_wide <- echino_wide %>%
  mutate(aphiaID = if_else(is_extinct_genus, "123423", as.character(aphiaID)))

# -----------------------------------------------------------------------------
# POST-AUDIT FIX: recompute taxonomic_resolution_level from accepted_name_final
# -----------------------------------------------------------------------------
# The fallback_rank logic in Section 6e reads genus__*/specificEpithet__* source
# columns, which do not always reflect the final resolved name in accepted_name_final.
# This produces two systematic errors:
#   - 3,721 genuine binomials flagged as "Genus" (underreports species-level ID)
#   - 177 single-word names flagged as "Species" (overreports species-level ID)
# Net effect: species-level count underreported by ~3,544 records (~12%).
# Fix: reassign resolution level directly from the word structure of
# accepted_name_final, which is authoritative at this point in the pipeline.
# -----------------------------------------------------------------------------

cat("== Recomputing taxonomic_resolution_level from accepted_name_final ==\n\n")

cat("Before fix:\n")
print(table(echino_wide$taxonomic_resolution_level))

echino_wide <- echino_wide %>%
  mutate(
    taxonomic_resolution_level = case_when(
      
      # Subspecies: 3 words, no parentheses e.g. "Oneirophanta mutabilis mutabilis"
      str_count(accepted_name_final, "\\S+") == 3 &
        !str_detect(coalesce(accepted_name_final, ""), "\\(") ~ "Subspecies",
      
      # Species: genuine binomial - 2+ words, not a subgenus-only name
      # Subgenus-only looks like: "Genus (Subgenus)" - exactly 2 tokens where
      # second token is entirely wrapped in parentheses
      str_count(accepted_name_final, "\\S+") >= 2 &
        !str_detect(coalesce(accepted_name_final, ""), "^\\S+ \\(\\S+\\)$") ~ "Species",
      
      # Subgenus-only: "Genus (Subgenus)" pattern
      str_detect(coalesce(accepted_name_final, ""), "^\\S+ \\(\\S+\\)$") ~ "Subgenus",
      
      # Single word = Genus (or higher, but higher ranks already have
      # accepted_name_final set to a class/order/family name - handled below)
      str_count(accepted_name_final, "\\S+") == 1 ~ "Genus",
      
      is.na(accepted_name_final) ~ taxonomic_resolution_level,
      
      TRUE ~ taxonomic_resolution_level
    )
  )

# Higher ranks: where accepted_name_final is a known class/order/family/etc.
# name, the single-word rule above assigns "Genus" which is wrong. Use
# worms_rank_master to correct these where available - this REPLACES an
# earlier version of this fix that only handled 4 rank groups (Class/
# Subclass/Superclass, Order/Suborder/Superorder, Family/Subfamily/
# Superfamily, Phylum) and silently left anything else (Subterclass,
# Infraclass, Infraorder, Parvorder, Cohort, etc.) as the wrong word-count
# guess. Confirmed cases this missed: "Articulata" (Subclass, was showing
# as "Class" - technically inside the old mapping but collapsing the rank
# distinction), "Echinacea" (Subterclass - not in the old list at all, so
# stayed wrongly bucketed as "Genus"), and an Infraorder-level name
# ("Echinidea") also wrongly bucketed as "Genus".
#
# This version maps every rank WoRMS is realistically going to return for
# echinoderm records to one of this pipeline's 8 target buckets, based on
# where that rank actually sits in the taxonomic hierarchy - not a hand
# picked shortlist that has to be extended every time a new rank shows up.
worms_rank_to_bucket <- c(
  # Phylum and above
  "Kingdom" = "Phylum", "Subkingdom" = "Phylum", "Infrakingdom" = "Phylum",
  "Superphylum" = "Phylum", "Phylum" = "Phylum", "Subphylum" = "Phylum",
  "Infraphylum" = "Phylum",
  # Class-level (includes the less common intermediate ranks WoRMS uses for
  # echinoderms specifically, e.g. Subterclass for Echinacea/Irregularia)
  "Superclass" = "Class", "Class" = "Class", "Subclass" = "Class",
  "Infraclass" = "Class", "Subterclass" = "Class", "Cohort" = "Class",
  "Supercohort" = "Class", "Subcohort" = "Class",
  # Order-level
  "Superorder" = "Order", "Order" = "Order", "Suborder" = "Order",
  "Infraorder" = "Order", "Parvorder" = "Order",
  # Family-level
  "Superfamily" = "Family", "Family" = "Family", "Subfamily" = "Family",
  "Tribe" = "Family", "Subtribe" = "Family",
  # Genus-level
  "Genus" = "Genus", "Subgenus" = "Subgenus",
  # Species-level and finer
  "Species" = "Species", "Subspecies" = "Subspecies",
  "Variety" = "Subspecies", "Forma" = "Subspecies"
)

# Diagnostic: flag any worms_rank_master value NOT covered by the mapping
# above, so a future WoRMS rank this pipeline hasn't seen before surfaces
# immediately as a visible warning instead of silently falling through to
# the word-count guess the way Subterclass/Infraorder did here.
unmapped_ranks <- echino_wide %>%
  filter(!is.na(worms_rank_master), !worms_rank_master %in% names(worms_rank_to_bucket)) %>%
  count(worms_rank_master, sort = TRUE)

if (nrow(unmapped_ranks) > 0) {
  cat("\nWARNING:", nrow(unmapped_ranks), "worms_rank_master value(s) not covered by",
      "worms_rank_to_bucket - these records will fall back to the word-count guess,",
      "which is exactly the bug this mapping exists to avoid. Add them above.\n")
  print(unmapped_ranks)
} else {
  cat("\nAll present worms_rank_master values are covered by worms_rank_to_bucket.\n")
}

echino_wide <- echino_wide %>%
  mutate(
    taxonomic_resolution_level = if_else(
      !is.na(worms_rank_master) & worms_rank_master %in% names(worms_rank_to_bucket),
      worms_rank_to_bucket[worms_rank_master],
      taxonomic_resolution_level
    )
  )

cat("\nAfter fix:\n")
print(table(echino_wide$taxonomic_resolution_level))

cat("\n=== Section 6 (including 6h, 6i) complete ===\n")
cat("Final dataset:", nrow(echino_wide), "records, 100% WoRMS-validated via accepted_name_final.\n")

# =============================================================================
# SECTION 6j: ECHINODERM CLASS RESOLUTION
# =============================================================================
# `class` is a conflict field (per-source class__* columns), so there is no
# single consensus `class` column produced upstream. This section builds one,
# resolving each record's class in two stages:
#   1. Direct resolution from the per-source class__* columns, in source-
#      priority order, case-insensitive (fixes 226 records that were
#      previously double-counted as separate lowercase/uppercase categories:
#      e.g. "asteroidea" vs "Asteroidea").
#   2. For the records still unresolved after stage 1, a taxon (genus/family/
#      order) -> class lookup is applied, but ONLY where accepted_name_final
#      gives a finer-than-Phylum identification (i.e. there IS a name to
#      infer class from). Records identified only to Phylum are correctly
#      left as NA - there is no class to know.
#   `echino_class_source` records which of the two paths resolved each
#   record, so the direct-vs-inferred distinction stays visible in the
#   final dataset rather than being silently merged.
# =============================================================================

cat("\n== Section 6j: Echinoderm class resolution ==\n\n")

# -----------------------------------------------------------------------------
# 6j-1: Direct resolution from source class__* columns (source-priority order,
# case-insensitive)
# -----------------------------------------------------------------------------

# Uses the single shared SOURCE_PRIORITY defined at the top of the script
# (now includes AM_direct, which this list was previously missing).
echino_class_cols <- paste0("class__", SOURCE_PRIORITY)
echino_class_cols <- intersect(echino_class_cols, names(echino_wide))

valid_classes <- c("asteroidea", "ophiuroidea", "echinoidea",
                   "holothuroidea", "crinoidea")

echino_wide <- echino_wide %>%
  mutate(
    echino_class = pmap_chr(select(., all_of(echino_class_cols)), function(...) {
      vals <- c(...)
      vals <- vals[!is.na(vals) & str_trim(vals) != ""]
      vals <- vals[str_to_lower(str_trim(vals)) %in% valid_classes]
      if (length(vals) == 0) return(NA_character_)
      str_to_title(str_trim(vals[1]))
    })
  )

cat("── Direct resolution from source class__* columns ──\n")
echino_wide %>% count(echino_class, sort = TRUE) %>% print()

# -----------------------------------------------------------------------------
# 6j-2 STEP 1: identify what the echino_class = NA records actually are
# -----------------------------------------------------------------------------

cat("\n── Step 1: identifying echino_class = NA records ──\n\n")

na_records <- echino_wide %>% filter(is.na(echino_class))
cat("Total records with echino_class = NA:", nrow(na_records), "\n\n")

na_records %>% count(taxonomic_resolution_level, sort = TRUE) %>% print()

phylum_only   <- na_records %>% filter(taxonomic_resolution_level == "Phylum")
finer_records <- na_records %>% filter(taxonomic_resolution_level != "Phylum")

cat("\nPhylum-only (nothing to recover):", nrow(phylum_only), "\n")
cat("Finer resolution (class inferable from name):", nrow(finer_records), "\n\n")

cat("── Sources contributing the finer, recoverable records ──\n")
finer_records %>% count(primary_source, sort = TRUE) %>% print()

taxa_needing_lookup <- finer_records %>%
  filter(!taxonomic_resolution_level %in% c("Class", "Phylum")) %>%
  mutate(lookup_key = word(accepted_name_final, 1)) %>%
  count(lookup_key, taxonomic_resolution_level, sort = TRUE)

cat("\nDistinct taxon keys requiring a lookup entry:",
    n_distinct(taxa_needing_lookup$lookup_key), "\n\n")

# -----------------------------------------------------------------------------
# 6j-2 STEP 2: taxon -> class lookup, covering the keys found in Step 1
# (standard echinoderm taxonomy; family/order assignments per WoRMS)
# -----------------------------------------------------------------------------

taxon_class_lookup <- tribble(
  ~taxon,             ~taxon_class,
  
  # --- Class ---
  "Ophiuroidea",    "Ophiuroidea",
  "Asteroidea",     "Asteroidea",
  "Echinoidea",     "Echinoidea",
  "Holothuroidea",  "Holothuroidea",
  "Crinoidea",      "Crinoidea",
  
  # --- Orders ---
  "Ophiurida",         "Ophiuroidea",
  "Comatulida",        "Crinoidea",
  
  # --- Families ---
  "Benthopectinidae",  "Asteroidea",
  "Echinasteridae",    "Asteroidea",
  "Goniasteridae",     "Asteroidea",
  "Astropectinidae",   "Asteroidea",
  "Ophiacanthidae",    "Ophiuroidea",
  "Amphiuridae",       "Ophiuroidea",
  "Ophiactidae",       "Ophiuroidea",
  "Cidaridae",         "Echinoidea",
  "Diadematidae",      "Echinoidea",
  
  # --- Genera ---
  "Ophiomusium",       "Ophiuroidea",
  "Ophiomusa",         "Ophiuroidea",
  "Hymenaster",        "Asteroidea",
  "Asteroschema",      "Ophiuroidea",
  "Ophiotreta",        "Ophiuroidea",
  "Araeosoma",         "Echinoidea",
  "Echinacea",         "Echinoidea",
  "Ophiactis",         "Ophiuroidea",
  "Metacrinus",        "Crinoidea",
  "Astropecten",       "Asteroidea",
  "Perissogonaster",   "Asteroidea",
  "Paracaudina",       "Holothuroidea",
  "Echinosigra",       "Echinoidea",
  "Astrobrachion",     "Ophiuroidea",
  "Cryptasterina",     "Asteroidea",
  "Phormosoma",        "Echinoidea",
  "Cheiraster",        "Asteroidea",
  "Neogymnocrinus",    "Crinoidea",
  "Echinaster",        "Asteroidea",
  "Chaetodiadema",     "Echinoidea",
  "Luidia",            "Asteroidea",
  "Nymphaster",        "Asteroidea",
  "Actinopyga",        "Holothuroidea",
  "Bohadschia",        "Holothuroidea",
  "Metrodira",         "Echinoidea",
  "Holothuria",        "Holothuroidea",
  "Linckia",          "Asteroidea",
  "Fromia",           "Asteroidea",
  "Gomophia",         "Asteroidea",
  "Ophidiaster",      "Asteroidea",
  "Aquilonastra",     "Asteroidea",
  "Stellaster",       "Asteroidea",
  "Goniodiscaster",   "Asteroidea",
  "Pentaceraster",    "Asteroidea",
  "Anthenea",         "Asteroidea",
  "Peronella",        "Echinoidea",
  "Salmaciella",      "Echinoidea",
  "Temnopleurus",     "Echinoidea",
  "Ophiothrix",       "Ophiuroidea",
  "Macrophiothrix",   "Ophiuroidea",
  "Iconometra",       "Crinoidea",
  "Himerometra",      "Crinoidea",
  "Euantedon",        "Crinoidea",
  "Comanthus",        "Crinoidea"
)

cat("Taxon -> class lookup entries:", nrow(taxon_class_lookup), "\n\n")

uncovered_keys <- setdiff(taxa_needing_lookup$lookup_key, taxon_class_lookup$taxon)
if (length(uncovered_keys) == 0) {
  cat("Lookup table covers all", n_distinct(taxa_needing_lookup$lookup_key),
      "taxon keys found in Step 1.\n\n")
} else {
  cat("WARNING: lookup table is missing", length(uncovered_keys), "taxon key(s):\n")
  print(uncovered_keys)
  cat("These records will remain unresolved until entries are added above.\n\n")
}

# -----------------------------------------------------------------------------
# 6j-2 STEP 3: apply the fix
#   - taxonomic_resolution_level == "Class"  -> accepted_name_final IS the
#     class name already, use it directly (no lookup needed)
#   - otherwise (Genus/Family/Order/Species) -> match the genus/family/order
#     token (first word of accepted_name_final) against the lookup
# -----------------------------------------------------------------------------

echino_wide <- echino_wide %>%
  mutate(
    .lookup_key = if_else(
      !is.na(worms_rank_master) & worms_rank_master == "Class",
      accepted_name_final,
      word(accepted_name_final, 1)
    )
  ) %>%
  left_join(taxon_class_lookup, by = c(".lookup_key" = "taxon")) %>%
  mutate(
    # Was: taxonomic_resolution_level == "Class" - correct while that bucket
    # only ever meant the exact rank "Class", but now that worms_rank_master
    # collapses Subclass/Superclass/Infraclass/Subterclass/Cohort into the
    # same "Class" bucket too (see the taxonomic_resolution_level recompute
    # above), accepted_name_final for those records is NOT a valid class name
    # - e.g. "Echinacea" (a Subterclass) is not one of the five echinoderm
    # classes, its parent class is "Echinoidea". Using the EXACT rank here
    # instead of the bucket avoids writing a Subterclass/Infraclass name into
    # echino_class in the first place, rather than relying on Section 6k's
    # cross-check to catch it after the fact (which it did, correctly, but
    # labelled it as a source field error - it wasn't one, this was).
    inferred_class = if_else(
      !is.na(worms_rank_master) & worms_rank_master == "Class",
      accepted_name_final,
      taxon_class
    ),
    echino_class_source = case_when(
      !is.na(echino_class) ~ "original_source_field",
      is.na(echino_class) & !is.na(inferred_class) &
        taxonomic_resolution_level != "Phylum" ~ "inferred_from_taxon_name",
      TRUE ~ NA_character_
    ),
    echino_class = if_else(
      is.na(echino_class) & taxonomic_resolution_level != "Phylum" & !is.na(inferred_class),
      inferred_class,
      echino_class
    )
  ) %>%
  select(-.lookup_key, -taxon_class, -inferred_class)

cat("── After class recovery fix ──\n")
echino_wide %>% count(echino_class, sort = TRUE) %>% print()
cat("\n── Recovery breakdown ──\n")
echino_wide %>% count(echino_class_source, sort = TRUE) %>% print()

# Handle AM's internal class-level coding (aa-Asteroidea etc.)
am_class_pattern <- "^aa-(.+)$"
echino_wide <- echino_wide %>%
  mutate(
    echino_class = if_else(
      is.na(echino_class) & str_detect(coalesce(family__AM_direct, ""), am_class_pattern),
      str_extract(coalesce(family__AM_direct, ""), "(?<=aa-).*"),
      echino_class
    ),
    echino_class_source = if_else(
      is.na(echino_class_source) & str_detect(coalesce(family__AM_direct, ""), am_class_pattern),
      "inferred_from_AM_class_code",
      echino_class_source
    )
  )

# -----------------------------------------------------------------------------
# 6j-2 STEP 4: verification - confirm every remaining NA is genuinely
# Phylum-only
# -----------------------------------------------------------------------------

n_na_total     <- sum(is.na(echino_wide$echino_class))
n_na_phylum    <- sum(is.na(echino_wide$echino_class) & echino_wide$taxonomic_resolution_level == "Phylum")
n_na_nonphylum <- sum(is.na(echino_wide$echino_class) & echino_wide$taxonomic_resolution_level != "Phylum")

cat("\n── Verification ──\n")
cat("Total remaining echino_class = NA:  ", n_na_total, "\n")
cat("  - Phylum-only (expected, OK):     ", n_na_phylum, "\n")
cat("  - Finer resolution (should be 0): ", n_na_nonphylum, "\n")

if (n_na_nonphylum == 0) {
  cat("VERIFIED: all remaining echino_class = NA records are genuine Phylum-only identifications.\n")
} else {
  cat("WARNING:", n_na_nonphylum, "finer-resolution records still have no class assigned -",
      "add missing taxon_class_lookup entries above.\n")
}

cat("\n=== Section 6j complete ===\n")

# =============================================================================
# SECTION 6k: HIGHER CLASSIFICATION FROM WoRMS (family, order)
# =============================================================================
# family and order were never resolved to a single consensus value (only
# per-source family__SOURCE/order__SOURCE columns exist, unreconciled - same
# situation basisOfRecord was in before Section 4's addition). Rather than
# coalescing those unreliable per-source columns, this pulls family and order
# straight from WoRMS via the aphiaID already resolved earlier in Section 6,
# batched by DISTINCT aphiaID (not by record) since many records share the
# same accepted name - keeps this to ~1,000 API calls instead of ~43,000.
#
# Records with no aphiaID (no confident WoRMS match, e.g. the Phylum-only
# records from 6j-2) get best_family/best_order = NA - expected, not a bug.
# =============================================================================

cat("\n== Section 6k: WoRMS higher classification (family, order) ==\n\n")

aphia_ids <- echino_wide$aphiaID %>% na.omit() %>% unique() %>% sort()
cat("Distinct aphiaIDs to resolve:", length(aphia_ids), "\n")

# wm_record_() (the batched/vectorized endpoint) was confirmed on a real
# run to return a row per AphiaID but with every classification field NA -
# WoRMS' batch endpoint doesn't return full record detail the way the
# single-record endpoint does. So this queries wm_record() once per
# distinct AphiaID instead. Slower (1,031 individual calls rather than
# ~21 batches) but it's the only endpoint that actually returns
# kingdom/phylum/class/order/family/genus. A short pause between calls
# avoids hammering the API; progress is printed every 100 so a run that's
# taking a while is visibly still working rather than looking hung.
query_worms_classification <- function(ids) {
  if (length(ids) == 0) return(tibble())
  n <- length(ids)
  results <- vector("list", n)
  for (i in seq_along(ids)) {
    r <- tryCatch(wm_record(id = ids[i]), error = function(e) NULL)
    results[[i]] <- if (is.null(r) || nrow(r) == 0) {
      tibble(AphiaID = ids[i], kingdom = NA_character_, phylum = NA_character_,
             class = NA_character_, order = NA_character_, family = NA_character_,
             genus = NA_character_)
    } else {
      r
    }
    if (i %% 100 == 0) cat("  ...", i, "of", n, "AphiaIDs queried\n")
    Sys.sleep(0.05)
  }
  bind_rows(results)
}

worms_classification <- query_worms_classification(aphia_ids) %>%
  select(aphiaID = AphiaID, worms_kingdom = kingdom, worms_phylum = phylum,
         worms_class = class, worms_order = order, worms_family = family,
         worms_genus = genus) %>%
  distinct(aphiaID, .keep_all = TRUE)

cat("Classification records returned:", nrow(worms_classification), "of",
    length(aphia_ids), "requested\n")

# The join below matched 0 rows on a real run despite 1031/1031 classification
# records coming back correctly - the cause was a silent type mismatch between
# echino_wide$aphiaID (double) and the AphiaID column WoRMS returns (worrms
# does not guarantee this comes back as double). Forcing both sides to the
# same integer type before joining, rather than trusting dplyr's implicit
# coercion, so this can't silently produce an all-NA result again.
cat("aphiaID type - echino_wide:", class(echino_wide$aphiaID),
    "| worms_classification:", class(worms_classification$aphiaID), "\n")

echino_wide <- echino_wide %>%
  mutate(aphiaID = as.integer(aphiaID))
worms_classification <- worms_classification %>%
  mutate(aphiaID = as.integer(aphiaID))

# If Section 6k has already been run once in this R session (e.g. while
# debugging), echino_wide will already carry best_family/best_order/worms_*
# columns from that earlier attempt. Re-joining without dropping them first
# produces two columns with the same name, and rename() below fails with
# "Names must be unique." Dropping any pre-existing copies here makes this
# section safe to re-run any number of times in the same session.
echino_wide <- echino_wide %>%
  select(-any_of(c("best_family", "best_order", "worms_kingdom", "worms_phylum",
                   "worms_class", "worms_genus")))

n_potential_matches <- sum(echino_wide$aphiaID %in% worms_classification$aphiaID, na.rm = TRUE)
cat("Records whose aphiaID has a match in worms_classification (pre-join check):",
    n_potential_matches, "of", sum(!is.na(echino_wide$aphiaID)), "\n")

echino_wide <- echino_wide %>%
  left_join(worms_classification, by = "aphiaID") %>%
  rename(best_family = worms_family, best_order = worms_order)

n_with_family <- sum(!is.na(echino_wide$best_family))
n_with_order  <- sum(!is.na(echino_wide$best_order))
cat(sprintf("Records with best_family resolved: %d of %d (%.1f%%)\n",
            n_with_family, nrow(echino_wide), 100 * n_with_family / nrow(echino_wide)))
cat(sprintf("Records with best_order resolved:  %d of %d (%.1f%%)\n",
            n_with_order, nrow(echino_wide), 100 * n_with_order / nrow(echino_wide)))

# -----------------------------------------------------------------------------
# SECTION 6k addendum: name-based WoRMS fallback for records with no aphiaID
# -----------------------------------------------------------------------------
# aphiaID is a raw pass-through field only present for sources that supply it
# natively (OBIS, ALA) - it covers only 23,045 of 43,470 records, regardless
# of whether the record's NAME was actually resolved against WoRMS. 43,352 of
# 43,470 records (99.7%) have a validated accepted_name_final. Rather than
# refactoring Sections 6c/6e/6h to capture and thread AphiaID through the
# name-resolution passes (a bigger job), this queries WoRMS directly BY NAME -
# using the same wm_records_names() function already used for name resolution
# elsewhere in this script - for every record still missing best_family after
# the aphiaID-based lookup above. Batched by DISTINCT accepted_name_final
# (not by record), same principle as the aphiaID batching above.
#
# A name can occasionally return more than one WoRMS record (rare homonyms);
# the "accepted" status match is preferred, using the same status_priority
# ranking already defined and used for name resolution in Section 6c.
#
# Existing aphiaID-based best_family/best_order values are never overwritten
# here - coalesce() only fills genuine gaps. aphiaID itself is also backfilled
# where this lookup finds one and it was previously missing, since that
# improves any other downstream use of the aphiaID column too.
# -----------------------------------------------------------------------------

cat("\n== Section 6k addendum: name-based WoRMS classification fallback ==\n\n")

names_needing_classification <- echino_wide %>%
  filter(is.na(best_family) & !is.na(accepted_name_final)) %>%
  distinct(accepted_name_final) %>%
  pull(accepted_name_final) %>%
  sort()

cat("Distinct names to query (records still missing best_family after aphiaID lookup):",
    length(names_needing_classification), "\n")

query_worms_classification_by_name <- function(names_vec) {
  if (length(names_vec) == 0) return(tibble())
  empty_result <- function(nm) {
    tibble(queried_name = nm, AphiaID = NA_integer_, kingdom = NA_character_,
           phylum = NA_character_, class = NA_character_, order = NA_character_,
           family = NA_character_, genus = NA_character_, status = NA_character_)
  }
  batches <- split(names_vec, ceiling(seq_along(names_vec) / 50))
  map_dfr(batches, function(b) {
    res <- tryCatch(wm_records_names(name = b, fuzzy = FALSE), error = function(e) NULL)
    if (is.null(res)) {
      return(map_dfr(b, function(nm) {
        r <- tryCatch(wm_records_names(name = nm, fuzzy = FALSE)[[1]], error = function(e2) NULL)
        if (is.null(r) || nrow(r) == 0) empty_result(nm) else r %>% mutate(queried_name = nm)
      }))
    }
    map2_dfr(b, res, function(nm, r) {
      if (is.null(r) || nrow(r) == 0) empty_result(nm) else r %>% mutate(queried_name = nm)
    })
  })
}

name_classification_raw <- query_worms_classification_by_name(names_needing_classification)
cat("Classification records returned:", nrow(name_classification_raw), "of",
    length(names_needing_classification), "names queried\n")

name_classification <- name_classification_raw %>%
  mutate(status_rank = coalesce(status_priority[status], 10)) %>%
  group_by(queried_name) %>%
  slice_min(status_rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(accepted_name_final = queried_name, name_aphiaID = AphiaID,
         name_kingdom = kingdom, name_phylum = phylum, name_class = class,
         name_order = order, name_family = family, name_genus = genus)

echino_wide <- echino_wide %>%
  select(-any_of(c("name_aphiaID", "name_kingdom", "name_phylum", "name_class",
                   "name_order", "name_family", "name_genus"))) %>%
  left_join(name_classification, by = "accepted_name_final") %>%
  mutate(
    aphiaID       = coalesce(aphiaID, name_aphiaID),
    worms_kingdom = coalesce(worms_kingdom, name_kingdom),
    worms_phylum  = coalesce(worms_phylum, name_phylum),
    worms_class   = coalesce(worms_class, name_class),
    best_family   = coalesce(best_family, name_family),
    best_order    = coalesce(best_order, name_order)
  ) %>%
  select(-name_aphiaID, -name_kingdom, -name_phylum, -name_class,
         -name_order, -name_family, -name_genus)

n_with_family_final <- sum(!is.na(echino_wide$best_family))
n_with_order_final  <- sum(!is.na(echino_wide$best_order))
cat(sprintf("Records with best_family resolved (aphiaID + name fallback combined): %d of %d (%.1f%%)\n",
            n_with_family_final, nrow(echino_wide), 100 * n_with_family_final / nrow(echino_wide)))
cat(sprintf("Records with best_order resolved (aphiaID + name fallback combined):  %d of %d (%.1f%%)\n",
            n_with_order_final, nrow(echino_wide), 100 * n_with_order_final / nrow(echino_wide)))

# -----------------------------------------------------------------------------
# Fill echino_class for records at a Class-bucket sub-rank (Subclass,
# Superclass, Infraclass, Subterclass, Cohort - e.g. "Echinacea", a
# Subterclass) that Section 6j deliberately left as NA rather than guessing,
# since accepted_name_final for these is not itself a valid class name (see
# the tightened worms_rank_master == "Class" condition in Section 6j-2 STEP
# 3). worms_class is now available from the aphiaID/name lookups above and
# gives the true parent class directly - this is a genuine fill of a real
# gap, not a "disagreement correction" (echino_class was NA here, not wrong),
# so it's kept separate from the mismatch-and-correct step below.
# -----------------------------------------------------------------------------

n_before_class_fill <- sum(is.na(echino_wide$echino_class) & echino_wide$taxonomic_resolution_level == "Class")

echino_wide <- echino_wide %>%
  mutate(
    echino_class_source = if_else(
      is.na(echino_class) & taxonomic_resolution_level == "Class" & !is.na(worms_class),
      "inferred_from_worms_class_subrank",
      echino_class_source
    ),
    echino_class = if_else(
      is.na(echino_class) & taxonomic_resolution_level == "Class" & !is.na(worms_class),
      worms_class,
      echino_class
    )
  )

n_after_class_fill <- sum(is.na(echino_wide$echino_class) & echino_wide$taxonomic_resolution_level == "Class")
cat(sprintf("echino_class filled from worms_class for Class-bucket sub-ranks (Subterclass etc.): %d of %d\n",
            n_before_class_fill - n_after_class_fill, n_before_class_fill))

# -----------------------------------------------------------------------------
# Cross-check: WoRMS' own class for each aphiaID/name should agree with
# echino_class (resolved independently in 6j from source class__* columns +
# the manual lookup table). Disagreement flags either a bad aphiaID/name
# match or a gap in the manual lookup table - not expected to fire often.
#
# MUST run here, after the addendum above, not right after the aphiaID-only
# join - worms_class gets filled two ways (the aphiaID join, then the
# name-based fallback for whatever the aphiaID join missed), and running
# this check between those two steps means records whose worms_class only
# arrives via the name fallback are invisible to it (worms_class is still
# NA at that point), so their echino_class never gets checked at all. This
# is exactly why a prior run reported "0 mismatches" while 24 known-bad
# records were confirmed still wrong by direct inspection of the saved CSV -
# their worms_class simply hadn't arrived yet when the check ran.
# -----------------------------------------------------------------------------

class_mismatch <- echino_wide %>%
  filter(!is.na(worms_class), !is.na(echino_class), worms_class != echino_class)

cat("\nRecords where WoRMS class disagrees with echino_class:", nrow(class_mismatch), "\n")
if (nrow(class_mismatch) > 0) {
  write_csv(
    class_mismatch %>% select(record_key, primary_source, accepted_name_final,
                              aphiaID, echino_class, worms_class),
    "echino_class_worms_mismatch.csv"
  )
  cat("Mismatches written to echino_class_worms_mismatch.csv for review.\n")
}

# Manual review of a prior batch of these mismatches (24 records, all
# ALA_echinodermata, record_key prefix NMR-SS...) confirmed accepted_name_final
# and worms_class agree with each other and are correct; echino_class was the
# wrong one, apparently carrying a static incorrect value from the source's
# raw class__ALA_echinodermata field for that whole record batch. Treating
# the two independently-agreeing, name-derived values as authoritative over
# the source's raw (evidently, at least sometimes, corrupted) class field.
if (nrow(class_mismatch) > 0) {
  echino_wide <- echino_wide %>%
    mutate(
      echino_class = if_else(
        record_key %in% class_mismatch$record_key,
        worms_class,
        echino_class
      ),
      echino_class_source = if_else(
        record_key %in% class_mismatch$record_key,
        "corrected_from_worms_aphiaID_source_class_field_error",
        echino_class_source
      )
    )
  cat("echino_class corrected via WoRMS cross-check for", nrow(class_mismatch),
      "records (source class__ALA_echinodermata field error, see\n",
      "echino_class_worms_mismatch.csv for the record_keys and original values).\n")
}

# Hemiaster (AphiaID 123423) still comes back with no family/order - or
# kingdom/phylum/class - from either lookup above: a fossil-adjacent genus,
# incomplete in WoRMS' classification tree for this specific AphiaID even
# though the aphiaID itself resolves fine (echino_class = Echinoidea was
# already correct via Section 6j's own resolution, independent of this gap).
# All five ranks confirmed by direct inspection of
# https://www.marinespecies.org/aphia.php?p=taxdetails&id=123423
# (Kingdom Animalia, Phylum Echinodermata, Class Echinoidea, Order
# Spatangoida, Family Hemiasteridae) and hardcoded here for this one record
# (is_extinct_genus, n=1) rather than left NA, since the correct values are
# now directly verified against the source rather than inferred.
echino_wide <- echino_wide %>%
  mutate(
    worms_kingdom = if_else(is_extinct_genus, "Animalia", worms_kingdom),
    worms_phylum  = if_else(is_extinct_genus, "Echinodermata", worms_phylum),
    worms_class   = if_else(is_extinct_genus, "Echinoidea", worms_class),
    best_family   = if_else(is_extinct_genus, "Hemiasteridae", best_family),
    best_order    = if_else(is_extinct_genus, "Spatangoida", best_order)
  )

cat("\n=== Section 6k complete ===\n")

# =============================================================================
# SECTION 7: FINAL SAVE + SUMMARY DIAGNOSTICS
# =============================================================================

cat("\n== Section 7: Final save + summary ==\n\n")

# -----------------------------------------------------------------------------
# 7a. Pipeline record-count audit trail
# -----------------------------------------------------------------------------

pipeline_audit <- tribble(
  ~stage, ~records, ~notes,
  "Pipeline integration output", 459956, "13 sources merged (ALA×3, CSIRO×2, GBIF×2, OBIS, OZCAM, iDigBio, MTQ_CIDARIS, CIDARIS_QMT, AM_direct); before catalogNumber case-insensitivity fix",
  "Section 0: case-duplicate record_key removal", 50027, "0 case duplicates found (already resolved upstream)",
  "Section 1: coordinate fixes + geographic scope", 50014, "13 records removed: 3 NSW/WA outliers (Fix 2) + 10 scope exclusions (Darwin NT, inland fungus, 4x Disaster genus inland NSW, 3x York Sound WA)",
  "Section 1b: coordinate recovery", 50014, "No records removed; 555 coordinates recovered — 48 idigbio:geoPoint (verified exact) + 507 locality gazetteer (estimated, spot-checked); 94 additional iDigBio recovery records skipped (already had coordinates from integration pipeline); 273 records remain without coordinates",
  "Section 2: country fixes", 50014, "No records removed; country field standardised",
  "Section 3: year fixes", 50014, "No records removed; 4 impossible years + 40 historical placeholders corrected",
  "Section 4: basis of record review", 43512, "6,502 removed: fossil/machine (45), living specimens (4,642), material samples (1,651) + 3 KX PreservedSpecimen GBIF duplicates of QM G22502, IMOS zooplankton (161); basisOfRecord standardised to Darwin Core camelCase",
  "Section 5: depth resolution", 43512, "No records removed; 24,669 records with depth (56.7%); 688 recovered via: event-sibling imputation (645) and reef flat/intertidal inference (43)",
  "Section 6.0: non-echinoderm contaminant removal", 43470, "42 removed: Asterina eupomatiae fungus (1), Anoplura lice (2), Rhodostoma Mollusca (1), ANWC mislabelled mammals (38)",
  "Section 6b fix applied (case-insensitive catalogNumber)", 50027, "37,741 case-mismatched duplicate rows collapsed",
  "Section 6e: name normalisation + WoRMS resolution (conflicted names only)", 43470, "accepted_name assigned to 43,352 records; three-criterion weighted vote (WoRMS status > taxonomic fineness > source trust) replaces old column-order tie-breaking — resolved via WoRMS (28,602) and manual (30)",
  "Section 6h: full WoRMS validation of every accepted_name", 43470, "100% WoRMS coverage; 28 names corrected via verified_valid_names (incl. Iconometra dormani — WoRMS-unresolvable, genus retained); aa- class-level AM records assigned accepted_name; accepted_name_final propagated",
  "Section 6i: post-audit fixes", 43470, "taxonomic_resolution_level recomputed from accepted_name_final word structure; is_extinct_genus flag added for Hemiaster",
  "Section 6j: echino_class resolution", 43470, "echino_class resolved for 42,642 records; 828 remain NA (Phylum-only, expected)",
  "Section 6k: higher classification from WoRMS (family, order)", 43470, "best_family/best_order resolved via aphiaID lookup, then a name-based WoRMS fallback for records with no aphiaID (final combined coverage printed at runtime - see Section 6k output); echino_class cross-checked against worms_class after both lookups complete, with any disagreement corrected (source class__ALA_echinodermata field error - see echino_class_worms_mismatch.csv)"
)

write_csv(pipeline_audit, "pipeline_summary.csv")
cat("Pipeline audit trail saved: pipeline_summary.csv\n\n")
print(pipeline_audit, n = Inf)


# -----------------------------------------------------------------------------
# 7b. Final dataset overview
# -----------------------------------------------------------------------------

cat("\n── Final dataset overview ──\n")
cat("Total records:", nrow(echino_wide), "\n")
cat("Total columns:", ncol(echino_wide), "\n")

year_cols <- names(echino_wide)[str_detect(names(echino_wide), "^year__")]
all_years_final <- echino_wide %>%
  select(all_of(year_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.integer(.x)))) %>%
  unlist() %>%
  na.omit()
cat("Date range (year), across all sources:", min(all_years_final), "-", max(all_years_final), "\n\n")

cat("── Geographic coverage ──\n")
echino_wide %>% count(country, sort = TRUE) %>% print()

cat("\n── Depth coverage ──\n")
echino_wide %>%
  summarise(
    n_total = n(),
    n_with_depth = sum(has_depth),
    pct_with_depth = round(mean(has_depth) * 100, 1)
  ) %>% print()

cat("\n── Depth zone distribution ──\n")
echino_wide %>%
  count(depth_zone, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = Inf)

cat("\n── Data quality flag distribution ──\n")
echino_wide %>%
  count(obs_quality_flag) %>%
  mutate(
    label = case_when(
      obs_quality_flag == 1L ~ "Flag 1 - Specimen-based + GBRSBD",
      obs_quality_flag == 2L ~ "Flag 2 - Formal scientific survey",
      obs_quality_flag == 3L ~ "Flag 3 - Citizen science",
      TRUE ~ "Unassigned"
    ),
    pct = round(n / sum(n) * 100, 1)
  ) %>% arrange(obs_quality_flag) %>% print()

cat("\n── Taxonomic resolution level ──\n")
echino_wide %>%
  count(taxonomic_resolution_level, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = Inf)

cat("\n── Final WoRMS validation coverage ──\n")
cat(sprintf("%d of %d records (%.1f%%) have a confirmed WoRMS status\n",
            sum(!is.na(echino_wide$worms_status_master)), nrow(echino_wide),
            100 * mean(!is.na(echino_wide$worms_status_master))))
echino_wide %>%
  count(worms_status_master, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = Inf)

cat("\n── Extinct genus flag ──\n")
echino_wide %>% count(is_extinct_genus) %>% print()


# -----------------------------------------------------------------------------
# 7c. Cross-tabulated completeness summary
# -----------------------------------------------------------------------------
# NOTE: these crosstabs need to exclude records whose depth range CROSSES
# a zone boundary (depth_zone == "Spans X to Y"), not records flagged
# is_straddler. is_straddler only measures depth-range width (> 100m) and
# is unrelated to whether a record's zone assignment is actually
# ambiguous — the two flags overlap on only 10 of 41 boundary-crossing
# records. Using is_straddler here would incorrectly keep 31 genuinely
# ambiguous records in a single zone while dropping 28 wide-but-unambiguous
# ones. See `column_documentation.csv` for both definitions.

echino_wide <- echino_wide %>%
  mutate(crosses_zone_boundary = str_starts(coalesce(depth_zone, ""), "Spans"))

cat("\nRecords excluded from zone-stratified crosstabs (crosses_zone_boundary):",
    sum(echino_wide$crosses_zone_boundary), "\n")

cat("\n── Cross-tab: depth_zone x taxonomic_resolution_level ──\n")

completeness_crosstab <- echino_wide %>%
  filter(!crosses_zone_boundary) %>%
  count(depth_zone, taxonomic_resolution_level) %>%
  group_by(depth_zone) %>%
  mutate(
    zone_total = sum(n),
    pct_of_zone = round(n / zone_total * 100, 1)
  ) %>%
  ungroup() %>%
  arrange(depth_zone, desc(n))

print(completeness_crosstab, n = Inf)

cat("\n── Species-level identification rate BY depth zone ──\n")
species_rate_by_zone <- echino_wide %>%
  filter(!crosses_zone_boundary) %>%
  group_by(depth_zone) %>%
  summarise(
    n = n(),
    n_species_level = sum(taxonomic_resolution_level == "Species"),
    pct_species_level = round(mean(taxonomic_resolution_level == "Species") * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(match(depth_zone, c("Continental Shelf", "Upper Slope", "Lower Slope",
                              "Abyssal", "No depth data")))
print(species_rate_by_zone, n = Inf)

cat("\n── Species-level identification rate BY data quality flag ──\n")
species_rate_by_flag <- echino_wide %>%
  group_by(obs_quality_flag) %>%
  summarise(
    n = n(),
    pct_species_level = round(mean(taxonomic_resolution_level == "Species") * 100, 1),
    .groups = "drop"
  )
print(species_rate_by_flag, n = Inf)

cat("\n── Full three-way cross-tab: depth_zone x obs_quality_flag x species-level rate ──\n")
full_crosstab <- echino_wide %>%
  filter(!crosses_zone_boundary) %>%
  group_by(depth_zone, obs_quality_flag) %>%
  summarise(
    n = n(),
    pct_species_level = round(mean(taxonomic_resolution_level == "Species") * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(match(depth_zone, c("Continental Shelf", "Upper Slope", "Lower Slope",
                              "Abyssal", "No depth data")),
          obs_quality_flag)
print(full_crosstab, n = Inf)

write_csv(completeness_crosstab, "completeness_summary_full_crosstab.csv")
write_csv(species_rate_by_zone, "completeness_summary_by_depth_zone.csv")
write_csv(species_rate_by_flag, "completeness_summary_by_quality_flag.csv")
write_csv(full_crosstab, "completeness_summary_depth_x_flag.csv")

cat("\nCompleteness summary tables saved.\n")


# -----------------------------------------------------------------------------
# 7d. Final column inventory
# -----------------------------------------------------------------------------
# FIX APPLIED: "27 corrected in 6h/6i" (was "21 corrected in 6h" - the two
# post-audit corrections, ochroleucus and intermedius, are now included in
# the count). Two rows added for columns that exist in the final dataset
# but were missing from this documentation table: depth_imputed_flag
# (created in Section 5g/5h) and is_extinct_genus (created in Section 6i).
# -----------------------------------------------------------------------------

key_columns <- tribble(
  ~column, ~description, ~section_created,
  "record_key", "Unique specimen identifier (case-insensitive, catalogNumber-reconciled)", "Pipeline 6/6b",
  "primary_source", "Source dataset assigned as this record's primary source, used as the tiebreaker/priority basis for source-priority resolution throughout the pipeline (coordinates, basisOfRecord, etc.)", "Pipeline integration",
  "basisOfRecord", "Raw, per-source basisOfRecord value prior to resolution; see best_basisOfRecord (Section 4) for the resolved consensus value used in analysis", "Pipeline integration",
  "completeness_score", "Composite score reflecting how many Darwin Core fields were populated for this record at the point of source integration (Layer A)", "Pipeline integration",
  "n_sources", "Number of distinct source datasets that reported this record (i.e. how many sources merged into this record_key)", "Pipeline integration",
  "sources_all", "List of all source datasets that reported this record, as a single delimited string", "Pipeline integration",
  "best_latitude", "Resolved decimal latitude", "Section 1a",
  "best_longitude", "Resolved decimal longitude", "Section 1a",
  "coord_conflict", "TRUE if sources disagreed on coordinates by >1 degree", "Section 1a",
  "coord_recovered_flag", "TRUE if coordinates came from recovery (iDigBio geoPoint or gazetteer) rather than original pipeline", "Section 1b",
  "coord_recovery_method", "idigbio_geopoint_recovery (verified exact) or locality_gazetteer_geocoding (estimated, spot-checked)", "Section 1b",
  "coord_uncertainty_m", "Assigned coordinate uncertainty radius in metres for gazetteer-estimated records; NA for exact geoPoint recoveries", "Section 1b",
  "coord_land_qc_flag", "Distance-to-coast + coordinateUncertaintyInMeters based QC flag for georeferenced records that fall on land. not_on_land = no flag (includes both confirmed marine coordinates and the 226 records with no coordinates at all); on_land_near_coast_or_island = within 2km of coast, consistent with genuine small-island/reef-station localities (kept as-is); on_land_high_uncertainty = coordinateUncertaintyInMeters >= 10km; on_land_far_inland_review = >5km inland with no uncertainty flag; on_land_moderate_distance_review = ambiguous, needs manual review", "Section 1d",
  "coord_qc_exclude_recommended", "TRUE only for on_land_high_uncertainty / on_land_far_inland_review (842 records, 1.94%). Convenience filter column - use this, not coord_land_qc_flag directly, to exclude low-confidence coordinates from an analysis", "Section 1d",
  "country", "Standardised country (Australia / Papua New Guinea)", "Section 2",
  "is_historical_collection", "TRUE for records with placeholder collection years flagged by the source database", "Section 3",
  "best_year", "Resolved collection year via majority-vote consensus across year__* source columns (mirrors best_latitude/best_longitude resolution). 4 records (0.01%) had genuine cross-source disagreement, each resolved by majority vote rather than arbitrary column order.", "Section 3",
  "best_eventDate", "Resolved collection event date via source-priority consensus across eventDate__* source columns (highest-priority source with a non-blank value wins, not majority vote - date strings rarely match verbatim across sources even when they agree). Kept to date-only (YYYY-MM-DD); time component stripped as largely placeholder data.", "Section 3",
  "obs_quality_flag", "1=specimen-based/GBRSBD, 2=formal survey, 3=citizen science", "Section 4",
  "dataResourceName", "Raw source-reported project/dataset name (e.g. a specific ALA or OBIS sub-collection); used to assign obs_quality_flag for human-observation records (Table 3/4)", "Section 4",
  "best_min_depth", "Resolved minimum depth (m); GBIF-flattening artifact and source min/max-swap corrected", "Section 5",
  "best_max_depth", "Resolved maximum depth (m); GBIF-flattening artifact and source min/max-swap corrected", "Section 5",
  "depth_zone", "Ecologically-grounded depth zone", "Section 5",
  "depth_median", "Midpoint of best_min_depth/best_max_depth - primary continuous depth variable", "Section 5",
  "depth_uncertainty", "best_max_depth - best_min_depth (range width); always non-negative after post-audit fix", "Section 5",
  "is_straddler", "TRUE if depth_uncertainty > 100m (depth-range WIDTH only - not zone-related; see crosses_zone_boundary for zone-stratified analysis)", "Section 5",
  "depth_imputed_flag", "TRUE if depth was recovered via event-sibling imputation or reef-flat/intertidal text inference, rather than directly measured/reported", "Section 5g/5h",
  "name_conflict_real", "TRUE if sources genuinely disagreed on scientific name after normalisation", "Section 6b",
  "accepted_name", "Resolved scientific name from Section 6e/6h (conflict-resolution pass)", "Section 6e",
  "accepted_name_worms", "TRUE if accepted_name was resolved via WoRMS in Section 6e ONLY (conflicted names only - NOT a complete validation flag, see worms_status_master instead)", "Section 6e",
  "accepted_name_manual", "TRUE if accepted_name was resolved via manual taxonomic correction in Section 6e", "Section 6e",
  "winning_source", "Source that won the three-criterion weighted vote for this record's accepted_name (status_rank > fineness_rank > src_rank). Useful for auditing which dataset drove the final taxonomic identification.", "Section 6e",
  "worms_status_master", "Complete WoRMS status for every record (accepted/unaccepted/synonym/etc) - covers 100% of records", "Section 6h",
  "worms_rank_master", "Exact WoRMS taxonomic rank for accepted_name_final (e.g. \"Subterclass\", \"Species\"), uncollapsed - covers 100% of records. Use this rather than taxonomic_resolution_level when the precise rank is needed, not the simplified category.", "Section 6h",
  "accepted_name_final", "FINAL authoritative name - current WoRMS-valid combination for ALL records, including the 135 records (31 distinct names) corrected via verified_valid_names. Use this column, not accepted_name, for analysis.", "Section 6h/6i",
  "taxonomic_resolution_level", "Rank bucketed into eight categories (Phylum/Class/Order/Family/Genus/Subgenus/Species/Subspecies), primarily from worms_rank_master via a hierarchy-based mapping (worms_rank_to_bucket) that collapses related WoRMS ranks into their parent category (e.g. Subclass, Superclass, Subterclass all map to \"Class\"); falls back to accepted_name_final word structure only where worms_rank_master is unavailable.", "Section 6e/6i/6i-fix",
  "is_extinct_genus", "TRUE for Hemiaster (extinct genus, intentionally retained - basisOfRecord=PreservedSpecimen not Fossil, genuine modern collection)", "Section 6i (post-audit)",
  "echino_class", "Resolved echinoderm class (Asteroidea/Ophiuroidea/Echinoidea/Holothuroidea/Crinoidea), cross-checked against worms_class every run and corrected where they disagreed (see echino_class_source). NA only for genuine Phylum-only identifications. worms_class (below) is the recommended column to cite as authoritative, since it traces directly to WoRMS rather than a source field.", "Section 6j/6k",
  "echino_class_source", "original_source_field if resolved directly from source class__* columns, inferred_from_taxon_name if recovered via taxon->class lookup, corrected_from_worms_aphiaID_source_class_field_error if corrected via the Section 6k WoRMS cross-check, NA if still unresolved (Phylum-only)", "Section 6j/6k",
  "worms_class", "Class resolved directly from WoRMS via aphiaID or accepted-name lookup, independent of any source's own class__* field. Identical to echino_class for every record after the Section 6k correction step - the recommended column to cite for class-level analysis.", "Section 6k",
  "worms_kingdom", "Kingdom (Animalia for all records) from the same WoRMS lookup as worms_class/best_family/best_order.", "Section 6k",
  "worms_phylum", "Phylum (Echinodermata for all records) from the same WoRMS lookup as worms_class/best_family/best_order.", "Section 6k",
  "best_family", "Family, resolved via WoRMS lookup on aphiaID (not from unreconciled per-source family__* columns). NA where aphiaID is missing or WoRMS has no family on record for that AphiaID.", "Section 6k",
  "best_order", "Order, resolved via WoRMS lookup on aphiaID (not from unreconciled per-source order__* columns). NA where aphiaID is missing or WoRMS has no order on record for that AphiaID.", "Section 6k",
  "crosses_zone_boundary", "TRUE if best_min_depth and best_max_depth fall in different depth_zone values (depth_zone == 'Spans X to Y'). Use THIS, not is_straddler, to exclude records from zone-stratified analyses (TIa, completeness crosstabs etc) - the two flags overlap on only 10/41 records", "Section 7c"
)

write_csv(key_columns, "column_documentation.csv")
print(key_columns, n = Inf)
cat("\nColumn documentation saved: column_documentation.csv\n")


# -----------------------------------------------------------------------------
# 7e. Final save
# -----------------------------------------------------------------------------

cat("\n== Saving final files ==\n\n")

write_csv(echino_wide, "echino_wide.csv")
cat("Saved: echino_wide.csv (", nrow(echino_wide), "records,", ncol(echino_wide), "columns )\n")

echino_wide_depth <- echino_wide %>% filter(has_depth)
write_csv(echino_wide_depth, "echino_wide_depth.csv")
cat("Saved: echino_wide_depth.csv (", nrow(echino_wide_depth), "records with resolved depth )\n")

cat("\n=== PHASE 1 POST-PROCESSING PIPELINE COMPLETE ===\n")
cat("Final dataset:", nrow(echino_wide), "records,", ncol(echino_wide), "columns\n")
cat("100% of records have a confirmed WoRMS taxonomic status (accepted_name_final / worms_status_master)\n")
cat("Files saved to working directory:\n")
cat("  echino_wide.csv                              - full final dataset\n")
cat("  echino_wide_depth.csv                        - depth-bearing subset\n")
cat("  pipeline_summary.csv                         - record-count audit trail\n")
cat("  name_corrections_for_review.csv              - Section 6d manual corrections\n")
cat("  section6h_newly_resolved_names.csv           - Section 6h newly-validated names\n")
cat("  all_non_accepted_names_summary.csv           - every non-accepted name, with current WoRMS-valid name\n")
cat("  worms_status_reference.csv                   - plain-language meaning of each WoRMS status\n")
cat("  worms_status_distribution.csv                - status distribution among Section 6e names\n")
cat("  completeness_summary_full_crosstab.csv       - depth x taxonomic resolution crosstab\n")
cat("  completeness_summary_by_depth_zone.csv       - species-ID rate by depth zone\n")
cat("  completeness_summary_by_quality_flag.csv     - species-ID rate by data quality flag\n")
cat("  completeness_summary_depth_x_flag.csv        - three-way crosstab\n")
cat("  column_documentation.csv                     - key column reference\n")

# =============================================================================
# SECTION 8: DARWIN CORE ARCHIVE EXPORT (occurrence core + eMoF extension)
#            PUBLIC RELEASE — OBIS / ALA / GBIF READY
# =============================================================================
# All source institutions have confirmed redistribution permission
# in writing (Australian Museum, Queensland Museum Tropics). This
# archive is the public deposit version.
#
# Outputs:
#   dwca_public_release/occurrence.txt        (occurrence core)
#   dwca_public_release/measurementorfact.txt (eMoF extension)
#   dwca_public_release/meta.xml
#   dwca_public_release/eml.xml
#   echinoderm_dwca_public_release.zip
#
# PLACEHOLDERS TO FILL IN BEFORE PUBLICATION (search for "<<< FILL IN >>>"):
#   * ORCID iDs for Alastair Birtles and Sue-Ann Watson
#   * Scientific Data DOI (in `references`, EML citation, and abstract)
#   * Dataset DOI          (in `datasetID`, EML citation)
#   * IPT / repository download URL (in EML <distribution>)
#   * GRSciColl institution/collection UUIDs (institutionID, collectionID)
#   * Grant / funding information (EML <project>)
# =============================================================================

cat("\n== Section 8: Darwin Core Archive export (public release) ==\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(purrr)
  library(tidyr)
})

dwc_dir <- "dwca_public_release"
dir.create(dwc_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 8.0  Helpers
# -----------------------------------------------------------------------------

coalesce_field_by_priority <- function(df, field_prefix) {
  cols <- names(df)[str_detect(names(df), paste0("^", field_prefix, "__"))]
  ordered_cols <- intersect(paste0(field_prefix, "__", SOURCE_PRIORITY), cols)
  ordered_cols <- c(ordered_cols, setdiff(cols, ordered_cols))
  if (length(ordered_cols) == 0) return(rep(NA_character_, nrow(df)))
  mat <- df %>% select(all_of(ordered_cols)) %>% as.matrix()
  apply(mat, 1, function(row_vals) {
    vals <- row_vals[!is.na(row_vals) & str_trim(coalesce(row_vals, "")) != ""]
    if (length(vals) == 0) return(NA_character_)
    vals[1]
  })
}

# OBIS strongly recommends WoRMS URN LSID form for taxonID / scientificNameID.
worms_urn <- function(aphia) {
  a <- suppressWarnings(as.integer(aphia))
  ifelse(is.na(a), "", paste0("urn:lsid:marinespecies.org:taxname:", a))
}

# Strip characters that would break a quote-less tab-delimited archive.
safe_text <- function(x) {
  x %>%
    str_replace_all("[\t\r\n]+", " ") %>%
    str_replace_all("\"", "'")
}

# ISO 8601 date check (accepts YYYY, YYYY-MM, YYYY-MM-DD, YYYY-MM-DD/YYYY-MM-DD)
is_iso_date <- function(x) {
  x <- coalesce(x, "")
  x == "" | str_detect(x, "^\\d{4}(-\\d{2}(-\\d{2})?)?(/\\d{4}(-\\d{2}(-\\d{2})?)?)?$")
}

# -----------------------------------------------------------------------------
# 8.1  Locate the AphiaID column (varies by pipeline section)
# -----------------------------------------------------------------------------
aphia_candidates <- c("accepted_aphiaID", "aphiaID", "worms_aphiaID")
aphia_col <- aphia_candidates[aphia_candidates %in% names(echino_wide)]
if (length(aphia_col) == 0) {
  warning("No AphiaID column found. scientificNameID / taxonID will be blank.")
  echino_wide$._aphia <- NA_integer_
} else {
  echino_wide$._aphia <- echino_wide[[aphia_col[1]]]
  cat("Using", aphia_col[1], "for scientificNameID / taxonID / acceptedNameUsageID.\n")
}

# -----------------------------------------------------------------------------
# 8.2  Coalesce per-source columns not otherwise resolved
# -----------------------------------------------------------------------------
cat("Coalescing locality / typeStatus / identifiedBy / dateIdentified...\n")
dwc_locality       <- coalesce_field_by_priority(echino_wide, "locality")
dwc_typeStatus     <- coalesce_field_by_priority(echino_wide, "typeStatus")
dwc_identifiedBy   <- coalesce_field_by_priority(echino_wide, "identifiedBy")
dwc_dateIdentified <- coalesce_field_by_priority(echino_wide, "dateIdentified")

resolution_lower <- str_to_lower(echino_wide$taxonomic_resolution_level)

# Parse "Genus (Subgenus) species subspecies" without corrupting epithet columns
name_match <- str_match(
  echino_wide$accepted_name_final,
  "^(\\S+)(?:\\s+\\(([^)]+)\\))?(?:\\s+(\\S+))?(?:\\s+(\\S+))?$"
)
# [,2] genus  [,3] subgenus  [,4] specificEpithet  [,5] infraspecificEpithet

# -----------------------------------------------------------------------------
# 8.3  Build the occurrence core
# -----------------------------------------------------------------------------
cat("Building occurrence core...\n")

occ <- tibble(
  # ---- record identifiers & provenance ----
  occurrenceID        = paste0("JCU-ECHINODERM:", echino_wide$record_key),
  modified            = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  language            = "eng",
  license             = "http://creativecommons.org/licenses/by/4.0/legalcode",
  rightsHolder        = "James Cook University",
  # <<< FILL IN: Scientific Data DOI once assigned >>>
  references          = "https://doi.org/[Scientific-Data-DOI]",
  # <<< FILL IN: Dataset DOI once minted (Zenodo, OBIS, or ALA) >>>
  datasetID           = "https://doi.org/[Dataset-DOI]",
  datasetName         = "Integrated echinoderm occurrence dataset for northeast Australia and the adjacent Coral Sea",
  
  # ---- institution / collection ----
  # <<< OPTIONAL: fill with GRSciColl UUIDs where known >>>
  institutionID       = "",
  collectionID        = "",
  institutionCode     = str_to_upper(str_trim(coalesce(echino_wide$institutionCode, ""))),
  collectionCode      = echino_wide$collectionCode,
  ownerInstitutionCode = str_to_upper(str_trim(coalesce(echino_wide$institutionCode, ""))),
  
  # ---- record type ----
  type = if_else(coalesce(echino_wide$best_basisOfRecord, "") == "PreservedSpecimen",
                 "PhysicalObject", ""),
  basisOfRecord       = echino_wide$best_basisOfRecord,
  informationWithheld = "",
  dataGeneralizations = if_else(
    coalesce(echino_wide$coord_recovered_flag, FALSE) &
      coalesce(echino_wide$coord_recovery_method, "") == "locality_gazetteer_geocoding",
    "Coordinates estimated from locality text via named-place gazetteer; see coordinateUncertaintyInMeters and georeferenceProtocol.",
    ""
  ),
  
  # ---- occurrence-level ----
  catalogNumber       = echino_wide$catalogNumber,
  recordNumber        = "",
  recordedBy          = echino_wide$recordedBy,
  individualCount     = as.character(echino_wide$individualCount),
  organismQuantity    = as.character(echino_wide$individualCount),
  organismQuantityType = if_else(
    !is.na(echino_wide$individualCount) &
      as.character(echino_wide$individualCount) != "",
    "individuals", ""
  ),
  occurrenceStatus    = coalesce(str_to_lower(echino_wide$occurrenceStatus), "present"),
  occurrenceRemarks   = paste0(
    "Integrated record: n_sources=", echino_wide$n_sources,
    " (", echino_wide$sources_all, ").",
    " Depth zone: ", echino_wide$depth_zone, ".",
    " Data-quality flag: ", echino_wide$obs_quality_flag, "."
  ),
  preparations = if_else(coalesce(echino_wide$best_basisOfRecord, "") == "PreservedSpecimen",
                         "preserved specimen", ""),
  
  # ---- event ----
  eventDate           = echino_wide$best_eventDate,
  year                = as.character(echino_wide$best_year),
  samplingProtocol    = "",
  
  # ---- location ----
  country             = echino_wide$country,
  countryCode         = echino_wide$countryCode,
  stateProvince       = str_to_title(str_trim(coalesce(echino_wide$stateProvince, ""))),
  locality            = dwc_locality,
  minimumDepthInMeters = as.character(echino_wide$best_min_depth),
  maximumDepthInMeters = as.character(echino_wide$best_max_depth),
  decimalLatitude     = as.character(echino_wide$best_latitude),
  decimalLongitude    = as.character(echino_wide$best_longitude),
  geodeticDatum       = "EPSG:4326",
  coordinateUncertaintyInMeters = as.character(
    coalesce(echino_wide$coord_uncertainty_m, echino_wide$coordinateUncertaintyInMeters)
  ),
  georeferenceSources = case_when(
    coalesce(echino_wide$coord_recovery_method, "") == "idigbio_geopoint_recovery"    ~ "iDigBio geoPoint JSON field",
    coalesce(echino_wide$coord_recovery_method, "") == "locality_gazetteer_geocoding" ~ "GeoNames gazetteer; locality-text lookup",
    TRUE ~ ""
  ),
  georeferenceProtocol = case_when(
    coalesce(echino_wide$coord_recovery_method, "") == "idigbio_geopoint_recovery"    ~ "Coordinate recovered from source geoPoint field",
    coalesce(echino_wide$coord_recovery_method, "") == "locality_gazetteer_geocoding" ~ "Locality-text geocoding against named-place gazetteer with assigned uncertainty radius",
    TRUE ~ ""
  ),
  georeferenceRemarks = if_else(coalesce(echino_wide$coord_recovered_flag, FALSE),
                                "Coordinate recovered during post-processing; not in original source record.",
                                ""),
  
  # ---- identification ----
  identifiedBy        = dwc_identifiedBy,
  dateIdentified      = dwc_dateIdentified,
  typeStatus          = dwc_typeStatus,
  
  # ---- taxon ----
  scientificName      = echino_wide$accepted_name_final,
  scientificNameAuthorship = echino_wide$scientificNameAuthorship,
  taxonID             = worms_urn(echino_wide$._aphia),
  scientificNameID    = worms_urn(echino_wide$._aphia),
  acceptedNameUsageID = worms_urn(echino_wide$._aphia),
  taxonRank           = resolution_lower,
  taxonomicStatus     = str_to_lower(coalesce(echino_wide$worms_status_master, "")),
  nomenclaturalCode   = "ICZN",
  kingdom             = "Animalia",
  phylum              = "Echinodermata",
  class               = coalesce(echino_wide$worms_class, echino_wide$echino_class),
  order               = coalesce(echino_wide$best_order, ""),
  family              = coalesce(echino_wide$best_family, ""),
  genus               = if_else(resolution_lower %in% c("species", "subspecies", "genus"),
                                coalesce(name_match[, 2], ""), ""),
  subgenus            = coalesce(name_match[, 3], ""),
  specificEpithet     = if_else(resolution_lower %in% c("species", "subspecies"),
                                coalesce(name_match[, 4], ""), ""),
  infraspecificEpithet = echino_wide$infraspecificEpithet,
  
  # ---- other ----
  fieldNumber          = "",
  associatedReferences = "",
  # dynamicProperties keeps only fields without a clean DwC mapping.
  # obs_quality_flag / depth_zone / n_sources / sources_all live in
  # occurrenceRemarks and in the eMoF extension below.
  dynamicProperties    = paste0(
    '{"depth_zone":"', coalesce(echino_wide$depth_zone, ""),
    '","is_straddler":', str_to_lower(as.character(coalesce(echino_wide$is_straddler, FALSE))),
    ',"crosses_zone_boundary":', str_to_lower(as.character(coalesce(echino_wide$crosses_zone_boundary, FALSE))),
    ',"coord_land_qc_flag":"', coalesce(echino_wide$coord_land_qc_flag, ""),
    '"}'
  )
) %>%
  mutate(across(everything(), ~ replace_na(as.character(.x), "")))

# -----------------------------------------------------------------------------
# 8.4  Pre-publication filters + validation
# -----------------------------------------------------------------------------
cat("Applying pre-publication filters...\n")

# No records are dropped here. coord_qc_exclude_recommended is preserved as
# metadata (dynamicProperties/coord_land_qc_flag) for downstream users to
# filter themselves - never used to remove records from the archive.
# The "duplicate registration" text check was investigated and removed:
# the note only ever appears in locality__iDigBio (never in an actual
# notes/remarks field, and never confirmed by the source institution's own
# system - see G12229/G12229.2 cross-check against QM's own collections
# database, which shows no duplicate flag and different preparations for
# each). Treated as unreliable and not used to drop or flag records.

echino_wide_kept <- echino_wide

cat("  Records kept: all", nrow(occ), "- no exclusion filters applied.\n")

# Blank out non-ISO event dates (keep the record; drop only the date)
occ$eventDate[!is_iso_date(occ$eventDate)] <- ""

# Strip tabs / newlines / stray double-quotes from every field
occ <- occ %>% mutate(across(everything(), safe_text))

# Final assertions
stopifnot(all(!is.na(occ$occurrenceID) & occ$occurrenceID != ""))
stopifnot(all(occ$basisOfRecord != ""))

cat("\n--- Section 8 QA ---\n")
cat("Records in archive           :", nrow(occ), "\n")
cat("Missing basisOfRecord        :", sum(occ$basisOfRecord == ""), "\n")
cat("Missing scientificName       :", sum(occ$scientificName == ""), "\n")
cat("Missing decimalLatitude      :", sum(occ$decimalLatitude == ""), "\n")
cat("Records with scientificNameID:", sum(occ$scientificNameID != ""), "\n")
cat("Records with recovered coords:", sum(occ$georeferenceSources != ""), "\n")
print(table(occ$basisOfRecord))
print(table(occ$taxonRank))

# -----------------------------------------------------------------------------
# 8.5  Write occurrence.txt
# -----------------------------------------------------------------------------
write_delim(occ, file.path(dwc_dir, "occurrence.txt"),
            delim = "\t", na = "", quote = "none", eol = "\n")
cat("Wrote", file.path(dwc_dir, "occurrence.txt"),
    ":", nrow(occ), "records x", ncol(occ), "fields\n")

# -----------------------------------------------------------------------------
# 8.6  Build MeasurementOrFact (eMoF) extension
# -----------------------------------------------------------------------------
cat("Building MeasurementOrFact (eMoF) extension...\n")

mof_rows <- bind_rows(
  tibble(
    occurrenceID     = occ$occurrenceID,
    measurementType  = "Observation quality flag",
    measurementValue = as.character(echino_wide_kept$obs_quality_flag),
    measurementUnit  = ""
  ),
  tibble(
    occurrenceID     = occ$occurrenceID,
    measurementType  = "Depth zone",
    measurementValue = as.character(echino_wide_kept$depth_zone),
    measurementUnit  = ""
  ),
  tibble(
    occurrenceID     = occ$occurrenceID,
    measurementType  = "Median depth",
    measurementValue = as.character(echino_wide_kept$depth_median),
    measurementUnit  = "m"
  ),
  tibble(
    occurrenceID     = occ$occurrenceID,
    measurementType  = "Depth range width",
    measurementValue = as.character(echino_wide_kept$depth_uncertainty),
    measurementUnit  = "m"
  ),
  tibble(
    occurrenceID     = occ$occurrenceID,
    measurementType  = "Number of contributing sources",
    measurementValue = as.character(echino_wide_kept$n_sources),
    measurementUnit  = ""
  )
) %>%
  filter(!is.na(measurementValue),
         str_trim(measurementValue) != "",
         measurementValue != "NA") %>%
  mutate(across(everything(), safe_text))

write_delim(mof_rows, file.path(dwc_dir, "measurementorfact.txt"),
            delim = "\t", na = "", quote = "none", eol = "\n")
cat("Wrote", file.path(dwc_dir, "measurementorfact.txt"),
    ":", nrow(mof_rows), "measurement rows for", nrow(occ), "occurrences\n")

# -----------------------------------------------------------------------------
# 8.7  Build meta.xml (core + eMoF extension)
# -----------------------------------------------------------------------------
cat("Building meta.xml...\n")

dwc_terms <- c(
  occurrenceID                  = "http://rs.tdwg.org/dwc/terms/occurrenceID",
  modified                      = "http://purl.org/dc/terms/modified",
  language                      = "http://purl.org/dc/terms/language",
  license                       = "http://purl.org/dc/terms/license",
  rightsHolder                  = "http://purl.org/dc/terms/rightsHolder",
  references                    = "http://purl.org/dc/terms/references",
  datasetID                     = "http://rs.tdwg.org/dwc/terms/datasetID",
  datasetName                   = "http://rs.tdwg.org/dwc/terms/datasetName",
  institutionID                 = "http://rs.tdwg.org/dwc/terms/institutionID",
  collectionID                  = "http://rs.tdwg.org/dwc/terms/collectionID",
  institutionCode               = "http://rs.tdwg.org/dwc/terms/institutionCode",
  collectionCode                = "http://rs.tdwg.org/dwc/terms/collectionCode",
  ownerInstitutionCode          = "http://rs.tdwg.org/dwc/terms/ownerInstitutionCode",
  type                          = "http://purl.org/dc/terms/type",
  basisOfRecord                 = "http://rs.tdwg.org/dwc/terms/basisOfRecord",
  informationWithheld           = "http://rs.tdwg.org/dwc/terms/informationWithheld",
  dataGeneralizations           = "http://rs.tdwg.org/dwc/terms/dataGeneralizations",
  catalogNumber                 = "http://rs.tdwg.org/dwc/terms/catalogNumber",
  recordNumber                  = "http://rs.tdwg.org/dwc/terms/recordNumber",
  recordedBy                    = "http://rs.tdwg.org/dwc/terms/recordedBy",
  individualCount               = "http://rs.tdwg.org/dwc/terms/individualCount",
  organismQuantity              = "http://rs.tdwg.org/dwc/terms/organismQuantity",
  organismQuantityType          = "http://rs.tdwg.org/dwc/terms/organismQuantityType",
  occurrenceStatus              = "http://rs.tdwg.org/dwc/terms/occurrenceStatus",
  occurrenceRemarks             = "http://rs.tdwg.org/dwc/terms/occurrenceRemarks",
  preparations                  = "http://rs.tdwg.org/dwc/terms/preparations",
  eventDate                     = "http://rs.tdwg.org/dwc/terms/eventDate",
  year                          = "http://rs.tdwg.org/dwc/terms/year",
  samplingProtocol              = "http://rs.tdwg.org/dwc/terms/samplingProtocol",
  country                       = "http://rs.tdwg.org/dwc/terms/country",
  countryCode                   = "http://rs.tdwg.org/dwc/terms/countryCode",
  stateProvince                 = "http://rs.tdwg.org/dwc/terms/stateProvince",
  locality                      = "http://rs.tdwg.org/dwc/terms/locality",
  minimumDepthInMeters          = "http://rs.tdwg.org/dwc/terms/minimumDepthInMeters",
  maximumDepthInMeters          = "http://rs.tdwg.org/dwc/terms/maximumDepthInMeters",
  decimalLatitude               = "http://rs.tdwg.org/dwc/terms/decimalLatitude",
  decimalLongitude              = "http://rs.tdwg.org/dwc/terms/decimalLongitude",
  geodeticDatum                 = "http://rs.tdwg.org/dwc/terms/geodeticDatum",
  coordinateUncertaintyInMeters = "http://rs.tdwg.org/dwc/terms/coordinateUncertaintyInMeters",
  georeferenceSources           = "http://rs.tdwg.org/dwc/terms/georeferenceSources",
  georeferenceProtocol          = "http://rs.tdwg.org/dwc/terms/georeferenceProtocol",
  georeferenceRemarks           = "http://rs.tdwg.org/dwc/terms/georeferenceRemarks",
  identifiedBy                  = "http://rs.tdwg.org/dwc/terms/identifiedBy",
  dateIdentified                = "http://rs.tdwg.org/dwc/terms/dateIdentified",
  typeStatus                    = "http://rs.tdwg.org/dwc/terms/typeStatus",
  scientificName                = "http://rs.tdwg.org/dwc/terms/scientificName",
  scientificNameAuthorship      = "http://rs.tdwg.org/dwc/terms/scientificNameAuthorship",
  taxonID                       = "http://rs.tdwg.org/dwc/terms/taxonID",
  scientificNameID              = "http://rs.tdwg.org/dwc/terms/scientificNameID",
  acceptedNameUsageID           = "http://rs.tdwg.org/dwc/terms/acceptedNameUsageID",
  taxonRank                     = "http://rs.tdwg.org/dwc/terms/taxonRank",
  taxonomicStatus               = "http://rs.tdwg.org/dwc/terms/taxonomicStatus",
  nomenclaturalCode             = "http://rs.tdwg.org/dwc/terms/nomenclaturalCode",
  kingdom                       = "http://rs.tdwg.org/dwc/terms/kingdom",
  phylum                        = "http://rs.tdwg.org/dwc/terms/phylum",
  class                         = "http://rs.tdwg.org/dwc/terms/class",
  order                         = "http://rs.tdwg.org/dwc/terms/order",
  family                        = "http://rs.tdwg.org/dwc/terms/family",
  genus                         = "http://rs.tdwg.org/dwc/terms/genus",
  subgenus                      = "http://rs.tdwg.org/dwc/terms/subgenus",
  specificEpithet               = "http://rs.tdwg.org/dwc/terms/specificEpithet",
  infraspecificEpithet          = "http://rs.tdwg.org/dwc/terms/infraspecificEpithet",
  fieldNumber                   = "http://rs.tdwg.org/dwc/terms/fieldNumber",
  associatedReferences          = "http://rs.tdwg.org/dwc/terms/associatedReferences",
  dynamicProperties             = "http://rs.tdwg.org/dwc/terms/dynamicProperties"
)

# Assert every column in occ has a term mapping (catches columns silently added later).
missing_terms <- setdiff(names(occ), names(dwc_terms))
if (length(missing_terms) > 0) {
  stop("meta.xml term mapping missing for columns: ",
       paste(missing_terms, collapse = ", "))
}

field_lines <- map_chr(seq_along(names(occ)) - 1, function(i) {
  col <- names(occ)[i + 1]
  sprintf('    <field index="%d" term="%s"/>', i, dwc_terms[[col]])
})

meta_xml <- c(
  '<?xml version="1.0" encoding="UTF-8"?>',
  '<archive xmlns="http://rs.tdwg.org/dwc/text/" metadata="eml.xml">',
  '  <core encoding="UTF-8"',
  '        fieldsTerminatedBy="\\t"',
  '        linesTerminatedBy="\\n"',
  '        fieldsEnclosedBy=""',
  '        ignoreHeaderLines="1"',
  '        rowType="http://rs.tdwg.org/dwc/terms/Occurrence">',
  '    <files>',
  '      <location>occurrence.txt</location>',
  '    </files>',
  '    <id index="0"/>',
  field_lines,
  '  </core>',
  '  <extension encoding="UTF-8"',
  '             fieldsTerminatedBy="\\t"',
  '             linesTerminatedBy="\\n"',
  '             fieldsEnclosedBy=""',
  '             ignoreHeaderLines="1"',
  '             rowType="http://rs.tdwg.org/dwc/terms/MeasurementOrFact">',
  '    <files>',
  '      <location>measurementorfact.txt</location>',
  '    </files>',
  '    <coreid index="0"/>',
  '    <field index="1" term="http://rs.tdwg.org/dwc/terms/measurementType"/>',
  '    <field index="2" term="http://rs.tdwg.org/dwc/terms/measurementValue"/>',
  '    <field index="3" term="http://rs.tdwg.org/dwc/terms/measurementUnit"/>',
  '  </extension>',
  '</archive>'
)
writeLines(meta_xml, file.path(dwc_dir, "meta.xml"), sep = "\n")
cat("Wrote", file.path(dwc_dir, "meta.xml"), "\n")

# -----------------------------------------------------------------------------
# 8.8  Build eml.xml (full GBIF EML 1.2 profile, no private comments)
# -----------------------------------------------------------------------------
cat("Building eml.xml...\n")

eml_xml <- c(
  '<?xml version="1.0" encoding="UTF-8"?>',
  '<eml:eml xmlns:eml="eml://ecoinformatics.org/eml-2.1.1"',
  '         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"',
  '         xsi:schemaLocation="eml://ecoinformatics.org/eml-2.1.1 http://rs.gbif.org/schema/eml-gbif-profile/1.2/eml.xsd"',
  '         packageId="JCU-ECHINODERM-DATASET/v1"',
  '         system="http://gbif.org" scope="system" xml:lang="eng">',
  '  <dataset>',
  '    <title xml:lang="eng">Integrated echinoderm occurrence dataset for northeast Australia and the adjacent Coral Sea</title>',
  '    <creator>',
  '      <individualName><givenName>Maria Camila</givenName><surName>Velez Diaz</surName></individualName>',
  '      <organizationName>James Cook University</organizationName>',
  '      <positionName>Postgraduate Researcher</positionName>',
  '      <electronicMailAddress>mariacamila.velezdiaz@my.jcu.edu.au</electronicMailAddress>',
  '      <userId directory="https://orcid.org/">0000-0003-4180-1077</userId>',
  '    </creator>',
  '    <metadataProvider>',
  '      <individualName><givenName>Maria Camila</givenName><surName>Velez Diaz</surName></individualName>',
  '      <organizationName>James Cook University</organizationName>',
  '      <electronicMailAddress>mariacamila.velezdiaz@my.jcu.edu.au</electronicMailAddress>',
  '    </metadataProvider>',
  '    <associatedParty>',
  '      <individualName><givenName>Alastair</givenName><surName>Birtles</surName></individualName>',
  '      <organizationName>James Cook University</organizationName>',
  #  <<< OPTIONAL: add <userId directory="https://orcid.org/">ORCID</userId> >>>
  '      <role>coPrincipalInvestigator</role>',
  '    </associatedParty>',
  '    <associatedParty>',
  '      <individualName><givenName>Sue-Ann</givenName><surName>Watson</surName></individualName>',
  '      <organizationName>James Cook University</organizationName>',
  #  <<< OPTIONAL: add <userId directory="https://orcid.org/">ORCID</userId> >>>
  '      <role>coPrincipalInvestigator</role>',
  '    </associatedParty>',
  '    <associatedParty>',
  '      <organizationName>Australian Museum</organizationName>',
  '      <role>contentProvider</role>',
  '    </associatedParty>',
  '    <associatedParty>',
  '      <organizationName>Queensland Museum Tropics</organizationName>',
  '      <role>contentProvider</role>',
  '    </associatedParty>',
  paste0('    <pubDate>', format(Sys.Date(), "%Y-%m-%d"), '</pubDate>'),
  '    <language>eng</language>',
  '    <abstract>',
  '      <para>An integrated occurrence dataset of echinoderms (Phylum Echinodermata) for northeast Australia and the adjacent Coral Sea, spanning shelf, slope and abyssal habitats between 10 degrees S to 29 degrees S and 142 degrees E to 154 degrees E. The dataset consolidates 43,470 curated occurrence records assembled from 13 source datasets covering national and international biodiversity aggregators, a national research infrastructure, and direct Collection Management System exports from the Queensland Museum Tropics and the Australian Museum. Records were integrated through a reproducible R workflow that harmonises Darwin Core fields, deduplicates records across sources, and resolves conflicts in coordinates, collection year, basis of record, depth and taxonomic name against the World Register of Marine Species (WoRMS). Per-record data-quality flags (coordinate quality, depth quality, taxonomic resolution) are provided as MeasurementOrFact rows and in dynamicProperties. Full methods are described in the accompanying Scientific Data descriptor.</para>',
  '    </abstract>',
  '    <keywordSet>',
  '      <keyword>Occurrence</keyword>',
  '      <keyword>Specimen</keyword>',
  '      <keywordThesaurus>GBIF Dataset Type Vocabulary</keywordThesaurus>',
  '    </keywordSet>',
  '    <keywordSet>',
  '      <keyword>Echinodermata</keyword>',
  '      <keyword>Asteroidea</keyword>',
  '      <keyword>Ophiuroidea</keyword>',
  '      <keyword>Echinoidea</keyword>',
  '      <keyword>Holothuroidea</keyword>',
  '      <keyword>Crinoidea</keyword>',
  '      <keyword>Coral Sea</keyword>',
  '      <keyword>Great Barrier Reef</keyword>',
  '      <keyword>Queensland</keyword>',
  '      <keyword>deep sea</keyword>',
  '    </keywordSet>',
  '    <intellectualRights>',
  '      <para>This work is licensed under a <ulink url="https://creativecommons.org/licenses/by/4.0/legalcode"><citetitle>Creative Commons Attribution (CC BY 4.0) License</citetitle></ulink>. Includes records supplied by the Australian Museum and by Queensland Museum Tropics, both used with written permission and acknowledged as contributing institutions.</para>',
  '    </intellectualRights>',
  '    <distribution scope="document">',
  '      <online>',
  # <<< FILL IN: IPT / Zenodo DwC-A URL once published >>>
  '        <url function="download">[Insert IPT or Zenodo DwC-A URL after publication]</url>',
  '      </online>',
  '    </distribution>',
  '    <coverage>',
  '      <geographicCoverage>',
  '        <geographicDescription>Northeast Australia and the adjacent Coral Sea, spanning shelf, slope and abyssal habitats from the Torres Strait to southeast Queensland, including the Gulf of Carpentaria and adjacent Coral Sea waters.</geographicDescription>',
  '        <boundingCoordinates>',
  '          <westBoundingCoordinate>142</westBoundingCoordinate>',
  '          <eastBoundingCoordinate>154</eastBoundingCoordinate>',
  '          <northBoundingCoordinate>-10</northBoundingCoordinate>',
  '          <southBoundingCoordinate>-29</southBoundingCoordinate>',
  '        </boundingCoordinates>',
  '      </geographicCoverage>',
  '      <temporalCoverage>',
  '        <rangeOfDates>',
  '          <beginDate><calendarDate>1861-01-01</calendarDate></beginDate>',
  '          <endDate><calendarDate>2026-07-30</calendarDate></endDate>',
  '        </rangeOfDates>',
  '      </temporalCoverage>',
  '      <taxonomicCoverage>',
  '        <taxonomicClassification><taxonRankName>phylum</taxonRankName><taxonRankValue>Echinodermata</taxonRankValue></taxonomicClassification>',
  '        <taxonomicClassification><taxonRankName>class</taxonRankName><taxonRankValue>Asteroidea</taxonRankValue></taxonomicClassification>',
  '        <taxonomicClassification><taxonRankName>class</taxonRankName><taxonRankValue>Ophiuroidea</taxonRankValue></taxonomicClassification>',
  '        <taxonomicClassification><taxonRankName>class</taxonRankName><taxonRankValue>Echinoidea</taxonRankValue></taxonomicClassification>',
  '        <taxonomicClassification><taxonRankName>class</taxonRankName><taxonRankValue>Holothuroidea</taxonRankValue></taxonomicClassification>',
  '        <taxonomicClassification><taxonRankName>class</taxonRankName><taxonRankValue>Crinoidea</taxonRankValue></taxonomicClassification>',
  '      </taxonomicCoverage>',
  '    </coverage>',
  '    <maintenance>',
  '      <description><para>Static snapshot corresponding to source download dates listed in the associated Scientific Data descriptor.</para></description>',
  '      <maintenanceUpdateFrequency>notPlanned</maintenanceUpdateFrequency>',
  '    </maintenance>',
  '    <contact>',
  '      <individualName><givenName>Maria Camila</givenName><surName>Velez Diaz</surName></individualName>',
  '      <organizationName>James Cook University</organizationName>',
  '      <electronicMailAddress>mariacamila.velezdiaz@my.jcu.edu.au</electronicMailAddress>',
  '      <userId directory="https://orcid.org/">0000-0003-4180-1077</userId>',
  '    </contact>',
  '    <methods>',
  '      <methodStep>',
  '        <description>',
  '          <para>Occurrence records were integrated from 13 source datasets, harmonised to Darwin Core, deduplicated using a hierarchical occurrenceID / catalogNumber / composite-key strategy, and cross-source conflicts in coordinates, collection year, depth, basis of record and taxonomic name were resolved through a scripted R workflow. All scientific names were validated against the World Register of Marine Species (WoRMS). Per-record data-quality flags are provided in the MeasurementOrFact extension and in dynamicProperties. Full methods are described in the associated Scientific Data descriptor.</para>',
  '        </description>',
  '      </methodStep>',
  '    </methods>',
  '    <project>',
  '      <title>Integrated echinoderm occurrence dataset for northeast Australia and the adjacent Coral Sea</title>',
  '      <personnel>',
  '        <individualName><givenName>Maria Camila</givenName><surName>Velez Diaz</surName></individualName>',
  '        <role>principalInvestigator</role>',
  '      </personnel>',
  # <<< OPTIONAL: add <funding> block if grant information applies >>>
  '    </project>',
  '  </dataset>',
  '  <additionalMetadata>',
  '    <metadata>',
  '      <gbif>',
  paste0('        <dateStamp>', format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), '</dateStamp>'),
  '        <hierarchyLevel>dataset</hierarchyLevel>',
  # <<< FILL IN: dataset DOI once minted >>>
  '        <citation>Velez Diaz, M.C., Birtles, A. and Watson, S.-A. (2026). Integrated echinoderm occurrence dataset for northeast Australia and the adjacent Coral Sea. James Cook University. https://doi.org/[Dataset-DOI]</citation>',
  '        <resourceLogoUrl></resourceLogoUrl>',
  '      </gbif>',
  '    </metadata>',
  '  </additionalMetadata>',
  '</eml:eml>'
)
writeLines(eml_xml, file.path(dwc_dir, "eml.xml"), sep = "\n")
cat("Wrote", file.path(dwc_dir, "eml.xml"), "\n")

# -----------------------------------------------------------------------------
# 8.9  Package ZIP
# -----------------------------------------------------------------------------
zip_target <- "echinoderm_dwca_public_release.zip"
if (file.exists(zip_target)) file.remove(zip_target)

old_wd <- getwd()
setwd(dwc_dir)
utils::zip(zipfile = file.path("..", zip_target),
           files = c("occurrence.txt", "measurementorfact.txt",
                     "meta.xml", "eml.xml"))
setwd(old_wd)

cat("\n=== Section 8 complete ===\n")
cat("Written:\n")
cat("  ", file.path(dwc_dir, "occurrence.txt"), "\n")
cat("  ", file.path(dwc_dir, "measurementorfact.txt"), "\n")
cat("  ", file.path(dwc_dir, "meta.xml"), "\n")
cat("  ", file.path(dwc_dir, "eml.xml"), "\n")
cat("  ", zip_target, "\n\n")

cat("Next steps before public deposit:\n")
cat("  1. Fill in the five placeholders marked '<<< FILL IN >>>' in this script:\n")
cat("       - Scientific Data DOI (references)\n")
cat("       - Dataset DOI (datasetID, EML citation)\n")
cat("       - IPT / Zenodo download URL (EML <distribution>)\n")
cat("       - Alastair Birtles / Sue-Ann Watson ORCID iDs (EML <associatedParty>)\n")
cat("       - Grant / funding information (EML <project>)\n")
cat("  2. Validate the archive at https://tools.gbif.org/dwca-validator/\n")
cat("  3. Deposit on Zenodo for an immediate citable DOI.\n")
cat("  4. Submit through the OBIS Australia IPT (or ALA IPT);\n")
cat("     federation to GBIF is automatic once accepted.\n")
