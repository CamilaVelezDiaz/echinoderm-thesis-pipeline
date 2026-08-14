# =============================================================================
# ECHINODERM DATASET INTEGRATION — TWO-LAYER ARCHITECTURE
# Camila Velez | JCU Marine Biology MSc | Phase 1
#
# LAYER A: echino_long.csv
#   - ALL records from ALL sources, no deduplication, ALL original columns
#
# LAYER B: echino_wide.csv
#   - One row per unique specimen (record_key)
#   - For fields with conflicts between sources: separate columns per source
#     e.g. depth_ALA_GBRSBD, depth_CSIRO_GBRSBD, depth_OBIS
#   - For fields with no conflict: single column with the value
#   - echino_wide_depth.csv: same but filtered to records with at least one depth value
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================

setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")

# =============================================================================
# Library
# =============================================================================

library(tidyverse)
library(vegan)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(stringr)
library(readxl)
library(worrms)

# =============================================================================
# SECTION 1: LOAD ALL SOURCES
# =============================================================================

# --- GBIF ---
gbif_other <- read_delim(
  "C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\GBIF Queensland Museum\\occurrence.txt",
  delim = "\t", show_col_types = FALSE
)
gbif_echi <- read_delim(
  "C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\GIBF - Echinodermata data_Au\\occurrence.txt",
  delim = "\t", show_col_types = FALSE
) %>%
  mutate(
    minimumDepthInMeters = as.numeric(depth),
    maximumDepthInMeters = as.numeric(depth)
  )

# --- ALA ---
ala_echi   <- read_csv("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\ALA - Echinodermata dat\\records-2026-04-10.csv",
                       show_col_types = FALSE)
ala_other  <- read_tsv("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\records Atlas - Other invertebrates collection\\records-2026-03-13.tsv",
                       show_col_types = FALSE)
ala_GBRSBD <- read_tsv("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\Record Atlas GBRBD\\records-2026-03-13.tsv",
                       show_col_types = FALSE)

# --- CSIRO ---
csiro_gbrsbd <- read_delim("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\Csiro- GBRSBD\\occurrence.txt",
                           delim = "\t", show_col_types = FALSE)
csiro_other  <- read_delim("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\Csiro-QLDM other invertebrates - marine records\\occurrence.txt",
                           delim = "\t", show_col_types = FALSE)

# --- IDigBio ---
idigbio_raw <- read_csv("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\iDigBio data\\occurrence.csv",
                        show_col_types = FALSE)
idigbio <- idigbio_raw %>%
  rename_with(~ str_remove(.x, "^dwc:"), starts_with("dwc:")) %>%
  mutate(
    decimalLatitude  = as.numeric(str_extract(`idigbio:geoPoint`, "^-?[0-9.]+")),
    decimalLongitude = as.numeric(str_extract(`idigbio:geoPoint`, "-?[0-9.]+$")),
    occurrenceID     = if_else(is.na(occurrenceID) | occurrenceID == "",
                               coreid, occurrenceID),
    year             = as.integer(str_extract(eventDate, "^[0-9]{4}"))
  )

# --- OZCAM ---
ozcam <- read_csv("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\OZCAM ALA data\\records-2026-04-03.csv",
                  show_col_types = FALSE)

# --- OBIS snapshot — downloaded 2026-05-25 ---
obis_echi <- read_csv("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\OBIS\\obis_echi_snapshot_2026-05-25.csv",
                      show_col_types = FALSE) %>%
  mutate(
    minimumDepthInMeters = if_else(
      is.na(as.numeric(minimumDepthInMeters)) & !is.na(depth),
      as.numeric(depth), as.numeric(minimumDepthInMeters)
    ),
    maximumDepthInMeters = if_else(
      is.na(as.numeric(maximumDepthInMeters)) & !is.na(depth),
      as.numeric(depth), as.numeric(maximumDepthInMeters)
    )
  )

# --- MICHELA CIDARIS (MTQ_CIDARIS) ---
michela_raw <- read_excel(
  "C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\Cidaris - Michela\\MTQ echinoderms - registered - June 2021 final.xlsx",
  sheet = "Raw data from Vernon"
)

# HELPER: Convert DMS to decimal degrees
dms_to_dd <- function(x) {
  x <- as.character(x)
  
  deg <- as.numeric(str_extract(x, "^\\d+"))
  min <- as.numeric(str_extract(x, "(?<=°\\s?)\\d+"))
  sec <- as.numeric(str_extract(x, "(?<=')\\s?[\\d.]+"))
  dir <- str_extract(x, "(North|South|East|West)")
  
  dd <- deg + min / 60 + replace_na(sec, 0) / 3600
  
  if_else(
    is.na(dir),
    -dd,                                        
    if_else(dir %in% c("South", "West"), -dd, dd)
  )
}

# HELPER: Parse depth strings like "32.3m" to "32.3"
parse_depth <- function(x) {
  suppressWarnings(as.numeric(str_remove_all(as.character(x), "[^0-9.]")))
}

# STANDARDISE
michela <- michela_raw %>%
  mutate(
    source        = "MTQ_CIDARIS",
    catalogNumber = as.character(`Registration Number`),
    occurrenceID  = NA_character_,
    scientificName = str_trim(`Taxonomic Classification`),
    phylum         = "Echinodermata",
    decimalLatitude  = dms_to_dd(`Field Coll Latitude from`),
    decimalLongitude = dms_to_dd(`Field Coll Longitude from`),
    minimumDepthInMeters = parse_depth(`Field Coll Depth from`),
    maximumDepthInMeters = parse_depth(`Field Coll Depth to`),
    eventDate = as.character(as.Date(as.numeric(`Field Coll Date`), 
                                     origin = "1899-12-30")),
    year      = as.integer(format(as.Date(as.numeric(`Field Coll Date`),
                                          origin = "1899-12-30"), "%Y")),
    locality         = coalesce(as.character(`Locality Name`),
                                as.character(`Locality Remarks`)),
    stateProvince    = "Queensland",
    country          = "Australia",
    recordedBy      = as.character(Identifier),
    basisOfRecord   = "PreservedSpecimen",
    institutionCode = "MTQ",
    collectionCode  = as.character(`Collecting Unit Name`),
    datasetName     = "MTQ Echinoderms registered June 2021",
    preparations    = as.character(`Field Coll Specimen Category`),
    verbatimLocality = as.character(`Field Coll Place`),
    fieldNumber     = as.character(`Field Coll Ref`)
  )

# --- STEFANO CIDARIS (CIDARIS_QMT) --- 
stefano_raw <- read_excel(
  "C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\MQT - Stefano\\CIDARIS_Echinoderms.xlsx"
)

stefano_dates <- dmy(stefano_raw$eventDate, locale = 'English')

stefano <- stefano_raw %>%
  mutate(
    source           = "CIDARIS_QMT",
    catalogNumber    = as.character(`Registration Number`),
    occurrenceID     = NA_character_,
    scientificName   = str_trim(`Taxonomic Classification`),
    phylum           = "Echinodermata",
    decimalLatitude  = dms_to_dd(`Field Coll Latitude from`),
    decimalLongitude = dms_to_dd(`Field Coll Longitude from`),
    minimumDepthInMeters = as.numeric(`Field Coll Depth from`),
    maximumDepthInMeters = as.numeric(`Field Coll Depth to`),
    eventDate        = as.character(stefano_dates),
    year             = as.integer(format(stefano_dates, "%Y")),
    locality         = str_trim(str_extract(
      `Object  (short  summary)`,
      "(?<=;\\s)(Coral Sea[^;]+|[^;]+(?=;\\s*M Pichon))"
    )),
    stateProvince    = "Queensland",
    country          = "Australia",
    recordedBy       = "M Pichon, A Birtles, P Arnold",
    basisOfRecord    = "PreservedSpecimen",
    institutionCode  = "QMT",
    collectionCode   = as.character(`Collecting Unit`),
    datasetName      = "CIDARIS Echinoderms QMT",
    fieldNumber      = as.character(`Field Coll Ref`),
    order            = as.character(`Taxon - Order`)
  )

# --- AM DIRECT EXPORT (AM_direct) ---
# Direct CMS export from Australian Museum supplied in response to data
# request (26-015_Velez.csv). 5,656 echinoderm specimen records from
# Queensland waters, all PreservedSpecimen, with depth coverage for
# 5,654/5,656 records (0-2,523m).
am_raw <- read_csv(
  "C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1\\AM ECHINO DEPTH\\26-015 Velez.csv",
  locale = locale(encoding = "latin1"),
  show_col_types = FALSE
)

am_direct <- am_raw %>%
  mutate(
    source           = "AM_direct",
    occurrenceID     = NA_character_,
    catalogNumber    = as.character(`Registration Number`),
    scientificName   = str_trim(paste(
      coalesce(as.character(Genus), ""),
      coalesce(as.character(Species), "")
    )),
    scientificName   = if_else(str_trim(scientificName) == "",
                               NA_character_, str_trim(scientificName)),
    phylum           = "Echinodermata",
    class            = as.character(Class),
    family           = as.character(Family),
    genus            = as.character(Genus),
    specificEpithet  = as.character(Species),
    identificationQualifier = as.character(`Identification Qualifier`),
    # Coordinates — midpoint if tow/transect, single point otherwise
    decimalLatitude  = if_else(
      !is.na(`Decimal latitude to`) &
        `Decimal latitude to` != `Decimal latitude from`,
      (`Decimal latitude from` + `Decimal latitude to`) / 2,
      `Decimal latitude from`
    ),
    decimalLongitude = if_else(
      !is.na(`Decimal longitude to`) &
        `Decimal longitude to` != `Decimal longitude from`,
      (`Decimal longitude from` + `Decimal longitude to`) / 2,
      `Decimal longitude from`
    ),
    minimumDepthInMeters = as.numeric(`Bottom Depth From`),
    maximumDepthInMeters = as.numeric(`Bottom Depth To`),
    # Date: parse "08 Feb 1987" to ISO
    eventDate        = format(
      suppressWarnings(
        as.Date(`Date Visited From`, format = "%d %b %Y")
      ), "%Y-%m-%d"
    ),
    year             = as.integer(format(
      suppressWarnings(
        as.Date(`Date Visited From`, format = "%d %b %Y")
      ), "%Y"
    )),
    country          = as.character(Country),
    stateProvince    = as.character(State),
    locality         = as.character(`Precise Location`),
    verbatimLocality = as.character(`Precise Location`),
    recordedBy       = as.character(`Collection Participant`),
    basisOfRecord    = "PreservedSpecimen",
    institutionCode  = "AM",
    collectionCode   = "INVERTEBRATES - MARINE & OTHER",
    datasetName      = "AM Direct Export 26-015_Velez",
    preparations     = as.character(`Kind of Object`)
  )

cat("all sources loaded.\n")


# =============================================================================
# SECTION 2: SOURCE PRIORITY ORDER
# Lower number = more authoritative / more QLD-specific
# =============================================================================

source_priority_levels <- c(
  "ALA_GBRSBD",         # 1 — most targeted: GBRSBD via ALA
  "ALA_other",          # 2 — QM other invertebrates via ALA
  "ALA_echinodermata",  # 3 — broad echinoderm pull via ALA
  "CSIRO_GBRSBD",       # 4 — GBRSBD via CSIRO
  "CSIRO_QM",           # 5 — QM other invertebrates via CSIRO
  "MTQ_CIDARIS",        # 6 — direct MTQ CMS export
  "CIDARIS_QMT",        # 7 — direct QMT CMS export (CIDARIS expeditions)
  "AM_direct",          # 8 — direct AM CMS export, high authority
  "GBIF_QM",            # 9 — QM other invertebrates via GBIF
  "GBIF_echinodermata", # 10 — broad echinoderm pull via GBIF
  "OZCAM",              # 11 — Australian collections aggregator
  "iDigBio",            # 12 — international aggregator
  "OBIS"                # 13 — OBIS (good depth, sometimes less taxonomic detail)
)


# =============================================================================
# SECTION 3: STANDARDISE — ADD SOURCE LABEL + SAFE TYPES
#             ALL original columns are preserved
# =============================================================================

standardise_source <- function(df, source_name) {
  
  # Numeric columns to preserve as numeric
  num_cols <- c(
    "decimalLatitude", "decimalLongitude",
    "minimumDepthInMeters", "maximumDepthInMeters",
    "minimumElevationInMeters", "maximumElevationInMeters",
    "coordinateUncertaintyInMeters", "coordinatePrecision",
    "individualCount", "organismQuantity",
    "depth", "elevation"
  )
  int_cols <- c("year", "month", "day", "startDayOfYear", "endDayOfYear")
  
  df %>%
    # Step 1: coerce ALL date/datetime columns to character (by type, not name)
    mutate(across(
      where(~ inherits(.x, "POSIXct") | inherits(.x, "POSIXt") | inherits(.x, "Date")),
      as.character
    )) %>%
    # Step 2: coerce logical, factor, and any double that isn't a known
    # numeric column to character — resolves type conflicts across sources
    mutate(across(
      !any_of(c(num_cols, int_cols)) &
        where(~ is.logical(.x) | is.factor(.x) | is.double(.x)),
      as.character
    )) %>%
    # Step 3: enforce correct types and add source label
    mutate(
      source = source_name,
      across(any_of(num_cols), ~ suppressWarnings(as.numeric(.x))),
      across(any_of(int_cols), ~ suppressWarnings(as.integer(.x)))
    )
}

# Bind all sources — columns missing from a source become NA
echino_long_raw <- bind_rows(
  standardise_source(ala_GBRSBD,   "ALA_GBRSBD"),
  standardise_source(ala_other,    "ALA_other"),
  standardise_source(ala_echi,     "ALA_echinodermata"),
  standardise_source(csiro_gbrsbd, "CSIRO_GBRSBD"),
  standardise_source(csiro_other,  "CSIRO_QM"),
  standardise_source(gbif_other,   "GBIF_QM"),
  standardise_source(gbif_echi,    "GBIF_echinodermata"),
  standardise_source(idigbio,      "iDigBio"),
  standardise_source(ozcam,        "OZCAM"),
  standardise_source(obis_echi,    "OBIS"),
  standardise_source(michela,      "MTQ_CIDARIS"),
  standardise_source(stefano,      "CIDARIS_QMT"),
  standardise_source(am_direct,    "AM_direct")   # NEW
) %>%
  # Fix any remaining logical columns (all-NA from missing sources)
  mutate(across(where(is.logical), as.character)) %>%
  mutate(across(where(~ inherits(.x, "POSIXct") | inherits(.x, "POSIXt")), as.character))

cat("✅ All sources combined. Total rows:", nrow(echino_long_raw), "\n")
cat("   Total columns:", ncol(echino_long_raw), "\n")


# =============================================================================
# SECTION 4: FILTER — ECHINODERMATA
# =============================================================================

echino_long_raw <- echino_long_raw %>%
  filter(
    str_to_lower(phylum) == "echinodermata" |
      class %in% c("Asteroidea", "Echinoidea", "Holothuroidea",
                   "Ophiuroidea", "Crinoidea") |
      str_detect(str_to_lower(coalesce(scientificName, "")), "echinoderm")
  )

cat("✅ After Echinodermata filter:", nrow(echino_long_raw), "rows\n")


# =============================================================================
# SECTION 5: FILTER — Northeast Australia
# =============================================================================
# Primary filter: bounding box (10°S–29°S, 142°E–154°E)
# Supplementary: stateProvince text match to capture Gulf of Carpentaria
# (west of 142°E) and Torres Strait (north of 10°S)
# Note: PNG records removed post-hoc in postprocessing script Section 1b

echino_long_raw <- echino_long_raw %>%
  mutate(
    in_qld_bbox = !is.na(decimalLatitude) & !is.na(decimalLongitude) &
      decimalLatitude  >= -29 & decimalLatitude  <= -10 &
      decimalLongitude >= 142 & decimalLongitude <= 154,
    in_qld_text = str_detect(str_to_lower(coalesce(stateProvince, "")),
                             "queensland|qld")
  ) %>%
  filter(in_qld_bbox | in_qld_text) %>%
  filter(coalesce(catalogNumber, "") != "G292")   # confirmed phantom crustacean record

cat("✅ After Northeast Australia filter:", nrow(echino_long_raw), "rows\n")


# =============================================================================
# SECTION 6: ASSIGN RECORD KEY
# =============================================================================
# Links the same physical specimen across sources.
# Hierarchy: occurrenceID > catalogNumber > composite fallback
#
# IMPORTANT: occurrenceID is normalised to uppercase before key assignment.
# iDigBio stores UUIDs and MCZ identifiers in lowercase (e.g. "0042b589-..."
# or "mcz:iz:ast-2321") while ALA/GBIF use uppercase for the same identifiers.
# Uppercasing prevents the same specimen from receiving two different record_keys.

keys_before_reconciliation <- n_distinct(echino_long_raw$record_key)

echino_long_raw <- echino_long_raw %>%
  mutate(
    occurrenceID = str_to_upper(occurrenceID),
    record_key = case_when(
      !is.na(occurrenceID) & occurrenceID != ""   ~ occurrenceID,
      !is.na(catalogNumber) & catalogNumber != "" ~ paste0("CAT_", catalogNumber),
      TRUE ~ paste(
        str_to_lower(str_trim(coalesce(scientificName, "unknown"))),
        round(coalesce(decimalLatitude,  -999), 4),
        round(coalesce(decimalLongitude, -999), 4),
        str_sub(coalesce(eventDate, "unknown"), 1, 10),
        sep = "|"
      )
    ),
    source_priority = coalesce(match(source, source_priority_levels), 99L)
  )

keys_before_reconciliation <- n_distinct(echino_long_raw$record_key)
cat("✅ Unique record keys:", keys_before_reconciliation, "\n")


# =============================================================================
# SECTION 6b: SECONDARY RECONCILIATION BY CATALOG NUMBER
# =============================================================================
# Problem: the same physical specimen can arrive from ALA, GBIF, and OBIS
# with three different occurrenceIDs — Section 6 assigns three record_keys.
# Fix: find all record_keys sharing the same catalogNumber and reassign to
# a single canonical key from the highest-priority source.
#
# IMPORTANT: catalogNumber is matched case-insensitively (uppercased before
# grouping) — iDigBio stores catalogNumber in lowercase (e.g. "f241951")
# while ALA/GBIF/OZCAM use uppercase ("F241951") for the same specimen.
# Without this fix, 18,387 specimens were retained as duplicate rows
# under different record_keys purely due to case mismatch.

catnum_lookup <- echino_long_raw %>%
  filter(!is.na(catalogNumber) & catalogNumber != "") %>%
  mutate(catalogNumber_upper = str_to_upper(catalogNumber)) %>%
  arrange(catalogNumber_upper, source_priority) %>%
  group_by(catalogNumber_upper) %>%
  slice(1) %>%
  ungroup() %>%
  select(catalogNumber_upper, canonical_key = record_key)

echino_long_raw <- echino_long_raw %>%
  mutate(catalogNumber_upper = str_to_upper(catalogNumber)) %>%
  left_join(catnum_lookup, by = "catalogNumber_upper") %>%
  mutate(
    record_key = if_else(!is.na(canonical_key), canonical_key, record_key)
  ) %>%
  select(-canonical_key, -catalogNumber_upper)

cat("✅ Unique record keys (after catalogNumber reconciliation):",
    n_distinct(echino_long_raw$record_key), "\n")
cat("   Rows collapsed by catalogNumber reconciliation:",
    keys_before_reconciliation - n_distinct(echino_long_raw$record_key), "\n")


# =============================================================================
# SECTION 7: LAYER A — SAVE LONG FORMAT
# All records, all columns, no deduplication. Permanent audit trail.
# =============================================================================

echino_long <- echino_long_raw %>%
  relocate(record_key, source, source_priority)

write_csv(echino_long, "echino_long.csv")
cat("\n✅ LAYER A saved: echino_long.csv\n")
cat("   Rows:", nrow(echino_long), "| Columns:", ncol(echino_long), "\n")


# =============================================================================
# SECTION 8: LAYER B — WIDE FORMAT (one row per record_key)
#
# conflict_fields: pivoted per source so values can be compared side by side
# consensus_fields: best available value from highest-priority source
# =============================================================================

# --- 8.1 Define conflict fields ---
conflict_fields <- c(
  "minimumDepthInMeters", "maximumDepthInMeters",
  "decimalLatitude", "decimalLongitude",
  "scientificName", "taxonRank",
  "kingdom", "phylum", "class", "order", "family", "genus", "specificEpithet",
  "year", "eventDate", "basisOfRecord",
  "locality", "typeStatus", "identifiedBy", "dateIdentified"
)
conflict_fields <- intersect(conflict_fields, names(echino_long_raw))

# --- 8.2 Per-source columns for conflict fields ---
wide_conflicts <- echino_long_raw %>%
  select(record_key, source, all_of(conflict_fields)) %>%
  pivot_longer(
    cols             = all_of(conflict_fields),
    names_to         = "field",
    values_to        = "value",
    values_transform = list(value = as.character)
  ) %>%
  mutate(col_name = paste0(field, "__", source)) %>%
  select(record_key, col_name, value) %>%
  pivot_wider(
    names_from  = col_name,
    values_from = value,
    values_fn   = ~ paste(unique(na.omit(.x)), collapse = ";")
  )

cat("✅ Per-source conflict columns:", ncol(wide_conflicts) - 1, "\n")

# --- 8.3 Consensus columns for non-conflict fields ---
consensus_fields <- setdiff(
  names(echino_long_raw),
  c(conflict_fields, "source", "source_priority",
    "in_qld_bbox", "in_qld_text", "record_key")
)

wide_consensus <- echino_long_raw %>%
  arrange(record_key, source_priority) %>%
  group_by(record_key) %>%
  summarise(
    across(
      all_of(consensus_fields),
      ~ { vals <- .x[!is.na(.x) & .x != ""]; if (length(vals) == 0) NA else vals[1] }
    ),
    sources_all    = paste(sort(unique(source)), collapse = ";"),
    primary_source = source[1],
    n_sources      = n_distinct(source),
    .groups = "drop"
  )

cat("✅ Consensus columns:", ncol(wide_consensus) - 4, "\n")

# --- 8.4 Join conflict + consensus ---
echino_wide <- wide_consensus %>%
  left_join(wide_conflicts, by = "record_key")

cat("✅ Wide table: ", nrow(echino_wide), "rows |", ncol(echino_wide), "columns\n")


# =============================================================================
# SECTION 9: CONFLICT FLAGS (NAME + COORDINATE ONLY) + COMPLETENESS SUMMARY
# =============================================================================
# NOTE: Depth resolution (min/max, conflict flagging, depth zones) is now
# handled entirely in postprocessing (Section 5), since it involves
# analytical decisions that should be transparent and documented there
# rather than baked into the integration pipeline.

lat_cols        <- names(echino_wide)[str_detect(names(echino_wide), "^decimalLatitude__")]
lon_cols        <- names(echino_wide)[str_detect(names(echino_wide), "^decimalLongitude__")]
sciname_cols    <- names(echino_wide)[str_detect(names(echino_wide), "^scientificName__")]
year_cols_tmp   <- names(echino_wide)[str_detect(names(echino_wide), "^year__")]
family_cols_tmp <- names(echino_wide)[str_detect(names(echino_wide), "^family__")]

lat_mat     <- echino_wide %>% select(all_of(lat_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
lon_mat     <- echino_wide %>% select(all_of(lon_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
sciname_mat <- echino_wide %>% select(all_of(sciname_cols))

n_distinct_vals <- function(mat) {
  apply(mat, 1, function(row) {
    vals <- unique(na.omit(row))
    length(vals) > 1
  })
}

coord_conflict <- n_distinct_vals(as.matrix(lat_mat))
name_conflict  <- n_distinct_vals(as.matrix(sciname_mat))

has_lat     <- apply(as.matrix(lat_mat),     1, function(r) any(!is.na(r)))
has_lon     <- apply(as.matrix(lon_mat),     1, function(r) any(!is.na(r)))
has_sciname <- apply(as.matrix(sciname_mat), 1, function(r) any(!is.na(r)))

echino_wide <- echino_wide %>%
  mutate(
    coord_conflict  = coord_conflict,
    name_conflict   = name_conflict,
    
    completeness_score = (
      as.integer(has_lat) +
        as.integer(has_lon) +
        as.integer(has_sciname) +
        as.integer(rowSums(!is.na(echino_wide[year_cols_tmp] %>%
                                    mutate(across(everything(), ~ suppressWarnings(as.integer(.x)))))) > 0) +
        as.integer(rowSums(!is.na(echino_wide[family_cols_tmp])) > 0)
    ) / 5
  ) %>%
  relocate(record_key, primary_source, sources_all, n_sources,
           coord_conflict, name_conflict,
           completeness_score)

cat("✅ Section 9 complete (depth handled in postprocessing).\n")

# =============================================================================
# SECTION 10: DIAGNOSTICS
# =============================================================================

cat("\n══ LAYER B DIAGNOSTICS ══════════════════════════════════\n")

cat("\n── Raw depth column coverage (pre-resolution) ──\n")
depth_min_cols_diag <- names(echino_wide)[str_detect(names(echino_wide), "^minimumDepthInMeters__")]
depth_max_cols_diag <- names(echino_wide)[str_detect(names(echino_wide), "^maximumDepthInMeters__")]

has_any_min <- rowSums(!is.na(echino_wide %>% select(all_of(depth_min_cols_diag)) %>%
                                mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))))) > 0
has_any_max <- rowSums(!is.na(echino_wide %>% select(all_of(depth_max_cols_diag)) %>%
                                mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))))) > 0

cat("Records with ANY minimumDepth value:", sum(has_any_min), "\n")
cat("Records with ANY maximumDepth value:", sum(has_any_max), "\n")
cat("Records with min OR max (raw, unresolved):", sum(has_any_min | has_any_max), "\n")
cat("(Depth resolution — min/max/best/conflicts/zones — done in postprocessing)\n")

cat("\n── Conflict summary (name + coordinate only) ──\n")
echino_wide %>%
  summarise(
    coord_conflicts = sum(coord_conflict, na.rm = TRUE),
    name_conflicts  = sum(name_conflict,  na.rm = TRUE),
    pct_coord = round(mean(coord_conflict, na.rm = TRUE) * 100, 1),
    pct_name  = round(mean(name_conflict,  na.rm = TRUE) * 100, 1)
  ) %>%
  pivot_longer(everything()) %>% print()

cat("\n── Multi-source records ──\n")
echino_wide %>%
  count(n_sources) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(desc(n_sources)) %>% print()

cat("\n── Primary source contribution ──\n")
echino_wide %>% count(primary_source) %>% arrange(desc(n)) %>% print()

cat("\n── Completeness summary ──\n")
spec_cols <- names(echino_wide)[str_detect(names(echino_wide), "^specificEpithet__")]
year_cols <- names(echino_wide)[str_detect(names(echino_wide), "^year__")]

has_year    <- rowSums(!is.na(echino_wide[year_cols] %>%
                                mutate(across(everything(), ~ suppressWarnings(as.integer(.x)))))) > 0
has_species <- rowSums(!is.na(echino_wide[spec_cols])) > 0
has_lat_vec <- rowSums(!is.na(echino_wide[lat_cols] %>%
                                mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))))) > 0
has_lon_vec <- rowSums(!is.na(echino_wide[lon_cols] %>%
                                mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))))) > 0

tibble(
  total             = nrow(echino_wide),
  with_coords       = sum(has_lat_vec & has_lon_vec),
  with_year         = sum(has_year),
  with_species      = sum(has_species),
  pct_coords        = round(with_coords  / total * 100, 1),
  pct_year          = round(with_year    / total * 100, 1),
  pct_species       = round(with_species / total * 100, 1),
  mean_completeness = round(mean(echino_wide$completeness_score, na.rm = TRUE) * 100, 1)
) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>% print()


# =============================================================================
# SECTION 11: SAVE OUTPUTS
# =============================================================================

write_csv(echino_wide, "echino_wide.csv")
cat("\n✅ LAYER B saved: echino_wide.csv\n")
cat("   Rows:", nrow(echino_wide), "| Columns:", ncol(echino_wide), "\n")

cat("\n🎉 Done! Two files ready:\n")
cat("   echino_long.csv  — all records, all columns, audit trail\n")
cat("   echino_wide.csv  — one row per specimen, conflicts visible (depth NOT yet resolved)\n")
cat("\nNext step: run echinoderm_postprocessing.R for full depth resolution,\n")
cat("name normalisation, and geographic/data-quality cleaning.\n")
