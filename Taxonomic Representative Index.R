# ==============================================================================
# TAXONOMIC REPRESENTATIVENESS INDEX (TRi) - FULL PIPELINE
# ==============================================================================
# TRi = CSi / CSe x 100
#   CSe = species documented by AFD as occurring in the study region
#         (Queensland/Torres Strait/Gulf of Carpentaria/Coral Sea),
#         WoRMS-resolved to current accepted names
#   CSi = distinct CSe species also present in the occurrence dataset
#
# This is the SINGLE, FINAL version of the whole TRi pipeline, consolidating
# everything built and debugged across multiple sessions:
#   PART A: retrieve species lists per class from ALA, scrape AFD for
#           distribution data, apply known manual corrections
#   PART B: build CSe using reviewed IMCRA regions
#   PART C: WoRMS-resolve every CSe species to its current accepted name
#   PART D: compute CSi and TRi (overall, by class, depth-stratified
#           secondary), plus gap-species outputs in both directions
#
# Run top to bottom. Each part is checkpointed to disk, so re-running after
# an interruption skips work already done. To force a completely fresh
# pull, delete: afd_echinoderm_species_list.csv, afd_echinoderm_master.csv,
# cse_afd_checklist.csv, cse_worms_resolved.csv
# ==============================================================================

# =============================================================================
# SETUP
# =============================================================================

setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")

# =============================================================================
# Library
# =============================================================================
library(httr)
library(jsonlite)
library(rvest)
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tidyr)

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
CLASSES            <- c("Asteroidea", "Ophiuroidea", "Echinoidea", "Holothuroidea", "Crinoidea")
SPECIES_LIST_CSV   <- "afd_echinoderm_species_list.csv"
MASTER_CSV         <- "afd_echinoderm_master.csv"
CSE_AFD_CSV        <- "cse_afd_checklist.csv"
CSE_RESOLVED_CSV   <- "cse_worms_resolved.csv"
REQUEST_DELAY      <- 1.0


# ==============================================================================
# PART A: RETRIEVE SPECIES LISTS + SCRAPE AFD DISTRIBUTION DATA
# ==============================================================================

resolve_class_guid <- function(class_name) {
  url <- paste0("https://bie-ws.ala.org.au/ws/search.json?q=", class_name, "&fq=rank:class")
  resp <- try(GET(url, timeout(15)), silent = TRUE)
  if (inherits(resp, "try-error") || status_code(resp) != 200) {
    warning(paste("GUID lookup failed for", class_name)); return(NA_character_)
  }
  parsed <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
  results <- parsed$searchResults$results
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("No GUID found for", class_name)); return(NA_character_)
  }
  hit <- results %>% filter(str_to_lower(name) == str_to_lower(class_name))
  if (nrow(hit) == 0) hit <- results[1, ]
  hit$guid[1]
}

EMPTY_SPECIES_TBL <- tibble(class = character(), species = character())

# Queries ALA n_attempts times and takes the UNION of all results. ALA's
# species-rank download for a class GUID was found to sometimes return an
# incomplete result set (e.g. Holothuroidea returning 246 species with ZERO
# from genus Holothuria, when the correct complete list has 335 species
# including 47 Holothuria) - querying repeatedly and unioning recovers the
# complete list reliably.
fetch_class_species_list <- function(class_name, guid, n_attempts = 5) {
  if (is.na(guid)) return(EMPTY_SPECIES_TBL)
  
  url <- paste0("https://bie-ws.ala.org.au/ws/download?q=",
                URLencode(paste0('rkid_class:"', guid, '"')),
                "&fq=rank:species")
  
  all_results <- vector("list", n_attempts)
  for (i in seq_len(n_attempts)) {
    resp <- try(GET(url, timeout(60)), silent = TRUE)
    if (inherits(resp, "try-error") || status_code(resp) != 200) next
    
    raw_text <- content(resp, as = "text", encoding = "UTF-8")
    df <- try(read_csv(raw_text, show_col_types = FALSE), silent = TRUE)
    if (inherits(df, "try-error") || !"scientificName" %in% names(df)) next
    
    all_results[[i]] <- df %>%
      filter(taxonRank == "species") %>%
      transmute(class = class_name, species = scientificName)
    
    if (i < n_attempts) Sys.sleep(1.0)
  }
  
  combined <- bind_rows(all_results) %>% distinct(class, species)
  message(sprintf("    %s: union of %d attempts = %d distinct species",
                  class_name, n_attempts, nrow(combined)))
  if (nrow(combined) == 0) warning(paste("All attempts failed for", class_name))
  combined
}

message("=== PART A: species lists + AFD distribution scrape ===\n")

if (file.exists(SPECIES_LIST_CSV)) {
  message("Found existing species list, loading from disk")
  species_master <- read_csv(SPECIES_LIST_CSV, show_col_types = FALSE)
} else {
  species_master <- map_dfr(CLASSES, function(cls) {
    guid <- resolve_class_guid(cls)
    message(sprintf("  %-15s guid: %s", cls, ifelse(is.na(guid), "NOT FOUND", guid)))
    Sys.sleep(REQUEST_DELAY)
    fetch_class_species_list(cls, guid)
  })
  
  if (nrow(species_master) == 0 || !"species" %in% names(species_master)) {
    stop("All class downloads failed - check ALA API status before re-running.")
  }
  
  # Keep anything that starts like a genus name and doesn't contain a known
  # bad pattern (unidentified/provisional qualifiers, digits, hybrid marker).
  # Deliberately permissive of subgenus parentheses and subspecies trinomials -
  # an earlier, stricter version of this filter silently discarded 177
  # species (89 from Holothuroidea alone) by requiring an exact two-word
  # "Genus species" shape that subgenus/trinomial names don't match.
  species_master <- species_master %>%
    filter(
      !is.na(species),
      str_detect(species, "^[A-Z][a-z]+"),
      !str_detect(species, regex("\\bsp\\.?\\b|\\bspp\\.?\\b|\\bcf\\.|\\baff\\.|\\bnr\\.|\\d|\u00d7", ignore_case = TRUE))
    ) %>%
    distinct(class, species)
  
  write_csv(species_master, SPECIES_LIST_CSV)
}

message(sprintf("\n-> %d species across %d classes", nrow(species_master), n_distinct(species_master$class)))
print(count(species_master, class))

# --- Scraping helpers ---
extract_afd_field <- function(txt, label) {
  m <- str_match(txt, paste0(label, "\\n+([^\\n]+)"))
  m[1, 2]
}

parse_depth_range <- function(depth_str) {
  if (is.na(depth_str)) return(c(min = NA_real_, max = NA_real_))
  rng <- str_match(depth_str, "([\\d.]+)\\s*-\\s*([\\d.]+)\\s*m")
  if (!is.na(rng[1, 1])) return(c(min = as.numeric(rng[1, 2]), max = as.numeric(rng[1, 3])))
  single <- str_match(depth_str, "([\\d.]+)\\s*m")
  if (!is.na(single[1, 1])) return(c(min = NA_real_, max = as.numeric(single[1, 2])))
  c(min = NA_real_, max = NA_real_)
}

scrape_afd_species <- function(binomial) {
  # str_replace_ALL, not str_replace - a single-replace version silently
  # broke multi-space names like "Holothuria (Halodeima) atra", converting
  # only the first space and leaving a literal space in the URL.
  slug <- str_replace_all(binomial, " ", "_")
  url  <- paste0("https://biodiversity.org.au/afd/taxa/", slug)
  
  resp <- try(GET(url, timeout(15), user_agent("JCU thesis research")), silent = TRUE)
  if (inherits(resp, "try-error") || status_code(resp) != 200) {
    return(tibble(species = binomial, afd_url = url, scrape_status = "fetch_failed"))
  }
  txt <- try(read_html(resp) %>% html_text2(), silent = TRUE)
  if (inherits(txt, "try-error") || is.na(txt)) {
    return(tibble(species = binomial, afd_url = url, scrape_status = "parse_failed"))
  }
  
  states_block <- extract_afd_field(txt, "States")
  extra_block  <- extract_afd_field(txt, "Extra Distribution Information")
  imcra_block  <- extract_afd_field(txt, "IMCRA")
  eco_block    <- extract_afd_field(txt, "Ecological Descriptors")
  depth        <- str_extract(extra_block %||% "", "Depth range[^.]*m")
  depth_range  <- parse_depth_range(depth)
  
  tibble(
    species        = binomial,
    afd_url        = url,
    states         = states_block,
    qld_flag       = str_detect(states_block %||% "", regex("Queensland", ignore_case = TRUE)),
    imcra_regions  = imcra_block,
    depth_text     = depth,
    depth_min_m    = depth_range["min"],
    depth_max_m    = depth_range["max"],
    ecological_descriptors = eco_block,
    scrape_status  = "ok"
  )
}

message("\nScraping AFD species pages...")

already_done <- if (file.exists(MASTER_CSV)) read_csv(MASTER_CSV, show_col_types = FALSE) else tibble(species = character())
to_scrape <- species_master %>% anti_join(already_done, by = "species") %>% pull(species)
message(sprintf("%d already scraped, %d remaining", nrow(already_done), length(to_scrape)))

if (length(to_scrape) > 0) {
  results_new <- vector("list", length(to_scrape))
  for (i in seq_along(to_scrape)) {
    results_new[[i]] <- scrape_afd_species(to_scrape[i])
    Sys.sleep(REQUEST_DELAY)
    if (i %% 50 == 0 || i == length(to_scrape)) {
      message(sprintf("  ...%d / %d scraped", i, length(to_scrape)))
      write_csv(bind_rows(already_done, bind_rows(results_new)), MASTER_CSV)
    }
  }
  species_detail <- bind_rows(already_done, bind_rows(results_new))
} else {
  species_detail <- already_done
}

afd_master <- species_master %>%
  dplyr::select(class, species) %>%   # class from Step 1's authoritative source, not species_detail
  left_join(species_detail %>% dplyr::select(-any_of("class")), by = "species") %>%
  dplyr::select(class, species, afd_url, states, qld_flag, imcra_regions,
                depth_text, depth_min_m, depth_max_m, ecological_descriptors, scrape_status)

# ------------------------------------------------------------------------------
# Known manual corrections for species that need special handling to scrape
# correctly (confirmed by hand: disambiguation pages, spelling differences,
# subgenus placements, author-citation text baked into the name)
# ------------------------------------------------------------------------------
known_fixes <- tribble(
  ~species,                                              ~correct_url,
  "Asterodiscides soelae",     "https://biodiversity.org.au/afd/taxa/Asterodiscides_soleae",
  "Stephanometra tenuispina",  "https://biodiversity.org.au/afd/taxa/Stephanometra_tenuipinna",
  "Cheiraster teres",          "https://biodiversity.org.au/afd/taxa/Cheiraster_(Luidiaster)_teres",
  "Antedon loveni",            "https://biodiversity.org.au/afd/taxa/Colobometra_perspinosa",
  "Astropecten zebra",         "https://biodiversity.org.au/afd/taxa/Astropecten_vappa",
  "Macrophiothrix microplax",  "https://biodiversity.org.au/afd/taxa/Macrophiothrix_longipeda",
  "Deima validum",             "https://biodiversity.org.au/afd/taxa/Deima_validum_validum",
  "Oneirophanta mutabilis",    "https://biodiversity.org.au/afd/taxa/Oneirophanta_mutabilis_mutabilis",
  "Actinopyga mauritiana",     "https://biodiversity.org.au/afd/taxa/Actinopyga_mauritiana;Actinopyga",
  "Araeosoma coriaceum",       "https://biodiversity.org.au/afd/taxa/Araeosoma_coriaceum;Araeosoma",
  "Nanometra johnstoni",       "https://biodiversity.org.au/afd/taxa/Nanometra_johnstoni;Nanometra",
  "Ophiocreas sibogae",        "https://biodiversity.org.au/afd/taxa/Ophiocreas_sibogae;Ophiocreas",
  "Ophiomastix elegans",       "https://biodiversity.org.au/afd/taxa/Ophiomastix_elegans;Ophiomastix"
)

# Genuinely absent from AFD - valid in WoRMS, no AFD record exists at all
excluded_from_afd <- tribble(
  ~species,                   ~reason,
  "Pseudoceramaster glasbyi", "Valid WoRMS name, no matching record in AFD - likely not yet reflected in AFD's database"
)

to_fix <- afd_master %>% filter(species %in% known_fixes$species, scrape_status != "ok")
if (nrow(to_fix) > 0) {
  message(sprintf("\nApplying %d known manual correction(s)...", nrow(to_fix)))
  fix_results <- vector("list", nrow(to_fix))
  for (i in seq_len(nrow(to_fix))) {
    sp  <- to_fix$species[i]
    url <- known_fixes$correct_url[known_fixes$species == sp]
    resp <- try(GET(url, timeout(20), user_agent("JCU thesis research")), silent = TRUE)
    if (inherits(resp, "try-error") || status_code(resp) != 200) {
      fix_results[[i]] <- tibble(species = sp, afd_url = url, scrape_status = "fetch_failed")
      next
    }
    txt <- read_html(resp) %>% html_text2()
    states_block <- extract_afd_field(txt, "States")
    extra_block  <- extract_afd_field(txt, "Extra Distribution Information")
    imcra_block  <- extract_afd_field(txt, "IMCRA")
    eco_block    <- extract_afd_field(txt, "Ecological Descriptors")
    depth        <- str_extract(extra_block %||% "", "Depth range[^.]*m")
    depth_range  <- parse_depth_range(depth)
    fix_results[[i]] <- tibble(
      species = sp, afd_url = url, states = states_block,
      qld_flag = str_detect(states_block %||% "", regex("Queensland", ignore_case = TRUE)),
      imcra_regions = imcra_block, depth_text = depth,
      depth_min_m = depth_range["min"], depth_max_m = depth_range["max"],
      ecological_descriptors = eco_block, scrape_status = "ok"
    )
    message(sprintf("  [%d/%d] %s -> %s", i, nrow(to_fix), sp, fix_results[[i]]$scrape_status))
    Sys.sleep(1.5)
  }
  fix_results <- bind_rows(fix_results)
  afd_master <- afd_master %>%
    filter(!(species %in% fix_results$species)) %>%
    bind_rows(fix_results %>% dplyr::select(any_of(names(afd_master)))) %>%
    arrange(class, species)
}

write_csv(afd_master, MASTER_CSV)

message("\n=== Scrape summary ===")
message(sprintf("Total species: %d | ok: %d | still failing: %d",
                nrow(afd_master), sum(afd_master$scrape_status == "ok"),
                sum(afd_master$scrape_status != "ok")))
still_failing <- afd_master %>% filter(scrape_status != "ok") %>% pull(species)
if (length(still_failing) > 0) {
  message("Still failing (expected: only genuinely AFD-absent species):")
  print(still_failing)
}


# ==============================================================================
# PART B: BUILD CSe USING REVIEWED IMCRA REGIONS
# ==============================================================================

message("\n\n=== PART B: building CSe ===\n")

# These 11 regions were determined by inspecting the actual IMCRA labels
# present in this dataset (run the inventory below to re-verify if the
# underlying AFD data changes materially in the future).
QLD_RELEVANT_IMCRA <- c(
  "Northeast Shelf Province", "Northeast Shelf Transition",
  "Northeast Province", "Northeast Transition",
  "Central Eastern Shelf Province", "Central Eastern Shelf Transition",
  "Central Eastern Province", "Central Eastern Transition",
  "Northern Shelf Province",       # Gulf of Carpentaria
  "Kenn Province", "Kenn Transition"  # Coral Sea (Kenn Plateau)
)

# Print inventory for reference/re-verification
imcra_inventory <- afd_master %>%
  filter(!is.na(imcra_regions)) %>%
  separate_rows(imcra_regions, sep = ",\\s*") %>%
  mutate(imcra_label = str_trim(str_remove(imcra_regions, "\\s*\\(\\d+\\)$"))) %>%
  count(imcra_label, sort = TRUE)
message("IMCRA labels currently in use (top 10, for reference):")
print(head(imcra_inventory, 10))

species_imcra_flag <- afd_master %>%
  filter(!is.na(imcra_regions)) %>%
  separate_rows(imcra_regions, sep = ",\\s*") %>%
  mutate(imcra_label = str_trim(str_remove(imcra_regions, "\\s*\\(\\d+\\)$"))) %>%
  group_by(species) %>%
  summarise(imcra_qld_match = any(imcra_label %in% QLD_RELEVANT_IMCRA), .groups = "drop")

species_summary <- afd_master %>%
  group_by(species) %>%
  summarise(
    class         = first(na.omit(class)),
    states        = first(na.omit(states)),
    qld_flag      = any(qld_flag, na.rm = TRUE),
    depth_min_m   = suppressWarnings(min(depth_min_m, na.rm = TRUE)),
    depth_max_m   = suppressWarnings(max(depth_max_m, na.rm = TRUE)),
    scrape_status = first(na.omit(scrape_status)),
    .groups = "drop"
  ) %>%
  mutate(
    depth_min_m = if_else(is.infinite(depth_min_m), NA_real_, depth_min_m),
    depth_max_m = if_else(is.infinite(depth_max_m), NA_real_, depth_max_m)
  )

cse_afd <- species_summary %>%
  left_join(species_imcra_flag, by = "species") %>%
  mutate(
    imcra_qld_match = coalesce(imcra_qld_match, FALSE),
    in_cse = qld_flag | imcra_qld_match
  ) %>%
  dplyr::select(class, species, states, qld_flag, imcra_qld_match, in_cse,
                depth_min_m, depth_max_m, scrape_status)

write_csv(cse_afd, CSE_AFD_CSV)

message(sprintf("\nCSe (AFD, pre-WoRMS-resolution): %d species", sum(cse_afd$in_cse)))
print(cse_afd %>% filter(in_cse) %>% count(class, name = "n_species"))


# ==============================================================================
# PART C: WoRMS-RESOLVE EVERY CSe SPECIES TO ITS CURRENT ACCEPTED NAME
# ==============================================================================
# AFD's own taxonomic currency doesn't always match WoRMS (the authority the
# occurrence dataset is validated against). Confirmed cases: outdated
# synonyms, spelling differences, subgenus notation. Every CSe species is
# resolved here so matching against occurrence data uses consistent
# taxonomic currency on both sides.

message("\n\n=== PART C: WoRMS resolution ===\n")

lookup_worms_accepted_name <- function(binomial) {
  url <- paste0("https://www.marinespecies.org/rest/AphiaRecordsByName/",
                URLencode(binomial), "?marine_only=false")
  resp <- try(GET(url, timeout(15)), silent = TRUE)
  if (inherits(resp, "try-error") || status_code(resp) != 200) {
    return(list(status = "worms_fetch_failed", accepted_name = NA_character_))
  }
  raw_text <- content(resp, as = "text", encoding = "UTF-8")
  parsed <- try(fromJSON(raw_text, flatten = TRUE), silent = TRUE)
  if (inherits(parsed, "try-error") || length(parsed) == 0 ||
      (is.data.frame(parsed) && nrow(parsed) == 0)) {
    return(list(status = "worms_no_match", accepted_name = NA_character_))
  }
  df <- if (is.data.frame(parsed)) parsed else bind_rows(parsed)
  if ("status" %in% names(df)) {
    accepted_row <- df %>% filter(status == "accepted")
    if (nrow(accepted_row) >= 1) df <- accepted_row
  }
  rec <- df[1, ]
  list(status = "ok", accepted_name = rec$valid_name %||% rec$scientificname %||% NA_character_)
}

cse <- read_csv(CSE_AFD_CSV, show_col_types = FALSE) %>% filter(in_cse)

already_resolved <- if (file.exists(CSE_RESOLVED_CSV)) {
  read_csv(CSE_RESOLVED_CSV, show_col_types = FALSE) %>%
    dplyr::select(species, worms_accepted_name, worms_match_status)  # lookup columns ONLY -
  # keeping other cached columns here risks a class/class collision on
  # the later left_join if the input CSe schema ever changes
} else {
  tibble(species = character())
}

to_resolve <- cse %>% anti_join(already_resolved, by = "species") %>% pull(species)
message(sprintf("%d already resolved, %d remaining", nrow(already_resolved), length(to_resolve)))

if (length(to_resolve) == 0 && nrow(already_resolved) != nrow(cse)) {
  warning("Cached WoRMS resolution size doesn't match current CSe size but 0 species ",
          "need resolving - the cache is likely stale. Delete ", CSE_RESOLVED_CSV, " and re-run.")
}

if (length(to_resolve) > 0) {
  results <- vector("list", length(to_resolve))
  for (i in seq_along(to_resolve)) {
    lookup <- lookup_worms_accepted_name(to_resolve[i])
    results[[i]] <- tibble(
      species = to_resolve[i],
      worms_accepted_name = if (lookup$status == "ok") lookup$accepted_name else NA_character_,
      worms_match_status = lookup$status
    )
    Sys.sleep(REQUEST_DELAY)
    if (i %% 50 == 0 || i == length(to_resolve)) {
      message(sprintf("  ...%d / %d resolved", i, length(to_resolve)))
      write_csv(bind_rows(already_resolved, bind_rows(results)), CSE_RESOLVED_CSV)
    }
  }
  resolved_lookup <- bind_rows(already_resolved, bind_rows(results))
} else {
  resolved_lookup <- already_resolved
}

cse_resolved <- cse %>%
  left_join(resolved_lookup, by = "species") %>%
  mutate(
    match_name = coalesce(worms_accepted_name, species),
    name_was_updated = !is.na(worms_accepted_name) & worms_accepted_name != species
  )
write_csv(cse_resolved, CSE_RESOLVED_CSV)

message(sprintf("\nWoRMS resolution: %d succeeded, %d name(s) updated, %d fell back to AFD name",
                sum(cse_resolved$worms_match_status == "ok", na.rm = TRUE),
                sum(cse_resolved$name_was_updated, na.rm = TRUE),
                sum(cse_resolved$worms_match_status != "ok" | is.na(cse_resolved$worms_match_status))))

dup_check <- cse_resolved %>% count(match_name) %>% filter(n > 1)
if (nrow(dup_check) > 0) {
  message(sprintf(
    "Note: %d match_name(s) appear twice post-resolution (two AFD entries resolving to the same current species) - CSe uses DISTINCT match_name.",
    nrow(dup_check)
  ))
}

cse_final_n <- n_distinct(cse_resolved$match_name)
message(sprintf("\nFINAL CSe (distinct species, WoRMS-resolved): %d", cse_final_n))


# ==============================================================================
# PART D: COMPUTE CSi AND TRi
# ==============================================================================

message("\n\n=== PART D: computing CSi and TRi ===\n")

normalize_name <- function(x) {
  x %>% str_replace_all("\\([A-Za-z]+\\)", " ") %>% str_squish() %>% str_to_lower()
}
normalize_name_binomial_only <- function(x) {
  x %>% str_replace_all("\\([A-Za-z]+\\)", " ") %>% str_squish() %>% str_to_lower() %>% word(1, 2)
}

# ------------------------------------------------------------------------------
# KNOWN SPELLING ALIASES: confirmed by direct cross-comparison of csi_gap_species
# (CSe not collected) against occurrence_only_species (collected not in CSe).
# Each pair represents the SAME real species, where AFD's spelling differs
# from WoRMS' current spelling by one or a few letters - close enough to be
# certainly the same taxon (verified manually), but different enough that
# neither normalize_name() nor the WoRMS cross-reference step could catch
# the mismatch automatically. Without this, both species appear as FALSE
# gaps on BOTH lists simultaneously. Verify manually before adding new pairs.
# ------------------------------------------------------------------------------
known_spelling_aliases <- tribble(
  ~variant_a,                  ~variant_b,
  "stichopus hermanni",        "stichopus herrmanni",
  "stephanometra tenuispina",  "stephanometra tenuipinna",
  "leiaster leachii",          "leiaster leachi",
  "ophiocentrus asperus",      "ophiocentrus aspera"
)

resolve_known_aliases <- function(key) {
  match_idx <- match(key, known_spelling_aliases$variant_b)
  if_else(!is.na(match_idx), known_spelling_aliases$variant_a[match_idx], key)
}

echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE)

if (!"rank6" %in% names(echino_wide)) {
  rank_map <- c("Species"="Species","Subspecies"="Species","Genus"="Genus","Subgenus"="Genus",
                "Family"="Family","Order"="Order","Class"="Class","Phylum"="Phylum")
  echino_wide <- echino_wide %>% mutate(rank6 = recode(taxonomic_resolution_level, !!!rank_map))
}
if (!"crosses_zone_boundary" %in% names(echino_wide)) {
  echino_wide <- echino_wide %>% mutate(crosses_zone_boundary = str_starts(coalesce(depth_zone, ""), "Spans"))
}

occurrence_species <- echino_wide %>%
  filter(rank6 == "Species", !is.na(accepted_name_final)) %>%
  distinct(accepted_name_final) %>%
  mutate(match_key = resolve_known_aliases(normalize_name(accepted_name_final)),
         match_key_binomial = resolve_known_aliases(normalize_name_binomial_only(accepted_name_final)))

cse_normalized <- cse_resolved %>%
  distinct(match_name, class) %>%
  mutate(match_key = resolve_known_aliases(normalize_name(match_name)),
         match_key_binomial = resolve_known_aliases(normalize_name_binomial_only(match_name)))

cse_total <- n_distinct(cse_normalized$match_key)

matched_full <- cse_normalized %>% inner_join(occurrence_species, by = "match_key")
cse_unmatched_full <- cse_normalized %>% anti_join(occurrence_species, by = "match_key")
matched_binomial_fallback <- cse_unmatched_full %>%
  dplyr::select(-match_key) %>%
  inner_join(occurrence_species %>% dplyr::select(-match_key), by = "match_key_binomial") %>%
  distinct(match_key_binomial, class, .keep_all = TRUE)

matched_species <- bind_rows(
  matched_full %>% mutate(matched_via = "full_name"),
  matched_binomial_fallback %>% mutate(matched_via = "binomial_fallback")
) %>% mutate(csi_id = coalesce(match_key, match_key_binomial))

csi_total <- n_distinct(matched_species$csi_id)
tri_overall <- round(csi_total / cse_total * 100, 1)

message(sprintf("CSe = %d | CSi = %d | TRi = %.1f%%", cse_total, csi_total, tri_overall))

write_csv(tibble(metric = c("CSe", "CSi", "TRi (%)"),
                 value = c(cse_total, csi_total, tri_overall)),
          "TRi_overall_summary.csv")

# By class
tri_by_class <- cse_normalized %>%
  group_by(class) %>% summarise(CSe = n_distinct(match_key), .groups = "drop") %>%
  left_join(matched_species %>% group_by(class) %>% summarise(CSi = n_distinct(csi_id), .groups = "drop"),
            by = "class") %>%
  mutate(CSi = coalesce(CSi, 0L), TRi = round(CSi / CSe * 100, 1)) %>%
  arrange(desc(TRi))
print(tri_by_class, n = Inf)
write_csv(tri_by_class, "TRi_by_class.csv")

# Gap species (both directions)
matched_cse_keys <- unique(c(matched_full$match_key, matched_binomial_fallback$match_key_binomial))

gap_species <- cse_normalized %>%
  filter(!(match_key %in% matched_cse_keys), !(match_key_binomial %in% matched_cse_keys)) %>%
  dplyr::select(class, species = match_name) %>% arrange(class, species)
write_csv(gap_species, "csi_gap_species.csv")

occurrence_record_counts <- echino_wide %>%
  filter(rank6 == "Species", !is.na(accepted_name_final)) %>%
  count(accepted_name_final, name = "n_records")
occurrence_only_species <- occurrence_species %>%
  filter(!(match_key %in% matched_cse_keys), !(match_key_binomial %in% matched_cse_keys)) %>%
  left_join(occurrence_record_counts, by = "accepted_name_final") %>%
  dplyr::select(species = accepted_name_final, n_records) %>% arrange(desc(n_records))
write_csv(occurrence_only_species, "occurrence_only_species.csv")

message(sprintf("\nGap species: %d in CSe not collected | %d collected not in CSe",
                nrow(gap_species), nrow(occurrence_only_species)))

# Secondary: depth-stratified TRi
zone_bounds <- tribble(
  ~depth_bin,           ~zone_min, ~zone_max,
  "Continental Shelf",  0,         200,
  "Upper Slope",        200,       1000,
  "Lower Slope",        1000,      4000,
  "Abyssal",    4000,      Inf
)

cse_with_depth <- cse_normalized %>%
  left_join(cse_resolved %>% distinct(match_name, depth_min_m, depth_max_m) %>%
              mutate(match_key = resolve_known_aliases(normalize_name(match_name))), by = "match_key") %>%
  filter(!is.na(depth_min_m) | !is.na(depth_max_m))

cse_zone_expected <- cse_with_depth %>%
  mutate(d_min = coalesce(depth_min_m, 0), d_max = coalesce(depth_max_m, 6000)) %>%
  dplyr::cross_join(zone_bounds) %>%
  filter(d_min <= zone_max, d_max >= zone_min) %>%
  distinct(match_key, class, depth_bin)

occurrence_species_by_zone <- echino_wide %>%
  filter(rank6 == "Species", !is.na(accepted_name_final), !crosses_zone_boundary,
         depth_zone %in% zone_bounds$depth_bin) %>%
  distinct(depth_zone, accepted_name_final) %>%
  rename(depth_bin = depth_zone) %>%
  mutate(match_key = resolve_known_aliases(normalize_name(accepted_name_final)))

tri_by_zone <- cse_zone_expected %>%
  group_by(depth_bin) %>% summarise(CSe_zone = n_distinct(match_key), .groups = "drop") %>%
  left_join(
    cse_zone_expected %>% inner_join(occurrence_species_by_zone, by = c("depth_bin", "match_key")) %>%
      group_by(depth_bin) %>% summarise(CSi_zone = n_distinct(match_key), .groups = "drop"),
    by = "depth_bin"
  ) %>%
  mutate(CSi_zone = coalesce(CSi_zone, 0L), TRi_zone = round(CSi_zone / CSe_zone * 100, 1),
         depth_bin = factor(depth_bin, levels = zone_bounds$depth_bin)) %>%
  arrange(depth_bin)

depth_coverage_pct <- round(100 * n_distinct(cse_with_depth$match_key) / cse_total, 1)
message(sprintf("\nDepth-stratified TRi (SECONDARY - %d/%d CSe species with depth data, %s%%):",
                n_distinct(cse_with_depth$match_key), cse_total, depth_coverage_pct))
print(tri_by_zone, n = Inf)
write_csv(tri_by_zone, "TRi_depth_stratified_secondary.csv")

message("\n=== TRi pipeline complete ===")
message("Outputs: TRi_overall_summary.csv, TRi_by_class.csv, csi_gap_species.csv,")
message("occurrence_only_species.csv, TRi_depth_stratified_secondary.csv")

# ==============================================================================
# TRi BY ECHINODERM CLASS - FIGURE
# ==============================================================================
library(dplyr)
library(ggplot2)

tri_by_class <- tibble::tribble(
  ~class,          ~CSe, ~CSi,
  "Ophiuroidea",    283,  235,
  "Crinoidea",       83,   67,
  "Echinoidea",     144,  115,
  "Holothuroidea",  171,  136,
  "Asteroidea",     192,  150
) %>%
  mutate(
    TRi = round(100 * CSi / CSe, 1),
    class = factor(class, levels = class[order(-TRi)]),
    label_y = TRi + 6   # offset increased so 2-line labels clear the dashed line even for bars close to it
  )

overall_TRi <- round(100 * sum(tri_by_class$CSi) / sum(tri_by_class$CSe), 1)

cat(sprintf("Overall TRi (reference line): %.1f%% (CSi=%d, CSe=%d)\n",
            overall_TRi, sum(tri_by_class$CSi), sum(tri_by_class$CSe)))

p_tri_class <- ggplot(tri_by_class, aes(x = class, y = TRi)) +
  geom_col(fill = "#2a78d6") +
  geom_hline(yintercept = overall_TRi, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_text(aes(y = label_y, label = sprintf("%.1f%%\n(%d/%d)", TRi, CSi, CSe)),
            size = 3, colour = "grey20", lineheight = 0.9) +
  scale_y_continuous(limits = c(0, max(tri_by_class$label_y) * 1.15),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Echinoderm class", y = "Taxonomic Representativeness Index (TRi, %)",
       title = "Taxonomic Representativeness Index by echinoderm class",
       subtitle = sprintf("Dashed line = overall TRi (%.1f%%, CSi=%d of CSe=%d); chi-square test: not significant (p = .737)",
                          overall_TRi, sum(tri_by_class$CSi), sum(tri_by_class$CSe))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("figure_TRi_by_class.png", p_tri_class, width = 9, height = 6, dpi = 300)
cat("Figure saved: figure_TRi_by_class.png\n")

print(tri_by_class %>% dplyr::select(class, CSe, CSi, TRi) %>% arrange(desc(TRi)))