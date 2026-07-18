# ==============================================================================
# TEMPORAL TRENDS
# ==============================================================================
# Addresses: has collecting effort changed over time, has identification
# quality kept pace, and is species discovery in this collection still
# actively growing or has it plateaued historically?
#
# Five parts:
#   A. Records per decade (raw collecting-effort trend)
#   B. Proportion identified to species level per decade, with a formal
#      trend test (logistic regression, decade as continuous predictor -
#      the temporal analog of the depth-based TIa analysis) and a GAM
#      non-linearity check, since the relationship turns out to be a
#      dip-and-recovery pattern rather than a simple trend
#   C. Cumulative species discovery curve - when was each species FIRST
#      recorded, and how has cumulative species count grown over the
#      collection's history. Distinct from rarefaction (sample-size based,
#      not calendar-time based) - this asks whether discovery is still
#      actively growing historically, not whether more sampling would
#      likely turn up more species right now.
#   D. Depth zone composition of collecting effort by decade - has deep-
#      water sampling grown as a SHARE of total effort over time, not
#      just in raw record count (which Part A already covers, pooled
#      across all depths).
#   E. Identification lag - time between collection and identification,
#      using dateIdentified__* source fields (confirmed present across
#      all 12 sources). Resolved via the same majority-vote approach as
#      best_year, with a non-negative-lag sanity check (a negative lag
#      would mean "identified before collected", a data error).
#
# Input:  echino_wide.csv (must include best_year - resolved centrally in
#         ECHINODERM_POST-PROCESSING.R Section 3, the same way
#         best_latitude/best_longitude are resolved in Section 1a. If this
#         script stops with a missing-column error below, re-run the
#         post-processing pipeline first.)
# Output: records_per_decade.csv, tia_by_decade.csv,
#         tia_by_decade_trend_test.csv, cumulative_species_by_decade.csv,
#         volume_vs_species_novelty_by_decade.csv, depth_zone_by_decade.csv,
#         identification_lag_by_decade.csv, figure_records_per_decade.png,
#         figure_tia_by_decade.png, figure_cumulative_species_discovery.png,
#         figure_new_species_by_decade.png, figure_depth_zone_by_decade.png,
#         figure_identification_lag_by_decade.png
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
library(tidyr)
library(ggplot2)
library(broom)

echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE)

if (!"rank6" %in% names(echino_wide)) {
  rank_map <- c("Species"="Species","Subspecies"="Species","Genus"="Genus","Subgenus"="Genus",
                "Family"="Family","Order"="Order","Class"="Class","Phylum"="Phylum")
  echino_wide <- echino_wide %>% mutate(rank6 = recode(taxonomic_resolution_level, !!!rank_map))
}

if (!"best_year" %in% names(echino_wide)) {
  stop("best_year not found - run the Section 3 addition in ",
       "ECHINODERM_POST-PROCESSING.R (majority-vote year resolution) before this script.")
}

n_with_year <- sum(!is.na(echino_wide$best_year))
cat(sprintf("best_year found: %d of %d records resolved (%.1f%%), range %d-%d\n\n",
            n_with_year, nrow(echino_wide), 100 * n_with_year / nrow(echino_wide),
            min(echino_wide$best_year, na.rm = TRUE), max(echino_wide$best_year, na.rm = TRUE)))


# ==============================================================================
# PART A: RECORDS PER DECADE
# ==============================================================================

cat("\n\n=== PART A: records per decade ===\n\n")

echino_wide <- echino_wide %>%
  mutate(decade = if_else(!is.na(best_year), floor(best_year / 10) * 10, NA_real_))

records_per_decade <- echino_wide %>%
  filter(!is.na(decade)) %>%
  count(decade) %>%
  arrange(decade)

print(records_per_decade, n = Inf)
write_csv(records_per_decade, "records_per_decade.csv")

n_no_year <- sum(is.na(echino_wide$best_year))
cat(sprintf("\nRecords with no resolvable year: %d (%.1f%%) - excluded from all temporal analyses below.\n",
            n_no_year, 100 * n_no_year / nrow(echino_wide)))

p_decade <- ggplot(records_per_decade, aes(x = factor(decade), y = n)) +
  geom_col(fill = "#2a78d6") +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3, colour = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Decade", y = "Number of records",
       title = "Echinoderm records by decade of collection",
       subtitle = sprintf("n = %d records with resolvable year (%.1f%% of dataset)",
                          sum(records_per_decade$n), 100 * sum(records_per_decade$n) / nrow(echino_wide))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figure_records_per_decade.png", p_decade, width = 9, height = 6, dpi = 300)
cat("Figure saved: figure_records_per_decade.png\n")


# ==============================================================================
# PART B: PROPORTION IDENTIFIED TO SPECIES LEVEL, BY DECADE
# ==============================================================================
# Temporal analog of TIa - same logic (species-level vs not, Wilson CI,
# logistic regression trend test), applied to decade instead of depth.

cat("\n\n=== PART B: species-level identification rate by decade ===\n\n")

decade_data <- echino_wide %>%
  filter(!is.na(decade)) %>%
  mutate(species_level = as.integer(rank6 == "Species"))

tia_by_decade <- decade_data %>%
  group_by(decade) %>%
  summarise(n = n(), n_species = sum(species_level), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    ci = list(suppressWarnings(prop.test(n_species, n)$conf.int)),
    tia = round(100 * n_species / n, 1),
    ci_lower = round(100 * ci[[1]], 1),
    ci_upper = round(100 * ci[[2]], 1)
  ) %>%
  ungroup() %>%
  dplyr::select(-ci) %>%
  arrange(decade)

print(tia_by_decade, n = Inf)
write_csv(tia_by_decade, "tia_by_decade.csv")

# Chi-square test: does species-level ID rate differ across decades?
decade_table <- decade_data %>% count(decade, species_level) %>%
  pivot_wider(names_from = species_level, values_from = n, values_fill = 0)
chisq_table <- as.matrix(decade_table[, -1])
rownames(chisq_table) <- decade_table$decade

cat("\n── Chi-square test: species-level ID rate across decades ──\n")
chisq_decade <- chisq.test(chisq_table)
print(chisq_decade)

# Logistic regression trend test: decade as continuous predictor
decade_model <- glm(species_level ~ decade, data = decade_data, family = binomial)
model_summary <- tidy(decade_model, conf.int = TRUE, exponentiate = TRUE)

cat("\n── Logistic regression: species-level identification vs. decade (continuous) ──\n")
cat("OR interpretable as: change in odds of species-level ID per 1-decade increase\n\n")
print(model_summary)
write_csv(model_summary, "tia_by_decade_trend_test.csv")

# ------------------------------------------------------------------------------
# GAM check of non-linearity, mirroring the depth-based GAM check in
# Taxonomic Completeness Section 3. The descriptive table above shows a
# clear dip-and-recovery pattern (high in early decades, crashing in the
# 1990s, recovering by the 2020s) rather than a smooth trend - a linear
# model averages the dip against the recovery and can badly understate
# what's actually happening. Fitted on best_year directly (finer
# resolution than decade) rather than the binned decade variable.
# ------------------------------------------------------------------------------
library(mgcv)

gam_year <- gam(species_level ~ s(best_year), data = decade_data, family = binomial)

cat("\n── GAM check: is the year trend actually linear? ──\n")
cat("(edf close to 1 = linear; higher = curved/non-monotonic)\n\n")
print(summary(gam_year))

aic_linear <- AIC(glm(species_level ~ best_year, data = decade_data, family = binomial))
aic_gam <- AIC(gam_year)
cat(sprintf("\nLinear model AIC: %.1f\n", aic_linear))
cat(sprintf("GAM model AIC:    %.1f\n", aic_gam))

if (aic_linear - aic_gam > 10) {
  cat(sprintf(
    "\nRESULT: the year trend is NOT well-approximated by a straight line (ΔAIC = %.1f\n",
    aic_linear - aic_gam
  ))
  cat("in favour of the GAM) - the linear OR above should NOT be reported as the\n")
  cat("main finding here. Report the shape (dip/recovery pattern visible in the\n")
  cat("descriptive table and figure) instead, backed by this GAM comparison.\n")
} else {
  cat("\nRESULT: linear and GAM fits are comparable - the linear trend is a\n")
  cat("reasonable summary of the year effect.\n")
}

# Diagnostic: is the 1990s dip attributable to a particular data source
# dominating that decade? (worth checking before interpreting the dip as
# a genuine identification-quality decline, rather than a source artefact)
cat("\n── Diagnostic: source composition by decade (checking the 1990s dip) ──\n")
decade_data %>%
  count(decade, primary_source) %>%
  group_by(decade) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  filter(decade %in% c(1980, 1990, 2000)) %>%
  arrange(decade, desc(n)) %>%
  print(n = Inf)

cat("\n── Diagnostic: obs_quality_flag composition by decade ──\n")
decade_data %>%
  count(decade, obs_quality_flag) %>%
  group_by(decade) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  filter(decade %in% c(1980, 1990, 2000)) %>%
  arrange(decade, obs_quality_flag) %>%
  print(n = Inf)

# ------------------------------------------------------------------------------
# GAM controlling for obs_quality_flag: does the year effect survive once
# identification-quality tier is accounted for? obs_quality_flag is wrapped
# in factor() here (categorical: separate effect per tier) to match how it
# was treated in the depth-based TIa/Taxonomic Completeness models - passing
# it as a raw integer would incorrectly force a linear assumption (flag 3
# exactly "3x" flag 1), which is not how it was modelled elsewhere in this
# thesis.
# ------------------------------------------------------------------------------
gam_controlled <- gam(species_level ~ s(best_year) + factor(obs_quality_flag),
                      data = decade_data, family = binomial)

cat("\n── GAM controlling for obs_quality_flag ──\n")
cat("(tests whether the year effect is genuine, or just reflects which\n")
cat("quality tier happened to dominate which decade)\n\n")
print(summary(gam_controlled))

cat("\nIf s(best_year) remains highly significant with a similar edf to the\n")
cat("uncontrolled model above, the year effect is real and not attributable\n")
cat("to obs_quality_flag composition alone.\n")

p_decade_tia <- ggplot(tia_by_decade, aes(x = decade, y = tia)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#2a78d6", alpha = 0.2) +
  geom_line(colour = "#2a78d6", linewidth = 1) +
  geom_point(colour = "#2a78d6", size = 2) +
  labs(x = "Decade", y = "Species-level identification (%)",
       title = "Species-level identification rate by decade",
       subtitle = "Shaded band = 95% CI") +
  theme_minimal(base_size = 12)

ggsave("figure_tia_by_decade.png", p_decade_tia, width = 9, height = 6, dpi = 300)
cat("\nFigure saved: figure_tia_by_decade.png\n")


# ==============================================================================
# PART C: CUMULATIVE SPECIES DISCOVERY CURVE
# ==============================================================================
# For each species, find its FIRST recorded year in this collection, then
# build a cumulative species count over calendar time. Distinct from
# rarefaction: this asks whether species discovery is still historically
# growing, not whether more CURRENT sampling effort would likely find more.

cat("\n\n=== PART C: cumulative species discovery over time ===\n\n")

first_year_by_species <- echino_wide %>%
  filter(rank6 == "Species", !is.na(accepted_name_final), !is.na(best_year)) %>%
  group_by(accepted_name_final) %>%
  summarise(first_year = min(best_year), .groups = "drop")

cat(sprintf("Species with a resolvable first-recorded year: %d\n", nrow(first_year_by_species)))

cumulative_by_decade <- first_year_by_species %>%
  mutate(decade = floor(first_year / 10) * 10) %>%
  count(decade, name = "n_new_species") %>%
  arrange(decade) %>%
  mutate(cumulative_species = cumsum(n_new_species))

print(cumulative_by_decade, n = Inf)
write_csv(cumulative_by_decade, "cumulative_species_by_decade.csv")

n_total_species_dated <- sum(cumulative_by_decade$n_new_species)
cat(sprintf("\nTotal species with a dated first record: %d\n", n_total_species_dated))

p_cumulative <- ggplot(cumulative_by_decade, aes(x = decade, y = cumulative_species)) +
  geom_line(colour = "#2a78d6", linewidth = 1) +
  geom_point(colour = "#2a78d6", size = 2) +
  labs(x = "Decade", y = "Cumulative species recorded in the collection",
       title = "Species accumulation in the collection",
       subtitle = sprintf("n = %d species with a dated first record", n_total_species_dated)) +
  theme_minimal(base_size = 12)

ggsave("figure_cumulative_species_discovery.png", p_cumulative, width = 9, height = 6, dpi = 300)
cat("Figure saved: figure_cumulative_species_discovery.png\n")

# ------------------------------------------------------------------------------
# COMBINED SUMMARY TABLE: joins record volume (Part A) against species newly
# recorded (this part) per decade, so claims like "the 2000s had huge record
# volume but low species novelty" are directly checkable in one place rather
# than requiring the reader to cross-reference two separate tables/figures
# that were never shown side by side.
# ------------------------------------------------------------------------------
volume_vs_novelty <- records_per_decade %>%
  rename(n_records = n) %>%
  full_join(cumulative_by_decade, by = "decade") %>%
  mutate(
    n_new_species = coalesce(n_new_species, 0L),
    pct_of_records = round(100 * n_records / sum(n_records, na.rm = TRUE), 1),
    pct_of_new_species = round(100 * n_new_species / n_total_species_dated, 1),
    records_per_new_species = if_else(n_new_species > 0, round(n_records / n_new_species, 1), NA_real_)
  ) %>%
  arrange(decade)

cat("\n── Combined table: record volume vs. species newly recorded, by decade ──\n")
cat("(records_per_new_species = how many records it took to yield one newly-\n")
cat("recorded species that decade - a rough index of collecting 'efficiency'\n")
cat("for novel taxa; higher = more repeated sampling of already-known taxa)\n\n")
print(volume_vs_novelty, n = Inf)
write_csv(volume_vs_novelty, "volume_vs_species_novelty_by_decade.csv")

# ------------------------------------------------------------------------------
# NEW FIGURE: species newly recorded per decade (not cumulative) - makes the
# "which decades contributed the most" claim directly visible, which the
# cumulative curve above cannot show on its own (a flattening cumulative
# curve shows slowing OVERALL growth, but not which earlier decades drove
# the growth that did happen).
# ------------------------------------------------------------------------------
p_new_species <- ggplot(cumulative_by_decade, aes(x = factor(decade), y = n_new_species)) +
  geom_col(fill = "#2a78d6") +
  geom_text(aes(label = n_new_species), vjust = -0.4, size = 3, colour = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Decade", y = "Species newly recorded in the collection",
       title = "Species newly recorded in the collection, by decade",
       subtitle = sprintf("n = %d species with a dated first record", n_total_species_dated)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figure_new_species_by_decade.png", p_new_species, width = 9, height = 6, dpi = 300)
cat("Figure saved: figure_new_species_by_decade.png\n")

# Simple diagnostic: is the curve still rising steeply in the most recent
# decade, or has it visibly plateaued? (descriptive only, no formal test -
# the shape is best assessed visually from the figure)
last_two <- tail(cumulative_by_decade, 2)
if (nrow(last_two) == 2) {
  pct_gain_last_decade <- round(100 * last_two$n_new_species[2] / last_two$cumulative_species[1], 1)
  cat(sprintf("\nMost recent decade added %d new species (%.1f%% increase over the prior cumulative total) -\n",
              last_two$n_new_species[2], pct_gain_last_decade))
  cat("a large value suggests discovery is still actively growing; a small value suggests plateauing.\n")
}


# ==============================================================================
# PART D: HAS DEEP-SEA SAMPLING EFFORT INCREASED OVER TIME?
# ==============================================================================
# Same logic as Part A, but broken down by depth zone rather than pooled -
# asks whether deeper zones make up a GROWING SHARE of collecting effort in
# recent decades, not just whether raw record counts have grown (which
# Part A already shows for the whole dataset regardless of depth).

cat("\n\n=== PART D: depth zone composition by decade ===\n\n")

depth_by_decade <- echino_wide %>%
  filter(!is.na(decade), !crosses_zone_boundary, depth_zone != "No depth data") %>%
  count(decade, depth_zone) %>%
  group_by(decade) %>%
  mutate(pct_of_decade = round(n / sum(n) * 100, 1)) %>%
  ungroup()

print(depth_by_decade, n = Inf)
write_csv(depth_by_decade, "depth_zone_by_decade.csv")

p_depth_decade <- ggplot(depth_by_decade, aes(x = factor(decade), y = pct_of_decade, fill = depth_zone)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = c("Continental Shelf" = "#2a78d6", "Upper Slope" = "#eda100",
               "Lower Slope" = "#e34948", "Abyssal" = "#4a3aa7")
  ) +
  labs(x = "Decade", y = "% of depth-zoned records that decade", fill = "Depth zone",
       title = "Depth zone composition of collecting effort by decade",
       subtitle = "Among records with a resolved depth zone only") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figure_depth_zone_by_decade.png", p_depth_decade, width = 9, height = 6, dpi = 300)
cat("\nFigure saved: figure_depth_zone_by_decade.png\n")


# ==============================================================================
# PART E: IDENTIFICATION LAG - TIME BETWEEN COLLECTION AND IDENTIFICATION
# ==============================================================================
# dateIdentified__SOURCE columns exist across all 12 sources (confirmed via
# grep against echino_wide). Resolves a best_year_identified consensus
# column the same way best_year was resolved (majority vote across sources,
# matching Section 3's approach), then computes identification lag =
# year_identified - best_year (collection year), restricted to non-negative
# lags (a negative lag would mean "identified before it was collected",
# which is impossible and indicates a data error rather than a genuine
# same-year rapid identification - see diagnostic below).

cat("\n\n=== PART E: identification lag (collection year to identification year) ===\n\n")

dateid_cols <- names(echino_wide)[str_detect(names(echino_wide), "^dateIdentified__")]
cat(sprintf("Found %d dateIdentified source columns\n", length(dateid_cols)))

# Extract just the year from each dateIdentified column - values may be
# full dates (YYYY-MM-DD) or bare years; suppressWarnings + regex extraction
# handles both without erroring on already-bare-year values
extract_year <- function(x) {
  yr <- str_extract(x, "^(\\d{4})")
  suppressWarnings(as.numeric(yr))
}

dateid_mat <- echino_wide %>%
  dplyr::select(all_of(dateid_cols)) %>%
  mutate(across(everything(), extract_year))

cat("\nCoverage per dateIdentified source column:\n")
for (col in dateid_cols) {
  n_filled <- sum(!is.na(dateid_mat[[col]]))
  cat(sprintf("  %-35s %6d non-missing (%.1f%%)\n", col, n_filled, 100 * n_filled / nrow(echino_wide)))
}

modal_year_id <- function(row_vals) {
  vals <- na.omit(row_vals)
  if (length(vals) == 0) return(NA_real_)
  tbl <- sort(table(vals), decreasing = TRUE)
  as.numeric(names(tbl)[1])
}

echino_wide$best_year_identified <- apply(as.matrix(dateid_mat), 1, modal_year_id)

n_with_id_year <- sum(!is.na(echino_wide$best_year_identified))
cat(sprintf("\nRecords with best_year_identified resolved: %d of %d (%.1f%%)\n",
            n_with_id_year, nrow(echino_wide), 100 * n_with_id_year / nrow(echino_wide)))

# Compute lag only where BOTH collection year and identification year exist
lag_data <- echino_wide %>%
  filter(!is.na(best_year), !is.na(best_year_identified)) %>%
  mutate(identification_lag_years = best_year_identified - best_year)

cat(sprintf("\nRecords with both collection and identification year: %d\n", nrow(lag_data)))

n_negative_lag <- sum(lag_data$identification_lag_years < 0, na.rm = TRUE)
cat(sprintf("Records with NEGATIVE lag (identified before collected - data error): %d (%.1f%%)\n",
            n_negative_lag, 100 * n_negative_lag / nrow(lag_data)))
if (n_negative_lag > 0) {
  cat("These are excluded below as data errors (impossible chronology), not genuine\n")
  cat("same-year-or-faster identifications.\n")
}

lag_data_clean <- lag_data %>% filter(identification_lag_years >= 0)

cat(sprintf("\nRecords with valid (non-negative) identification lag: %d\n", nrow(lag_data_clean)))
cat("\n── Identification lag summary (years) ──\n")
lag_summary <- lag_data_clean %>%
  summarise(
    n = n(),
    median_lag = median(identification_lag_years),
    mean_lag = round(mean(identification_lag_years), 1),
    q25 = quantile(identification_lag_years, 0.25),
    q75 = quantile(identification_lag_years, 0.75),
    max_lag = max(identification_lag_years)
  )
print(lag_summary)

# Has lag changed over time? Group by COLLECTION decade (not identification
# decade) - asks "for specimens collected in decade X, how long did they
# typically wait to be identified"
cat("\n── Median identification lag BY collection decade ──\n")
lag_by_decade <- lag_data_clean %>%
  mutate(collection_decade = floor(best_year / 10) * 10) %>%
  group_by(collection_decade) %>%
  summarise(
    n = n(),
    median_lag = median(identification_lag_years),
    mean_lag = round(mean(identification_lag_years), 1),
    .groups = "drop"
  ) %>%
  arrange(collection_decade)

print(lag_by_decade, n = Inf)
write_csv(lag_by_decade, "identification_lag_by_decade.csv")

# GAM check: has the lag itself changed non-linearly over time, mirroring
# the same non-linearity check used for TIa above
gam_lag <- gam(identification_lag_years ~ s(best_year), data = lag_data_clean)
cat("\n── GAM check: has identification lag changed shape over time? ──\n")
print(summary(gam_lag))

p_lag <- ggplot(lag_by_decade, aes(x = collection_decade, y = median_lag)) +
  geom_line(colour = "#e34948", linewidth = 1) +
  geom_point(colour = "#e34948", size = 2) +
  labs(x = "Collection decade", y = "Median identification lag (years)",
       title = "Time between collection and identification, by collection decade") +
  theme_minimal(base_size = 12)

ggsave("figure_identification_lag_by_decade.png", p_lag, width = 9, height = 6, dpi = 300)
cat("\nFigure saved: figure_identification_lag_by_decade.png\n")

# ------------------------------------------------------------------------------
# COVERAGE CHECK: is having a dateIdentified value AT ALL biased by collection
# era? This is a different question from lag (which only looks at records
# where BOTH dates exist) - this asks whether older specimens are simply
# less likely to have this metadata recorded in the first place, which
# would itself be informative (a structural/provenance pattern, not a
# genuine change in identification speed).
# ------------------------------------------------------------------------------
cat("\n── Coverage: is dateIdentified availability biased by collection decade? ──\n\n")

coverage_by_decade <- echino_wide %>%
  filter(!is.na(decade)) %>%
  group_by(decade) %>%
  summarise(
    n = n(),
    n_with_id_date = sum(!is.na(best_year_identified)),
    pct_with_id_date = round(100 * mean(!is.na(best_year_identified)), 1),
    .groups = "drop"
  )
print(coverage_by_decade, n = Inf)
write_csv(coverage_by_decade, "identification_date_coverage_by_decade.csv")

cat("\n── Coverage: which sources actually populate dateIdentified? ──\n\n")
source_coverage <- tibble(source_col = dateid_cols) %>%
  mutate(
    n_filled = sapply(dateid_cols, function(col) sum(!is.na(dateid_mat[[col]]))),
    pct_filled = round(100 * n_filled / nrow(echino_wide), 1)
  ) %>%
  arrange(desc(n_filled))
print(source_coverage, n = Inf)

cat("\n── Does dateIdentified coverage vary by primary_source overall? ──\n\n")
echino_wide %>%
  count(primary_source, has_id_date = !is.na(best_year_identified)) %>%
  group_by(primary_source) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  filter(has_id_date) %>%
  arrange(desc(pct)) %>%
  ungroup() %>%
  dplyr::select(primary_source, n_with_date = n, pct) %>%
  print(n = Inf)



cat("\n=== Temporal trends analysis complete ===\n")
cat("Outputs: records_per_decade.csv, tia_by_decade.csv, tia_by_decade_trend_test.csv,\n")
cat("cumulative_species_by_decade.csv, volume_vs_species_novelty_by_decade.csv,\n")
cat("depth_zone_by_decade.csv, identification_lag_by_decade.csv,\n")
cat("figure_records_per_decade.png, figure_tia_by_decade.png,\n")
cat("figure_cumulative_species_discovery.png, figure_new_species_by_decade.png,\n")
cat("figure_depth_zone_by_decade.png, figure_identification_lag_by_decade.png\n")