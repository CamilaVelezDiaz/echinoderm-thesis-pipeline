# =============================================================================
# TAXONOMIC IDENTIFICATION INDEX (TIa) - full analysis (Methods section 2.x.2)
# Camila Velez | JCU Marine Biology MSc | Thesis Chapter — Analysis
#
# Scope: TIa (% of records identified to species level) - the collapsed
# binary index, distinct from the six-level taxonomic completeness
# hierarchy covered in section 2.x.1 (taxonomic_completeness_full.R).
#
# NOTE: the log-depth specification used in the regression below was
# already justified in 2.x.1 via a GAM linearity check (edf = 8.8, non-
# linear; breakpoints at 22.5m and 84.5m) - that check is NOT repeated
# here, only referenced, to avoid duplicating work across subsections.
#
# Order (matches the 2.x.1 structure for consistency across subsections):
#   1. Descriptive: TIa by zone, by class, by class x zone (with Wilson
#      95% CIs), + figures
#   2. Chi-square test: TIa vs depth zone (simple, global)
#   3. Multivariable logistic regression (log-depth + class + quality flag)
#   4. Depth x class interaction test (likelihood ratio test)
#   5. Pairwise depth-zone comparisons with Holm correction + compact
#      letter groupings (relabelled so "a" = highest TIa)
#   6. Sensitivity analysis: robustness to excluding Flag 3 / Flags 2+3
#
# Packages needed beyond tidyverse: broom, multcompView (both optional -
# script degrades gracefully with clear messages if missing).
# Install if needed: install.packages(c("broom", "multcompView"))
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================

setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")

# =============================================================================
# Library
# =============================================================================
library(tidyverse)

# -----------------------------------------------------------------------------
# SETUP - always loads fresh from disk
# -----------------------------------------------------------------------------

echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE)

if (!"rank6" %in% names(echino_wide)) {
  rank_map <- c(
    "Species" = "Species", "Subspecies" = "Species",
    "Genus"   = "Genus",   "Subgenus"   = "Genus",
    "Family"  = "Family", "Order" = "Order",
    "Class"   = "Class",  "Phylum" = "Phylum"
  )
  echino_wide <- echino_wide %>%
    mutate(rank6 = recode(taxonomic_resolution_level, !!!rank_map))
}

if (!"crosses_zone_boundary" %in% names(echino_wide)) {
  echino_wide <- echino_wide %>%
    mutate(crosses_zone_boundary = str_starts(coalesce(depth_zone, ""), "Spans"))
}

zone_order_shallow_first <- c("Continental Shelf", "Upper Slope", "Lower Slope", "Abyssal")

wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(lower = NA_real_, upper = NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- (z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / denom
  c(lower = max(0, (center - half) * 100), upper = min(100, (center + half) * 100))
}

model_data <- echino_wide %>%
  filter(has_depth, !crosses_zone_boundary, !is.na(echino_class)) %>%
  mutate(
    species_level = as.integer(rank6 == "Species"),
    log_depth = log10(depth_median + 1),
    depth_bin = factor(depth_zone, levels = zone_order_shallow_first),
    echino_class = factor(echino_class),
    echino_class = relevel(echino_class, ref = "Ophiuroidea"),
    obs_quality_flag = factor(obs_quality_flag)
  )

# Separate, UNRESTRICTED dataset for the overall (non-class-stratified)
# zone-level TIa table below. model_data excludes the 828 depth-bearing
# records with unknown echino_class - a filter that's necessary for the
# class-stratified table and the regression models (which use class as a
# covariate) but has no bearing on the overall zone-level TIa, since
# whether a record reached species level doesn't depend on whether its
# class happens to be known. Using the unrestricted set here makes this
# table EXACTLY reproduce the Species% column in Taxonomic Completeness's
# Table 3 (2.x.1) - a deliberate cross-validation between the two sections,
# not a coincidence.
zone_data_all <- echino_wide %>%
  filter(has_depth, !crosses_zone_boundary) %>%
  mutate(
    species_level = as.integer(rank6 == "Species"),
    depth_bin = factor(depth_zone, levels = zone_order_shallow_first)
  )


# =============================================================================
# SECTION 1: DESCRIPTIVE - TIa by zone, by class, by class x zone (Wilson CIs)
# =============================================================================

cat("=== SECTION 1: TIa descriptive tables ===\n\n")

# --- TIa by depth zone (ALL depth-bearing records, not class-restricted -
# see note above on zone_data_all) ---
tia_by_zone <- zone_data_all %>%
  group_by(depth_bin) %>%
  summarise(n = n(), n_species = sum(species_level), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    TIa = round(n_species / n * 100, 1),
    ci = list(wilson_ci(n_species, n)),
    ci_lower = round(ci[["lower"]], 1),
    ci_upper = round(ci[["upper"]], 1)
  ) %>%
  ungroup() %>%
  dplyr::select(-ci)

cat("── TIa by depth zone (with 95% Wilson CI) - ALL depth-bearing records ──\n")
cat("(should exactly match Taxonomic Completeness Table 3's Species% column)\n\n")
print(tia_by_zone, n = Inf)
write_csv(tia_by_zone, "table3_tia_by_depth_with_CI.csv")

# --- TIa by class x depth zone ---
tia_class_ci <- model_data %>%
  group_by(echino_class, depth_bin) %>%
  summarise(n = n(), n_species = sum(species_level), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    TIa = round(n_species / n * 100, 1),
    ci = list(wilson_ci(n_species, n)),
    ci_lower = round(ci[["lower"]], 1),
    ci_upper = round(ci[["upper"]], 1)
  ) %>%
  ungroup() %>%
  dplyr::select(-ci)

cat("\n── TIa by class x depth zone (with 95% Wilson CI) ──\n")
print(tia_class_ci, n = Inf)
write_csv(tia_class_ci, "table2_tia_class_by_depth_with_CI.csv")

# --- Figure: TIa by zone, horizontal bar + CI (depth-ordered, shelf at top) ---
p_tia_h <- tia_by_zone %>%
  mutate(depth_bin_plot = factor(depth_bin, levels = rev(zone_order_shallow_first))) %>%
  ggplot(aes(x = depth_bin_plot, y = TIa)) +
  geom_col(fill = "#2a78d6", width = 0.65) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 0.6, color = "#0c447c") +
  geom_text(aes(label = sprintf("%.1f%%", TIa)), hjust = -0.4, vjust = -0.8, size = 3.2) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "TIa - records identified to species level (%, with 95% CI)",
       title = "Taxonomic identification completeness across the depth gradient") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

ggsave("figure_tia_by_depth_horizontal.png", p_tia_h, width = 7, height = 4, dpi = 300)
cat("\nFigure saved: figure_tia_by_depth_horizontal.png\n")

# --- Figure: TIa by class across depth zones, line + point + CI ---
p_tia_line <- tia_class_ci %>%
  mutate(depth_bin = factor(depth_bin, levels = zone_order_shallow_first)) %>%
  ggplot(aes(x = depth_bin, y = TIa, color = echino_class, group = echino_class)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.12,
                position = position_dodge(width = 0.3), linewidth = 0.6) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.3), size = 2.5) +
  scale_color_manual(values = c(
    "Asteroidea" = "#2a78d6", "Crinoidea" = "#1baf7a", "Echinoidea" = "#eda100",
    "Holothuroidea" = "#4a3aa7", "Ophiuroidea" = "#e34948"
  )) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Depth zone", y = "TIa (%, with 95% Wilson CI)", color = "Class",
       title = "TIa by echinoderm class across the depth gradient") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("figure_tia_by_class_lineplot_CI.png", p_tia_line, width = 9, height = 5.5, dpi = 300)
cat("Figure saved: figure_tia_by_class_lineplot_CI.png\n")


# =============================================================================
# SECTION 2: CHI-SQUARE TEST - TIa vs depth zone
# =============================================================================

cat("\n\n=== SECTION 2: Chi-square test (TIa vs depth zone) ===\n")
cat("H0: species-level identification rate is the same across all depth zones\n\n")
cat("Uses the unrestricted zone-level dataset (zone_data_all), matching\n")
cat("Table 3 above - this test doesn't involve class, so there's no reason\n")
cat("to exclude the 828 records with unknown echino_class here.\n\n")

zone_table <- table(zone_data_all$depth_bin, zone_data_all$species_level)
print(zone_table)

chisq_result <- chisq.test(zone_table)
print(chisq_result)

cat("\nNOTE: does not control for echino_class or obs_quality_flag, and\n")
cat("Abyssal's small n (20) carries little weight in this test -\n")
cat("that's what Section 3's regression is for (which DOES need to restrict\n")
cat("to known-class records, since class enters as a covariate there).\n")


# =============================================================================
# SECTION 3: LOGISTIC REGRESSION (continuous depth, controlling for class
# and quality flag). Log-depth specification justified in Section 2.x.1.
# =============================================================================

cat("\n\n=== SECTION 3: Logistic regression ===\n\n")

model_main <- glm(
  species_level ~ log_depth + echino_class + obs_quality_flag,
  data = model_data, family = binomial(link = "logit")
)
print(summary(model_main))

if (requireNamespace("broom", quietly = TRUE)) {
  or_table <- broom::tidy(model_main, exponentiate = TRUE, conf.int = TRUE) %>%
    rename(OR = estimate, ci_lower = conf.low, ci_upper = conf.high) %>%
    mutate(across(c(OR, ci_lower, ci_upper), ~ round(.x, 3)))
} else {
  or_table <- as_tibble(exp(cbind(OR = coef(model_main), confint(model_main))),
                        rownames = "term") %>%
    rename(ci_lower = `2.5 %`, ci_upper = `97.5 %`) %>%
    mutate(across(c(OR, ci_lower, ci_upper), ~ round(.x, 3)))
}

cat("\n── Odds ratios (95% CI) ──\n")
print(or_table, n = Inf)
write_csv(or_table, "logistic_regression_odds_ratios.csv")

# --- Model diagnostics: multicollinearity (VIF) and pseudo-R2 ---
cat("\n── Multicollinearity check (VIF) ──\n")
cat("Rule of thumb: VIF > 5 warrants attention, VIF > 10 indicates a\n")
cat("serious multicollinearity problem undermining individual coefficients.\n\n")

if (requireNamespace("car", quietly = TRUE)) {
  vif_result <- car::vif(model_main)
  print(vif_result)
  if (is.matrix(vif_result)) {
    # car::vif returns a matrix (with df-adjusted GVIF) for factors with >2 levels
    cat("\n(echino_class and obs_quality_flag have >2 levels, so this shows\n")
    cat("GVIF^(1/(2*Df)) - the standard adjusted metric, comparable to VIF.)\n")
  }
} else {
  cat("Package 'car' not installed - run install.packages('car') to compute VIF.\n")
}

# --- McFadden's pseudo-R2 ---
# 1 - (residual deviance / null deviance). NOTE: unlike OLS R2, McFadden's
# pseudo-R2 is typically much lower even for models with strong, genuine
# predictive power - McFadden (1974) suggested 0.2-0.4 already represents
# an EXCELLENT fit for this metric, not a weak one. Don't compare this
# number directly to an OLS R2 intuition.
pseudo_r2 <- 1 - (model_main$deviance / model_main$null.deviance)
cat("\n── McFadden's pseudo-R² ──\n")
cat("Pseudo-R² =", round(pseudo_r2, 3), "\n")
cat("(McFadden's pseudo-R2 is NOT directly comparable to OLS R2 - values of\n")
cat("0.2-0.4 are considered an excellent fit by this metric's own convention,\n")
cat("per McFadden 1974, not a weak one. Don't apply OLS-R2 intuition here.)\n")


# =============================================================================
# SECTION 4: DEPTH x CLASS INTERACTION TEST
# =============================================================================

cat("\n\n=== SECTION 4: Depth x class interaction (likelihood ratio test) ===\n")
cat("H0: the rate of decline in species-level ID with depth is the same for\n")
cat("    every echinoderm class\n\n")

model_interaction <- glm(
  species_level ~ log_depth * echino_class + obs_quality_flag,
  data = model_data, family = binomial(link = "logit")
)

lr_test <- anova(model_main, model_interaction, test = "LRT")
print(lr_test)

cat("\nIf p < .05: classes decline at genuinely different rates, not just\n")
cat("different starting points - backs up the class-specific slopes visible\n")
cat("in figure_tia_by_class_lineplot_CI.png.\n")

if (requireNamespace("broom", quietly = TRUE)) {
  interaction_or_table <- broom::tidy(model_interaction, exponentiate = TRUE, conf.int = TRUE) %>%
    rename(OR = estimate, ci_lower = conf.low, ci_upper = conf.high) %>%
    mutate(across(c(OR, ci_lower, ci_upper), ~ round(.x, 3)))
  write_csv(interaction_or_table, "logistic_regression_interaction_odds_ratios.csv")
}


# =============================================================================
# SECTION 5: PAIRWISE DEPTH-ZONE COMPARISONS + compact letter groupings
# (relabelled so "a" = highest TIa, since multcompLetters' default order
# is arbitrary and reads confusingly in a figure)
# =============================================================================

cat("\n\n=== SECTION 5: Pairwise depth-zone comparisons ===\n\n")
cat("Uses zone_data_all (unrestricted) to match Table 3 / Figure 1 exactly -\n")
cat("these letters are meant to annotate that figure, so they need the same n.\n\n")

zone_summary <- zone_data_all %>%
  group_by(depth_bin) %>%
  summarise(n = n(), n_species = sum(species_level), .groups = "drop") %>%
  mutate(TIa = round(n_species / n * 100, 1))

pairwise_result <- pairwise.prop.test(
  x = zone_summary$n_species, n = zone_summary$n, p.adjust.method = "holm"
)
dimnames(pairwise_result$p.value) <- list(
  zone_summary$depth_bin[-1], zone_summary$depth_bin[-nrow(zone_summary)]
)

cat("── Pairwise comparison p-values (Holm-adjusted) ──\n")
print(pairwise_result)

if (requireNamespace("multcompView", quietly = TRUE)) {
  pmat <- pairwise_result$p.value
  pvec <- c()
  for (i in seq_len(nrow(pmat))) {
    for (j in seq_len(ncol(pmat))) {
      if (!is.na(pmat[i, j])) {
        pvec[paste0(rownames(pmat)[i], "-", colnames(pmat)[j])] <- pmat[i, j]
      }
    }
  }
  raw_letters <- multcompView::multcompLetters(pvec)$Letters
  
  # Relabel: "a" = highest TIa, "b" = next, etc. (does not change which
  # zones are/aren't significantly different - only which letter is used)
  zone_summary_sorted <- zone_summary %>% arrange(desc(TIa))
  new_letters <- setNames(letters[seq_len(nrow(zone_summary_sorted))], zone_summary_sorted$depth_bin)
  
  letter_table <- tibble(depth_bin = names(raw_letters)) %>%
    mutate(letter = new_letters[depth_bin]) %>%
    left_join(zone_summary, by = "depth_bin") %>%
    arrange(desc(TIa))
  
  cat("\n── Compact letter groupings (a = highest TIa; zones sharing a\n")
  cat("    letter are NOT significantly different) ──\n")
  print(letter_table %>% dplyr::select(depth_bin, n, TIa, letter))
  write_csv(letter_table, "pairwise_zone_letters_corrected.csv")
} else {
  cat("\nPackage 'multcompView' not installed - run\n")
  cat("install.packages('multcompView') for automatic letter groupings.\n")
}


# =============================================================================
# SECTION 6: SENSITIVITY ANALYSIS - robustness to Flag 3 / Flags 2+3 exclusion
# =============================================================================

cat("\n\n=== SECTION 6: Sensitivity analysis ===\n\n")

fit_depth_model <- function(data, label, drop_flag = FALSE) {
  if (drop_flag) {
    m <- glm(species_level ~ log_depth + echino_class, data = data, family = binomial)
  } else {
    m <- glm(species_level ~ log_depth + echino_class + obs_quality_flag,
             data = data, family = binomial)
  }
  or <- exp(coef(m)["log_depth"])
  ci <- exp(confint(m, "log_depth", trace = FALSE))
  p <- summary(m)$coefficients["log_depth", "Pr(>|z|)"]
  tibble(model = label, n = nrow(data), log_depth_OR = round(or, 3),
         ci_lower = round(ci[1], 3), ci_upper = round(ci[2], 3), p_value = signif(p, 3))
}

sensitivity_results <- bind_rows(
  fit_depth_model(model_data, "All flags (1+2+3)"),
  fit_depth_model(model_data %>% filter(obs_quality_flag != 3), "Excluding Flag 3"),
  fit_depth_model(model_data %>% filter(obs_quality_flag == 1),
                  "Flag 1 only (excl. 2 and 3)", drop_flag = TRUE)
)

cat("Odds ratio for log_depth across model specifications:\n\n")
print(sensitivity_results, n = Inf)
write_csv(sensitivity_results, "sensitivity_depth_effect_by_flag_exclusion.csv")

cat("\nIf log_depth_OR stays similar across all rows, the depth effect is\n")
cat("robust to how citizen-science/survey records are handled.\n")

cat("\n\n=== TIa analysis complete ===\n")
cat("Outputs: table2_tia_class_by_depth_with_CI.csv, table3_tia_by_depth_with_CI.csv,\n")
cat("figure_tia_by_depth_horizontal.png, figure_tia_by_class_lineplot_CI.png,\n")
cat("logistic_regression_odds_ratios.csv, logistic_regression_interaction_odds_ratios.csv,\n")
cat("pairwise_zone_letters_corrected.csv, sensitivity_depth_effect_by_flag_exclusion.csv\n")