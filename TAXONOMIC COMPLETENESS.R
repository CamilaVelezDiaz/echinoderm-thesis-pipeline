# =============================================================================
# TAXONOMIC COMPLETENESS - full analysis (Methods section 2.x.1)
# Camila Velez | JCU Marine Biology MSc | Thesis Chapter — Analysis
#
# Scope: the SIX-LEVEL RESOLUTION HIERARCHY (Species/Genus/Family/Order/
# Class/Phylum) as its own outcome - NOT the TIa binary index, which is a
# separate subsection (2.x.2) with its own script.
#
# Order (matches the Methods subsection structure):
#   1. Descriptive: rank composition by depth zone (table + figure)
#   2. Chi-square test on the full rank x depth-zone table (global test,
#      no assumptions about depth's functional form)
#   3. GAM check of the depth relationship's shape - justifies how depth
#      enters the regression models below - + breakpoint estimation if
#      non-linearity is found (pooled across all classes)
#   3b. Class-specific version of the same check: does the curve SHAPE
#      differ by echinoderm class? Uses a constrained smooth (k=4) to
#      avoid overfitting sparse classes, plus per-class breakpoints via
#      simpler single-breakpoint segmented models.
#   4. Ordinal logistic regression (proportional odds model), using the
#      depth specification justified in step 3
#   5. Proportional odds assumption check (Brant test); given a violation,
#      cutpoint-specific binary models reported as the corrected result
#
# NOTE: the GAM/breakpoint result from step 3 also justifies the log-depth
# specification used in the TIa logistic regression (section 2.x.2) - it
# is deliberately only run ONCE, here, and referenced (not repeated) there.
#
# Packages needed beyond tidyverse: mgcv, segmented, brant (all optional -
# script degrades gracefully with clear messages if any are missing).
# Install if needed: install.packages(c("mgcv", "segmented", "brant"))
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
# SETUP - always loads fresh from disk (see earlier scripts for why this
# matters: reusing an in-memory echino_wide can silently be a filtered copy
# left behind by a previous script in the same session)
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
rank_order <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")

echino_wide_depth <- echino_wide %>%
  filter(!crosses_zone_boundary, depth_zone %in% zone_order_shallow_first) %>%
  mutate(depth_bin_shallow_first = factor(depth_zone, levels = zone_order_shallow_first))


# =============================================================================
# SECTION 1: DESCRIPTIVE - rank composition by depth zone (table + figure)
# =============================================================================

cat("=== SECTION 1: Rank composition by depth zone (descriptive) ===\n\n")

table1_rank_by_depthbin <- echino_wide_depth %>%
  count(depth_bin_shallow_first, rank6) %>%
  group_by(depth_bin_shallow_first) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup() %>%
  dplyr::select(depth_bin_shallow_first, rank6, pct) %>%
  pivot_wider(names_from = rank6, values_from = pct, values_fill = 0) %>%
  dplyr::select(depth_bin_shallow_first, any_of(rank_order)) %>%
  left_join(echino_wide_depth %>% count(depth_bin_shallow_first, name = "n_records"),
            by = "depth_bin_shallow_first") %>%
  rename(depth_bin = depth_bin_shallow_first)

cat("── Table 1: % records by taxonomic rank, per depth zone ──\n")
print(table1_rank_by_depthbin, n = Inf)
write_csv(table1_rank_by_depthbin, "table1_rank_by_depthbin.csv")

# Figure: 100% stacked horizontal bar, depth-ordered shelf-at-top
rank_by_zone <- echino_wide_depth %>%
  count(depth_bin_shallow_first, rank6) %>%
  group_by(depth_bin_shallow_first) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    rank6 = factor(rank6, levels = rev(rank_order)),
    depth_bin_plot = factor(depth_bin_shallow_first, levels = rev(zone_order_shallow_first))
  )

p_rank <- rank_by_zone %>%
  ggplot(aes(x = depth_bin_plot, y = pct, fill = rank6)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "Phylum" = "#898781", "Class" = "#e34948", "Order" = "#7b3294",
    "Family" = "#eda100", "Genus" = "#1baf7a", "Species" = "#2a78d6"
  ), breaks = rank_order) +
  labs(x = NULL, y = "% of records", fill = "Rank",
       title = "Taxonomic rank composition across the depth gradient") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

ggsave("figure_rank_composition_by_depth.png", p_rank, width = 8, height = 4, dpi = 300)
cat("\nFigure saved: figure_rank_composition_by_depth.png\n")


# =============================================================================
# SECTION 2: CHI-SQUARE TEST on the full rank x depth-zone table
# =============================================================================

cat("\n\n=== SECTION 2: Chi-square test (full rank x depth-zone table) ===\n")
cat("H0: the distribution across all 6 ranks is the same across depth zones\n\n")

zone_data <- echino_wide %>%
  filter(has_depth, !crosses_zone_boundary, depth_zone %in% zone_order_shallow_first) %>%
  mutate(depth_bin = factor(depth_zone, levels = zone_order_shallow_first))

full_table <- table(zone_data$depth_bin, zone_data$rank6)
cat("── Full contingency table ──\n")
print(full_table)

chisq_standard <- chisq.test(full_table)
n_low_expected <- sum(chisq_standard$expected < 5)

cat("\nCells with expected count < 5:", n_low_expected, "of", length(chisq_standard$expected), "\n")

if (n_low_expected > 0) {
  cat("\n\u26A0\uFE0F  Standard chi-square not valid (sparse cells) - using Monte Carlo\n")
  cat("   simulated p-value instead:\n\n")
  set.seed(42)
  chisq_result <- chisq.test(full_table, simulate.p.value = TRUE, B = 10000)
  print(chisq_result)
  cat("\nUSE THIS p-value in the thesis (cite: Monte Carlo simulation, B=10000).\n")
} else {
  chisq_result <- chisq_standard
  print(chisq_result)
  cat("\nAll expected counts >= 5 - standard chi-square is valid as-is.\n")
}


# =============================================================================
# SECTION 3: GAM CHECK of the depth relationship's shape + breakpoint
# estimation. This justifies the log-depth specification used in BOTH the
# ordinal model below AND the TIa logistic regression in section 2.x.2 -
# run once here, referenced (not repeated) there.
# =============================================================================

cat("\n\n=== SECTION 3: GAM check of depth relationship shape ===\n\n")

model_data <- echino_wide %>%
  filter(has_depth, !crosses_zone_boundary, !is.na(echino_class)) %>%
  mutate(
    species_level = as.integer(rank6 == "Species"),
    log_depth = log10(depth_median + 1),
    echino_class = factor(echino_class),
    echino_class = relevel(echino_class, ref = "Ophiuroidea"),
    obs_quality_flag = factor(obs_quality_flag)
  )

if (requireNamespace("mgcv", quietly = TRUE)) {
  library(mgcv)
  
  model_linear <- glm(
    species_level ~ log_depth + echino_class + obs_quality_flag,
    data = model_data, family = binomial
  )
  model_gam <- gam(
    species_level ~ s(log_depth) + echino_class + obs_quality_flag,
    data = model_data, family = binomial
  )
  
  cat("── Smooth term summary (edf close to 1 = linear; higher = curved) ──\n\n")
  print(summary(model_gam))
  
  cat("\n── Model comparison ──\n")
  cat("Linear model AIC:", round(AIC(model_linear), 1), "\n")
  cat("GAM model AIC:   ", round(AIC(model_gam), 1), "\n")
  
  edf <- summary(model_gam)$s.table[1, "edf"]
  delta_aic <- AIC(model_linear) - AIC(model_gam)
  
  if (edf > 2 && delta_aic > 2) {
    cat("\nRESULT: relationship is NOT well-approximated by a straight line\n")
    cat("(edf =", round(edf, 1), ", \u0394AIC =", round(delta_aic, 1), "in favour of the GAM).\n")
    cat("Proceeding to estimate the breakpoint location.\n")
  } else {
    cat("\nRESULT: linear log-depth specification appears reasonable\n")
    cat("(edf close to 1 and/or AIC difference small) - no breakpoint needed.\n")
  }
  
  # Save the fitted curve for the figure
  depth_range <- seq(min(model_data$log_depth), max(model_data$log_depth), length.out = 200)
  pred_data <- tibble(
    log_depth = depth_range,
    echino_class = factor("Ophiuroidea", levels = levels(model_data$echino_class)),
    obs_quality_flag = factor("1", levels = levels(model_data$obs_quality_flag))
  )
  pred <- predict(model_gam, newdata = pred_data, type = "response", se.fit = TRUE)
  pred_data <- pred_data %>%
    filter(is.finite(log_depth)) %>%
    mutate(
      depth_m = 10^log_depth - 1,
      fit = pred$fit[is.finite(log_depth)],
      lower = pmax(0, pred$fit[is.finite(log_depth)] - 1.96 * pred$se.fit[is.finite(log_depth)]),
      upper = pmin(1, pred$fit[is.finite(log_depth)] + 1.96 * pred$se.fit[is.finite(log_depth)])
    )
  
  p_gam <- pred_data %>%
    filter(depth_m > 0) %>%
    ggplot(aes(x = depth_m, y = fit * 100)) +
    geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100), fill = "#2a78d6", alpha = 0.2) +
    geom_line(color = "#2a78d6", linewidth = 1) +
    scale_x_log10(labels = scales::comma) +
    labs(
      x = "Depth (m, log scale)", y = "Predicted probability of species-level ID (%)",
      title = "GAM-fitted depth relationship (Ophiuroidea, Flag 1 reference)"
    ) +
    theme_minimal(base_size = 12)
  
  ggsave("figure_gam_depth_shape.png", p_gam, width = 8, height = 5, dpi = 300)
  cat("\nFigure saved: figure_gam_depth_shape.png\n")
  
  # --- Breakpoint estimation, only if non-linearity was found ---
  if (edf > 2 && delta_aic > 2 && requireNamespace("segmented", quietly = TRUE)) {
    library(segmented)
    
    cat("\n── Breakpoint estimation (segmented regression) ──\n\n")
    
    seg_model <- tryCatch(
      segmented(model_linear, seg.Z = ~log_depth, npsi = 2),
      error = function(e) {
        cat("2-breakpoint model failed to converge, trying 1 breakpoint\n")
        segmented(model_linear, seg.Z = ~log_depth, npsi = 1)
      }
    )
    print(summary(seg_model))
    
    breakpoints_log <- seg_model$psi[, "Est."]
    breakpoints_depth <- round(10^breakpoints_log - 1, 1)
    breakpoint_ci <- confint(seg_model)
    
    cat("\nRaw breakpoint CI matrix (check column order matches Est./lower/upper):\n")
    print(breakpoint_ci)
    
    cat("\n── Estimated breakpoint(s), back-transformed to metres ──\n")
    breakpoint_table <- tibble()
    for (i in seq_along(breakpoints_depth)) {
      ci_log <- breakpoint_ci[i, ]
      ci_depth_lo <- round(10^ci_log[2] - 1, 1)
      ci_depth_hi <- round(10^ci_log[3] - 1, 1)
      cat(sprintf("Breakpoint %d: %.1f m (95%% CI: %.1f - %.1f m)\n",
                  i, breakpoints_depth[i], ci_depth_lo, ci_depth_hi))
      breakpoint_table <- bind_rows(breakpoint_table, tibble(
        breakpoint_number = i, depth_m = breakpoints_depth[i],
        ci_lower_m = ci_depth_lo, ci_upper_m = ci_depth_hi
      ))
    }
    write_csv(breakpoint_table, "depth_breakpoints.csv")
    cat("\nCompare to the ecological zone boundaries (200m, 1000m, 4000m).\n")
    
    # --- Check whether a 3rd breakpoint further out improves the fit ---
    # (the 2-breakpoint model forces everything past breakpoint 2 into a
    # single segment - if the steep 85m+ decline actually flattens out
    # again at greater depth, e.g. near 1000m or 4000m, a 2-breakpoint
    # model can't capture that. This checks explicitly.)
    cat("\n── Checking whether a 3rd breakpoint improves the fit ──\n\n")
    
    seg_model_3 <- tryCatch(
      segmented(model_linear, seg.Z = ~log_depth, npsi = 3),
      error = function(e) {
        cat("3-breakpoint model failed to converge:", conditionMessage(e), "\n")
        NULL
      }
    )
    
    if (!is.null(seg_model_3)) {
      aic_2 <- AIC(seg_model)
      aic_3 <- AIC(seg_model_3)
      cat("2-breakpoint model AIC:", round(aic_2, 1), "\n")
      cat("3-breakpoint model AIC:", round(aic_3, 1), "\n")
      
      if (aic_3 < aic_2 - 2) {
        cat("\nAIC suggests the 3-breakpoint model fits better (\u0394AIC =",
            round(aic_2 - aic_3, 1), ") - but AIC alone is not enough here.\n")
        cat("Checking for degeneracy (breakpoints collapsing onto each other,\n")
        cat("which produces a spurious AIC improvement from local noise-fitting\n")
        cat("rather than a genuine third structural change):\n\n")
        
        breakpoints_log_3 <- seg_model_3$psi[, "Est."]
        
        # Degeneracy check 1: are any two breakpoints suspiciously close
        # together relative to the full depth range?
        min_gap <- min(diff(sort(breakpoints_log_3)))
        full_range <- diff(range(model_data$log_depth))
        gap_is_tiny <- min_gap < 0.05 * full_range
        
        # Degeneracy check 2: are the jump-size (U.) coefficients wildly
        # unstable (SE many times larger than the estimate)?
        seg_coefs <- summary(seg_model_3)$coefficients
        u_rows <- grepl("^U[0-9]+\\.log_depth$", rownames(seg_coefs))
        u_unstable <- any(abs(seg_coefs[u_rows, "Std. Error"]) >
                            10 * abs(seg_coefs[u_rows, "Estimate"]))
        
        if (gap_is_tiny || u_unstable) {
          cat("\u26A0\uFE0F  DEGENERACY DETECTED - the 3-breakpoint model is NOT reliable:\n")
          if (gap_is_tiny) {
            cat("  - Two breakpoints have collapsed onto nearly the same location\n")
            cat("    (smallest gap:", round(min_gap, 3), "on the log-depth scale,\n")
            cat("    versus a full data range of", round(full_range, 2), ").\n")
          }
          if (u_unstable) {
            cat("  - At least one jump-size coefficient has a standard error\n")
            cat("    more than 10x its estimate - the fit is numerically unstable.\n")
          }
          cat("\nREJECTING the 3-breakpoint model despite its lower AIC - this is\n")
          cat("overfitting to local noise, not a genuine third regime change.\n")
          cat("STICK WITH the 2-breakpoint result saved above (depth_breakpoints.csv)\n")
          cat("as the final answer.\n")
        } else {
          cat("No degeneracy detected - the 3-breakpoint model appears genuine.\n\n")
          
          print(summary(seg_model_3))
          
          breakpoints_depth_3 <- round(10^breakpoints_log_3 - 1, 1)
          breakpoint_ci_3 <- confint(seg_model_3)
          
          cat("\n── 3-breakpoint estimates, back-transformed to metres ──\n")
          breakpoint_table_3 <- tibble()
          for (i in seq_along(breakpoints_depth_3)) {
            ci_log <- breakpoint_ci_3[i, ]
            ci_depth_lo <- round(10^ci_log[2] - 1, 1)
            ci_depth_hi <- round(10^ci_log[3] - 1, 1)
            cat(sprintf("Breakpoint %d: %.1f m (95%% CI: %.1f - %.1f m)\n",
                        i, breakpoints_depth_3[i], ci_depth_lo, ci_depth_hi))
            breakpoint_table_3 <- bind_rows(breakpoint_table_3, tibble(
              breakpoint_number = i, depth_m = breakpoints_depth_3[i],
              ci_lower_m = ci_depth_lo, ci_upper_m = ci_depth_hi
            ))
          }
          write_csv(breakpoint_table_3, "depth_breakpoints_3psi.csv")
          cat("\nUse THESE 3 breakpoints instead of the 2-breakpoint set above -\n")
          cat("saved to depth_breakpoints_3psi.csv.\n")
        }
        
      } else {
        cat("\nRESULT: the 3-breakpoint model does NOT fit meaningfully better\n")
        cat("(\u0394AIC =", round(aic_2 - aic_3, 1), ", below the >2 threshold) - the\n")
        cat("2-breakpoint model (already saved above) remains the best choice.\n")
      }
    }
  } else if (!requireNamespace("segmented", quietly = TRUE)) {
    cat("\nPackage 'segmented' not installed - run install.packages('segmented')\n")
    cat("to estimate the breakpoint location precisely.\n")
  }
  
} else {
  cat("Package 'mgcv' not installed - run install.packages('mgcv').\n")
}


# =============================================================================
# SECTION 3b: DOES THE DEPTH RELATIONSHIP'S SHAPE DIFFER BY CLASS?
# The Section 3 GAM/breakpoint above used a SHARED-SHAPE model - echino_class
# entered as a simple additive term, so every class was forced to share one
# curve shape (same breakpoints), just shifted up/down. This checks that
# assumption directly, given the interaction test (see TIa_full.R) already
# found classes decline at significantly different RATES with depth.
# =============================================================================

cat("\n\n=== SECTION 3b: Class-specific depth relationship shape ===\n\n")

cat("── Sample size per class x depth zone ──\n")
print(table(model_data$echino_class, model_data$depth_zone))
cat("\nNote: Crinoidea and Echinoidea have ZERO records at Abyssal-\n")
cat("their per-class models below may fail or be unstable for this reason.\n\n")

if (requireNamespace("mgcv", quietly = TRUE)) {
  
  # k = 4 caps each class's smooth at a modest effective degrees of freedom,
  # since sparse classes cannot support a highly flexible curve without
  # overfitting to noise (an unconstrained version produced visibly
  # implausible, wildly oscillating curves for the sparsest classes).
  model_gam_byclass <- gam(
    species_level ~ echino_class + s(log_depth, by = echino_class, k = 4) + obs_quality_flag,
    data = model_data, family = binomial
  )
  
  cat("── Per-class smooth term summary (k=4 constrained) ──\n\n")
  print(summary(model_gam_byclass))
  
  aic_shared_gam <- AIC(model_gam)
  aic_byclass_gam <- AIC(model_gam_byclass)
  cat("\nShared-shape GAM AIC:  ", round(aic_shared_gam, 1), "\n")
  cat("Class-specific GAM AIC:", round(aic_byclass_gam, 1), "\n")
  cat("\u0394AIC:", round(aic_shared_gam - aic_byclass_gam, 1), "in favour of class-specific model\n")
  cat("\nNOTE: before trusting this AIC comparison, inspect\n")
  cat("figure_gam_depth_shape_by_class.png for any class whose curve still\n")
  cat("oscillates implausibly even at k=4 - that class's AIC contribution is\n")
  cat("likely inflated by overfitting to sparse/heterogeneous data, not\n")
  cat("evidence of a genuine shape difference. Cross-check any such class\n")
  cat("against the per-class breakpoint table below, which uses a much\n")
  cat("simpler, less overfit-prone model.\n")
  
  depth_range_bc <- seq(min(model_data$log_depth), max(model_data$log_depth), length.out = 200)
  pred_grid <- expand_grid(
    log_depth = depth_range_bc, echino_class = levels(model_data$echino_class)
  ) %>%
    mutate(
      echino_class = factor(echino_class, levels = levels(model_data$echino_class)),
      obs_quality_flag = factor("1", levels = levels(model_data$obs_quality_flag))
    )
  pred_bc <- predict(model_gam_byclass, newdata = pred_grid, type = "response", se.fit = TRUE)
  pred_grid <- pred_grid %>%
    mutate(
      depth_m = 10^log_depth - 1, fit = pred_bc$fit,
      lower = pmax(0, pred_bc$fit - 1.96 * pred_bc$se.fit),
      upper = pmin(1, pred_bc$fit + 1.96 * pred_bc$se.fit)
    ) %>%
    filter(is.finite(depth_m), depth_m > 0)
  
  p_byclass <- pred_grid %>%
    ggplot(aes(x = depth_m, y = fit * 100, color = echino_class, fill = echino_class)) +
    geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.9) +
    scale_x_log10(labels = scales::comma) +
    scale_color_manual(values = c(
      "Asteroidea" = "#2a78d6", "Crinoidea" = "#1baf7a", "Echinoidea" = "#eda100",
      "Holothuroidea" = "#4a3aa7", "Ophiuroidea" = "#e34948"
    )) +
    scale_fill_manual(values = c(
      "Asteroidea" = "#2a78d6", "Crinoidea" = "#1baf7a", "Echinoidea" = "#eda100",
      "Holothuroidea" = "#4a3aa7", "Ophiuroidea" = "#e34948"
    )) +
    labs(x = "Depth (m, log scale)", y = "Predicted probability of species-level ID (%)",
         color = "Class", fill = "Class",
         title = "Class-specific GAM-fitted depth relationships (Flag 1 reference)") +
    theme_minimal(base_size = 12)
  
  ggsave("figure_gam_depth_shape_by_class.png", p_byclass, width = 9, height = 5.5, dpi = 300)
  cat("\nFigure saved: figure_gam_depth_shape_by_class.png\n")
  
  # --- Per-class breakpoint estimation (simpler, single-breakpoint model) ---
  if (requireNamespace("segmented", quietly = TRUE)) {
    cat("\n── Per-class breakpoint estimation ──\n\n")
    
    classes <- levels(model_data$echino_class)
    breakpoint_by_class <- tibble()
    
    for (cls in classes) {
      cat("── ", cls, " ──\n")
      d <- model_data %>% filter(echino_class == cls)
      
      m_linear_cls <- tryCatch(
        glm(species_level ~ log_depth + obs_quality_flag, data = d, family = binomial),
        error = function(e) NULL
      )
      if (is.null(m_linear_cls)) { cat("  Linear model failed to fit - skipping.\n\n"); next }
      
      seg_m <- tryCatch(
        segmented(m_linear_cls, seg.Z = ~log_depth, npsi = 1),
        error = function(e) { cat("  Segmented model failed to converge.\n"); NULL }
      )
      if (is.null(seg_m)) {
        cat("  No stable breakpoint found (likely sparse deep records).\n\n")
        next
      }
      
      bp_log <- seg_m$psi[, "Est."]
      bp_depth <- round(10^bp_log - 1, 1)
      ci_cls <- tryCatch(confint(seg_m), error = function(e) NULL)
      
      seg_coefs_cls <- summary(seg_m)$coefficients
      u_row_cls <- grepl("^U1\\.log_depth$", rownames(seg_coefs_cls))
      u_unstable_cls <- any(abs(seg_coefs_cls[u_row_cls, "Std. Error"]) >
                              10 * abs(seg_coefs_cls[u_row_cls, "Estimate"]))
      
      if (u_unstable_cls) {
        cat("  \u26A0\uFE0F Unstable jump coefficient - NOT reporting for", cls, ".\n\n")
        next
      }
      
      if (!is.null(ci_cls)) {
        ci_lo <- round(10^ci_cls[1, 2] - 1, 1)
        ci_hi <- round(10^ci_cls[1, 3] - 1, 1)
        cat(sprintf("  Breakpoint: %.1f m (95%% CI: %.1f - %.1f m)\n\n", bp_depth, ci_lo, ci_hi))
        breakpoint_by_class <- bind_rows(breakpoint_by_class, tibble(
          echino_class = cls, n = nrow(d), breakpoint_m = bp_depth,
          ci_lower_m = ci_lo, ci_upper_m = ci_hi
        ))
      }
    }
    
    cat("── Summary: per-class breakpoints ──\n")
    if (nrow(breakpoint_by_class) > 0) {
      print(breakpoint_by_class, n = Inf)
      write_csv(breakpoint_by_class, "depth_breakpoints_by_class.csv")
      cat("\nCompare to the pooled estimate (22.5m, 84.5m). Overlapping CIs\n")
      cat("across classes indicate the LOCATION of the collapse is broadly\n")
      cat("consistent even though the earlier interaction test showed the\n")
      cat("RATE of decline differs significantly by class.\n")
    } else {
      cat("No class produced a stable, trustworthy breakpoint estimate.\n")
    }
  } else {
    cat("Package 'segmented' not installed - run install.packages('segmented').\n")
  }
  
} else {
  cat("Package 'mgcv' not installed - skipping class-specific shape check.\n")
}


# =============================================================================
# SECTION 4: ORDINAL LOGISTIC REGRESSION (full resolution hierarchy)
# Uses the log-depth specification justified in Section 3.
# =============================================================================

cat("\n\n=== SECTION 4: Ordinal logistic regression ===\n\n")

model_data_ord <- model_data %>%
  filter(rank6 %in% c("Species", "Genus", "Family", "Order", "Class")) %>%
  mutate(rank_ordered = factor(rank6, levels = c("Species", "Genus", "Family", "Order", "Class"),
                               ordered = TRUE))

cat("rank6 distribution in this subset (Phylum excluded - zero records with depth):\n")
print(table(model_data_ord$rank6))
cat("\nRecords in ordinal model:", nrow(model_data_ord), "\n\n")

ordinal_model <- MASS::polr(
  rank_ordered ~ log_depth + echino_class + obs_quality_flag,
  data = model_data_ord, Hess = TRUE
)
print(summary(ordinal_model))

coef_table <- coef(summary(ordinal_model))
p_values <- pnorm(abs(coef_table[, "t value"]), lower.tail = FALSE) * 2
coef_table <- cbind(coef_table, p_value = p_values)
cat("\n── Coefficients with p-values ──\n")
print(round(coef_table, 4))

or_ci <- exp(cbind(OR = coef(ordinal_model), confint(ordinal_model)))
cat("\n── Odds ratios (95% CI) - odds of a WORSE resolution category ──\n")
print(round(or_ci, 3))

ordinal_results <- as_tibble(round(or_ci, 3), rownames = "term") %>%
  rename(ci_lower = `2.5 %`, ci_upper = `97.5 %`) %>%
  left_join(
    as_tibble(round(coef_table, 4), rownames = "term") %>% dplyr::select(term, p_value),
    by = "term"
  )
write_csv(ordinal_results, "ordinal_regression_taxonomic_completeness.csv")


# =============================================================================
# SECTION 5: PROPORTIONAL ODDS CHECK + cutpoint-specific models
# (the corrected result, given a violation is expected)
# =============================================================================

cat("\n\n=== SECTION 5: Proportional odds assumption check ===\n\n")
cat("H0: each predictor's effect is the SAME at every cutpoint\n\n")

if (requireNamespace("brant", quietly = TRUE)) {
  brant_result <- brant::brant(ordinal_model)
  print(brant_result)
  cat("\nA significant (p < .05) result means that predictor's effect genuinely\n")
  cat("differs across cutpoints - the single ordinal OR is an oversimplification\n")
  cat("for that term. See cutpoint-specific models below for the corrected view.\n")
} else {
  cat("Package 'brant' not installed - run install.packages('brant').\n")
  cat("Proceeding directly to cutpoint-specific models below, which are\n")
  cat("informative regardless of whether the formal test is run.\n")
}

cat("\n\n── Cutpoint-specific binary models ──\n")
cat("(report these INSTEAD of, or alongside, the single ordinal OR if\n")
cat("proportional odds is violated - depth's effect size at each cutpoint\n")
cat("is itself the finding, not noise to average away)\n\n")

cutpoint_defs <- list(
  list(better = c("Species"), label = "Species | worse"),
  list(better = c("Species", "Genus"), label = "\u2264Genus | worse"),
  list(better = c("Species", "Genus", "Family"), label = "\u2264Family | worse"),
  list(better = c("Species", "Genus", "Family", "Order"), label = "\u2264Order | Class")
)

cutpoint_results <- map_dfr(cutpoint_defs, function(cp) {
  d <- model_data_ord %>% mutate(y = as.integer(rank_ordered %in% cp$better))
  m <- glm(y ~ log_depth + echino_class + obs_quality_flag, data = d, family = binomial)
  ci <- confint(m, "log_depth", trace = FALSE)
  tibble(
    cutpoint = cp$label,
    n = nrow(d),
    OR_worse_resolution = round(exp(-coef(m)["log_depth"]), 3),
    ci_lower = round(exp(-ci[2]), 3),
    ci_upper = round(exp(-ci[1]), 3),
    p_value = signif(summary(m)$coefficients["log_depth", "Pr(>|z|)"], 3)
  )
})

cat("── OR (per 10-fold depth increase) of being on the WORSE side of each cutpoint ──\n")
print(cutpoint_results, n = Inf)
write_csv(cutpoint_results, "cutpoint_specific_odds_ratios.csv")

cat("\nIf these ORs increase from top to bottom, depth has a progressively\n")
cat("STRONGER effect on the worst resolution outcomes - the concrete finding\n")
cat("behind the proportional-odds violation.\n")

cat("\n\n=== Taxonomic completeness analysis complete ===\n")
cat("Outputs: table1_rank_by_depthbin.csv, figure_rank_composition_by_depth.png,\n")
cat("figure_gam_depth_shape.png, depth_breakpoints.csv,\n")
cat("figure_gam_depth_shape_by_class.png, depth_breakpoints_by_class.csv,\n")
cat("ordinal_regression_taxonomic_completeness.csv, cutpoint_specific_odds_ratios.csv\n")