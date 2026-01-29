# =====================================================================
# 10/12/2025, C. A. Higham
# Linear mixed-effects model for normalised deposition rate (Ddot)
#
# This script implements the method described in Section 2.5.1:
# - Build normalised deposition rate Ddot 
# - Log-transform Ddot
# - Assess normality and homoscedasticity (Shapiro–Wilk + Q–Q plots)
# - Fit linear mixed-effects models with ventilation and interval
# - Use ML-based likelihood ratio tests for model simplification
# - Refit final reduced model with REML
# - Use Type III Wald tests for fixed-effects significance
# - Check residual and random-effect assumptions
# =====================================================================

# =====================================================================
# Packages and theme
# =====================================================================
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(tidyr)
library(purrr)
library(readr)

# Base plotting theme used for all diagnostic plots
theme_ampas <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

# =====================================================================
# Paths
# =====================================================================

# Base project directory
base_dir <- "Your_Directory"

# ICP-MS data directory
data_dir <- file.path(base_dir, "results_csv", "deposition_data")

# Save plot data directory
save_dir <- file.path(base_dir, "output_plots")
# ----------------------------
# Plot save helpers
# ----------------------------
# Narrow, tall plots (e.g. Q–Q plots for each ACH)
save_plot_narrow <- function(plot, name) {
  ggsave(
    filename = file.path(save_dir, paste0(name, ".png")),
    plot = plot,
    width = 8.5, height = 10, units = "cm", dpi = 400,
    bg = "transparent"
  )
}

# Narrow, shorter plots (e.g. residuals vs fitted)
save_plot_narrow_short <- function(plot, name) {
  ggsave(
    filename = file.path(save_dir, paste0(name, ".png")),
    plot = plot,
    width = 8.5, height = 8, units = "cm", dpi = 400,
    bg = "transparent"
  )
}

# Wide, tall plots
save_plot_wide <- function(plot, name) {
  ggsave(
    filename = file.path(save_dir, paste0(name, ".png")),
    plot = plot,
    width = 18.5, height = 10, units = "cm", dpi = 400,
    bg = "transparent"
  )
}

# Wide, short plots (e.g. combined Q–Q for residuals + random effects)
save_plot_wide_short <- function(plot, name) {
  ggsave(
    filename = file.path(save_dir, paste0(name, ".png")),
    plot = plot,
    width = 18.5, height = 6, units = "cm", dpi = 400,
    bg = "transparent"
  )
}

# Constants: plate geometry and experimental settings
plate_radius_cm <- 5.5 / 2
plate_area_m2   <- pi * (plate_radius_cm^2) / 1e4  # convert cm^2 to m^2
extract_vol_L   <- 0.006                           # extraction volume (L)
vent_rates      <- c(1.5, 10)                      # ACH values used in model

# =====================================================================
# Build normalised deposition rate per hour (Ddot = norm_dep_rate)
# =====================================================================

# Load the AMPAS plate data and impinger airborne data
ampas_55mm    <- read_csv(file.path(data_dir, "ampas-55mm.csv"),    show_col_types = FALSE)
impinger_conc <- read_csv(file.path(data_dir, "impinger_conc.csv"), show_col_types = FALSE)

# Compute airborne LiCl concentration per hour from impinger
#     (control_ppb_per_h is the airborne reference used for normalisation)
impinger_summary_all <- impinger_conc %>%
  filter(
    Ventilation_Rate %in% vent_rates,
    is.finite(conc_ppb),
    is.finite(sample_time),
    sample_time > 0
  ) %>%
  mutate(
    control_ppb_per_h = conc_ppb / sample_time * 60
  ) %>%
  select(Experiment_Index, Ventilation_Rate, control_ppb_per_h)

# For each ventilation rate, find the maximum airborne concentration per hour.
# This is used to re-scale the plate deposition so that Ddot is scaled by max impinger conc
impinger_max_by_vr <- impinger_summary_all %>%
  group_by(Ventilation_Rate) %>%
  summarise(
    max_control_ppb_per_h = max(control_ppb_per_h, na.rm = TRUE),
    .groups = "drop"
  )

# Convert plate concentrations (ppb in extract) to deposition rate (µg m^-2 h^-1)
ampas_dep_rate <- ampas_55mm %>%
  filter(
    Ventilation_Rate %in% vent_rates,
    AMPAS_Interval %in% c(30, 60),   # only 30 and 60 min intervals used here
    is.finite(conc_ppb)
  ) %>%
  mutate(
    mass_ug     = conc_ppb * extract_vol_L,          # µg collected on plate
    dep_ug_m2   = mass_ug / plate_area_m2,          # µg m^-2
    dep_ug_m2_h = dep_ug_m2 / (AMPAS_Interval / 60) # convert to µg m^-2 h^-1
  ) %>%
  select(
    Experiment_Index, Ventilation_Rate, AMPAS_Interval,
    Location, dep_ug_m2_h
  )

# Build the final normalised deposition rate dataset:
# This corresponds to Ddot in the method section.
norm_dep_rate <- ampas_dep_rate %>%
  left_join(impinger_summary_all,
            by = c("Experiment_Index", "Ventilation_Rate")) %>%
  left_join(impinger_max_by_vr, by = "Ventilation_Rate") %>%
  filter(
    is.finite(dep_ug_m2_h),
    !is.na(control_ppb_per_h),
    control_ppb_per_h > 0
  ) %>%
  mutate(
    dep_norm_raw  = dep_ug_m2_h / control_ppb_per_h,
    norm_dep_rate = dep_norm_raw * max_control_ppb_per_h,   # µg m^-2 h^-1
    Ventilation_f = factor(Ventilation_Rate, levels = vent_rates),
    Interval_f    = factor(AMPAS_Interval, levels = c(30, 60)),
    Location_f    = factor(Location),
    Experiment_f  = factor(Experiment_Index)
  )

# =====================================================================
# Define Ddot and assess normality (raw vs log)
# This corresponds to checking model assumptions before fitting LMM.
# =====================================================================

# Define Ddot and log(Ddot) for modelling, and remove non-positive values
log_norm_dep_rate <- norm_dep_rate %>%
  mutate(
    Ddot     = norm_dep_rate,       # Ddot in the method (µg m^-2 h^-1)
    log_Ddot = log(Ddot)            # response used in the mixed model
  ) %>%
  filter(is.finite(Ddot), Ddot > 0)

# Shapiro–Wilk tests on raw Ddot by ACH × interval
# (these show the right-skew / non-normality motivating the log transform)
shapiro_raw_by_cell <- log_norm_dep_rate %>%
  group_by(Ventilation_f, Interval_f) %>%
  summarise(
    n      = length(Ddot),
    W      = shapiro.test(Ddot)$statistic,
    pvalue = shapiro.test(Ddot)$p.value,
    .groups = "drop"
  )

cat("\n===== Shapiro–Wilk tests for raw Ddot (by ACH × interval) =====\n")
print(shapiro_raw_by_cell)

# Shapiro–Wilk tests on log(Ddot) by ACH × interval
# (these support the use of log-transformed Ddot in the LMM)
shapiro_log_by_cell <- log_norm_dep_rate %>%
  group_by(Ventilation_f, Interval_f) %>%
  summarise(
    n      = length(log_Ddot),
    W      = shapiro.test(log_Ddot)$statistic,
    pvalue = shapiro.test(log_Ddot)$p.value,
    .groups = "drop"
  )

cat("\n===== Shapiro–Wilk tests for log(Ddot) (by ACH × interval) =====\n")
print(shapiro_log_by_cell)

# =====================================================================
# Q–Q plots of raw vs log Ddot (by ACH)
# Visual assessment of normality, as described in the methods.
# =====================================================================

# Long-format data for Q–Q plots (raw vs log on separate facets)
qq_data <- log_norm_dep_rate %>%
  select(Ventilation_f, Ddot, log_Ddot) %>%
  pivot_longer(
    cols      = c(Ddot, log_Ddot),
    names_to  = "Scale",
    values_to = "Value"
  ) %>%
  mutate(
    Scale = dplyr::recode(
      Scale,
      Ddot     = "Raw D\u0307",
      log_Ddot = "log(D\u0307)"
    )
  )

# Colour by ventilation rate
ach_cols <- c("1.5" = "#009aDE", "10" = "#ff1f5b")

# Q–Q plots for 1.5 ACH
qq_data_15 <- qq_data %>%
  filter(Ventilation_f == "1.5")

p_qq_15 <- ggplot(
  qq_data_15,
  aes(sample = Value, colour = Ventilation_f)
) +
  stat_qq(size = 1.4) +
  stat_qq_line(colour = "black") +
  facet_grid(Scale ~ ., scales = "free_y") +
  scale_colour_manual(
    values = ach_cols,
    name = "Ventilation rate\n(ACH)"
  ) +
  labs(
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  theme_ampas +
  theme(
    legend.position      = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", colour = "grey70"),
    legend.title         = element_text(size = 10, hjust = 0.5),
    legend.text          = element_text(size = 10),
    legend.key.size      = unit(0.7, "lines")
  )

print(p_qq_15)
save_plot_narrow(p_qq_15, "qq_plot_15")

# Q–Q plots for 10 ACH
qq_data_10 <- qq_data %>%
  filter(Ventilation_f == "10")

p_qq_10 <- ggplot(
  qq_data_10,
  aes(sample = Value, colour = Ventilation_f)
) +
  stat_qq(size = 1.4) +
  stat_qq_line(colour = "black") +
  facet_grid(Scale ~ ., scales = "free_y") +
  scale_colour_manual(
    values = ach_cols,
    name = "Ventilation rate\n(ACH)"
  ) +
  labs(
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  theme_ampas +
  theme(
    legend.position      = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", colour = "grey70"),
    legend.title         = element_text(size = 10, hjust = 0.5),
    legend.text          = element_text(size = 10),
    legend.key.size      = unit(0.7, "lines")
  )

print(p_qq_10)
save_plot_narrow(p_qq_10, "qq_plot_10")

# =====================================================================
# Linear mixed-effects models: raw vs log(Ddot)
# Full model with ventilation, interval and their interaction.
# =====================================================================

# Full model on raw Ddot (used mainly to illustrate poor fit / skewness)
model_raw <- lmer(
  Ddot ~ Ventilation_f * Interval_f +
    (1 | Location_f) + (1 | Experiment_f),
  data = log_norm_dep_rate
)

# Full model on log Ddot (this is the model structure described in the methods)
model_log <- lmer(
  log_Ddot ~ Ventilation_f * Interval_f +
    (1 | Location_f) + (1 | Experiment_f),
  data = log_norm_dep_rate
)

cat("\n===== Mixed model on raw Ddot =====\n")
print(summary(model_raw))

cat("\n===== Mixed model on log(Ddot) =====\n")
print(summary(model_log))

# From here on, all inference is based on the log-transformed model.

# =====================================================================
# Model selection and inference (log-transformed Ddot)
# Maximum Likelihood Likelehood Ratio Test for simplification
# Then REML.
# =====================================================================

# FULL model (REML) — used for Type III ANOVA (Wald tests)
model_log_full <- lmer(
  log_Ddot ~ Ventilation_f * Interval_f +
    (1 | Location_f) + (1 | Experiment_f),
  data = log_norm_dep_rate,
  REML = TRUE
)

cat("\n===== Type III ANOVA for FULL model (REML) =====\n")
anova_full <- anova(model_log_full, type = 3)
print(anova_full)

# FULL vs REDUCED (Maximum Likelihood) — likelihood ratio test for model simplification
# This checks whether interval and interaction terms can be removed
# without a significant loss of fit.
model_log_full_ML <- lmer(
  log_Ddot ~ Ventilation_f * Interval_f +
    (1 | Location_f) + (1 | Experiment_f),
  data = log_norm_dep_rate,
  REML = FALSE
)

model_log_reduced_ML <- lmer(
  log_Ddot ~ Ventilation_f +
    (1 | Location_f) + (1 | Experiment_f),
  data = log_norm_dep_rate,
  REML = FALSE
)

cat("\n===== Likelihood ratio test (FULL vs REDUCED) — both ML =====\n")
lrt_full_reduced <- anova(model_log_reduced_ML, model_log_full_ML)
print(lrt_full_reduced)

# FINAL reduced model (ventilation only, REML) — reported model
# This is the model used for effect size estimation and interpretation.
model_log_reduced_REML <- lmer(
  log_Ddot ~ Ventilation_f +
    (1 | Location_f) + (1 | Experiment_f),
  data = log_norm_dep_rate,
  REML = TRUE
)

cat("\n===== FINAL reduced model (ventilation only, REML) =====\n")
print(summary(model_log_reduced_REML))

# Type III ANOVA on the final reduced model (ventilation only)
# This provides Wald tests for the significance of ventilation.
cat("\n===== Type III ANOVA for FINAL reduced model (REML) =====\n")
anova_reduced <- anova(model_log_reduced_REML, type = 3)
print(anova_reduced)

# =====================================================================
# Diagnostics for the final log model
# Residuals + random effects: normality and homoscedasticity checks.
# =====================================================================

# Extract residuals and fitted values from the final model
resid_final  <- resid(model_log_reduced_REML)
fitted_final <- fitted(model_log_reduced_REML)

# Extract random intercepts for experiments and locations
re_exp <- ranef(model_log_reduced_REML)$Experiment_f[, "(Intercept)"]
re_loc <- ranef(model_log_reduced_REML)$Location_f[, "(Intercept)"]

# Combined Q–Q plot for residuals and random effects
# This visually assesses the normality assumptions for:
# - residuals,
# - experiment-level random intercepts,
# - location-level random intercepts.
re_df <- tibble(
  Residual = c(resid_final, re_exp, re_loc),
  Type = c(
    rep("Residuals", length(resid_final)),
    rep("Experiment random effects", length(re_exp)),
    rep("Location random effects", length(re_loc))
  )
)

re_df$Type <- factor(
  re_df$Type,
  levels = c("Residuals", "Experiment random effects", "Location random effects")
)

p_resid_qq_combined <- ggplot(re_df, aes(sample = Residual)) +
  stat_qq(size = 1) +
  stat_qq_line() +
  facet_wrap(~ Type, scales = "free") +
  labs(
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  theme_ampas

print(p_resid_qq_combined)
save_plot_wide_short(p_resid_qq_combined, "qq-residuals-random")

# Shapiro–Wilk tests for residuals and random effects
# These provide a numerical check of the normality assumption.
shapiro_resid <- shapiro.test(resid_final)
shapiro_exp   <- shapiro.test(re_exp)
shapiro_loc   <- shapiro.test(re_loc)

cat("\n===== Shapiro–Wilk: Residuals =====\n")
print(shapiro_resid)

cat("\n===== Shapiro–Wilk: Experiment random effects =====\n")
print(shapiro_exp)

cat("\n===== Shapiro–Wilk: Location random effects =====\n")
print(shapiro_loc)

# Homoscedasticity check (residuals vs fitted values)
log_norm_dep_rate <- log_norm_dep_rate %>%
  mutate(
    resid_final = resid_final,
    fit_final   = fitted_final
  )

p_homo_fitted <- ggplot(log_norm_dep_rate, aes(x = fit_final, y = resid_final)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Fitted values",
    y = "Residuals"
  ) +
  theme_ampas

print(p_homo_fitted)
save_plot_narrow_short(p_homo_fitted, "homo_resid_vs_fitted_final")
