# =====================================================================
# 10/12/2025, C. A. Higham, LiCl deposition output code
# Processes lithium chloride (LiCl) tracer data collected using the
# AMPAS device and impinger samples in a controlled aerobiology
# chamber. It performs normalisation, summarises spatial and temporal
# variability, and generates figures.
# =====================================================================

# =====================================================================
# Load libraries
# =====================================================================
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)
library(grid)

# =====================================================================
# Folders for data & output
# =====================================================================

# Base project directory
base_dir <- "Your_Directory"

# ICP-MS data directory
data_dir <- file.path(base_dir, "results_csv", "deposition_data")

# Save plot data directory
save_dir <- file.path(base_dir, "output_plots")

# =====================================================================
# Plot saving helpers
# =====================================================================
save_plot_narrow <- function(plot, name) {
  ggsave(
    filename = file.path(save_dir, paste0(name, ".png")),
    plot     = plot,
    width    = 8.5,
    height   = 10,
    units    = "cm",
    dpi      = 400,
    bg       = "transparent"
  )
}

save_plot_wide <- function(plot, name) {
  ggsave(
    filename = file.path(save_dir, paste0(name, ".png")),
    plot     = plot,
    width    = 18.5,
    height   = 10,
    units    = "cm",
    dpi      = 400,
    bg       = "transparent"
  )
}

# =====================================================================
# Base theme
# =====================================================================
theme_ampas <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

# =====================================================================
# Load data
# =====================================================================

# AMPAS plate data
ampas_55mm            <- read_csv(file.path(data_dir, "ampas-55mm.csv"),            show_col_types = FALSE)
# Impinger data
impinger_conc         <- read_csv(file.path(data_dir, "impinger_conc.csv"),         show_col_types = FALSE)
# 90 mm vs 55 mm plate data
plate_area_comparison <- read_csv(file.path(data_dir, "plate_area_comparison.csv"), show_col_types = FALSE)

# =====================================================================
# Constants and helper objects
# =====================================================================

# Ventilation rates
vent_rates <- c(1.5, 10)

# Plate area ratio (constant)
area_ratio <- (90^2) / (55^2)

# Position dodge
pd <- position_dodge(width = 0.5)

# Colours / shapes for ACH
scale_ach_colour <- scale_colour_manual(values = c("10" = "#ff1f5b", "1.5" = "#009aDE"))
scale_ach_shape  <- scale_shape_manual(values = c("10" = 15, "1.5" = 16))

ach_legend_tweaks <- list(
  theme(legend.position = "top", legend.title = element_text(hjust = 0.5)),
  guides(
    colour = guide_legend(title.position = "top"),
    shape  = guide_legend(title.position = "top")
  )
)

# ACH aesthetics
aes_ach_location <- function() {
  aes(
    x = factor(Location),
    colour = Ventilation_f,
    shape  = Ventilation_f,
    group  = Ventilation_f
  )
}

# Ratio error propagation
ratio_with_sd <- function(mean_num, sd_num, mean_den, sd_den) {
  ratio <- mean_num / mean_den
  ratio * sqrt(
    (sd_num / mean_num)^2 +
      (sd_den / mean_den)^2
  )
}

# =====================================================================
# Normalise concentrations by the impinger concentration
# Multiply them by (C_{max} / C_{exp})
# See Section 2.5
# =====================================================================

ampas_subset_all <- ampas_55mm %>%
  filter(Ventilation_Rate %in% vent_rates, is.finite(conc_ppb)) %>%
  select(
    Experiment_Index, Type, Ventilation_Rate, AMPAS_Interval,
    Location, Consecutive_Sample_ID, conc_ppb
  )

impinger_summary_all <- impinger_conc %>%
  filter(
    Ventilation_Rate %in% vent_rates,
    is.finite(conc_ppb),
    is.finite(sample_time),
    sample_time > 0
  ) %>%
  mutate(control_ppb_per_h = conc_ppb / sample_time * 60) %>%
  select(Experiment_Index, Ventilation_Rate, control_ppb_per_h)

impinger_max_by_vr <- impinger_summary_all %>%
  group_by(Ventilation_Rate) %>%
  summarise(max_control_ppb_per_h = max(control_ppb_per_h), .groups = "drop")

ampas_norm_all <- ampas_subset_all %>%
  left_join(impinger_summary_all, by = c("Experiment_Index", "Ventilation_Rate")) %>%
  left_join(impinger_max_by_vr, by = "Ventilation_Rate") %>%
  filter(control_ppb_per_h > 0) %>%
  mutate(
    conc_norm_imp_raw = conc_ppb / control_ppb_per_h,
    conc_norm_imp     = conc_norm_imp_raw * max_control_ppb_per_h
  )

# =====================================================================
# Coefficient of Variability summary
# =====================================================================

ampas_variability <- ampas_norm_all %>%
  group_by(Ventilation_Rate, Experiment_Index, AMPAS_Interval, Location, Type) %>%
  summarise(
    n              = n(),
    mean_conc_norm = mean(conc_norm_imp),
    sd_conc_norm   = sd(conc_norm_imp),
    .groups = "drop"
  ) %>%
  mutate(
    cv_percent    = if_else(mean_conc_norm > 0, 100 * sd_conc_norm / mean_conc_norm, NA_real_),
    Ventilation_f = factor(Ventilation_Rate, levels = vent_rates)
  )

# =====================================================================
# Normalised deposition rate plot (scaled)
# Has AMPAS interval as two facets
# =====================================================================

ampas_dep_rate <- ampas_variability %>%
  mutate(
    interval_h    = AMPAS_Interval / 60,
    # Convert deposition into deposition rate
    dep_rate_mean = mean_conc_norm / interval_h,
    dep_rate_sd   = sd_conc_norm   / interval_h,
    Interval_f    = factor(AMPAS_Interval, levels = c(30, 60))
  )

# Global maximum max(mean + sd) across all ventilations and AMPAS intervals
global_max <- max(ampas_dep_rate$dep_rate_mean + ampas_dep_rate$dep_rate_sd, na.rm = TRUE)

# Scale the deposition rates
ampas_dep_rate <- ampas_dep_rate %>%
  mutate(
    dep_rate_mean_scaled = dep_rate_mean / global_max,
    dep_rate_sd_scaled   = dep_rate_sd   / global_max
  )

p_deposition_ampas_scaled <- ggplot(ampas_dep_rate, aes_ach_location()) +
  geom_point(aes(y = dep_rate_mean_scaled), size = 3, position = pd) +
  geom_errorbar(
    aes(
      ymin = pmax(dep_rate_mean_scaled - dep_rate_sd_scaled, 0),
      ymax = dep_rate_mean_scaled + dep_rate_sd_scaled
    ),
    width = 0.3, linewidth = 0.6, position = pd
  ) +
  facet_wrap2(
    ~ Interval_f,
    labeller = as_labeller(function(x) paste("AMPAS interval", x, "min")),
    axes = "all",
    nrow = 1
  ) +
  labs(
    x = "Plate location",
    y = expression(frac(dot(D), max(bar(dot(D)) + sigma(dot(D))))),
    colour = "Ventilation rate (ACH)",
    shape  = "Ventilation rate (ACH)"
  ) +
  theme_ampas +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_ach_colour + scale_ach_shape + ach_legend_tweaks

# =====================================================================
# CV plots
# =====================================================================

p_cv <- ampas_variability %>%
  filter(is.finite(cv_percent)) %>%
  ggplot(aes(
    x = factor(Location), y = cv_percent,
    colour = Ventilation_f, shape = Ventilation_f, group = Ventilation_f
  )) +
  geom_point(size = 3, position = pd) +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(
    ~ AMPAS_Interval,
    labeller = as_labeller(function(x) paste("AMPAS sampling interval:", x, "min"))
  ) +
  labs(
    x = "Plate location",
    y = "Coefficient of Variation (%)",
    colour = "Ventilation rate (ACH)",
    shape  = "Ventilation rate (ACH)"
  ) +
  theme_ampas + scale_ach_colour + scale_ach_shape + ach_legend_tweaks

# CV by interval across the locations
cv_by_interval <- ampas_variability %>%
  group_by(Ventilation_Rate, AMPAS_Interval) %>%
  summarise(
    mean_cv = mean(cv_percent, na.rm = TRUE),
    sd_cv   = sd(cv_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Ventilation_f = factor(Ventilation_Rate, levels = vent_rates))

p_cv_interval <- ggplot(cv_by_interval,
                        aes(
                          x = factor(AMPAS_Interval), y = mean_cv,
                          colour = Ventilation_f, shape = Ventilation_f, group = Ventilation_f
                        )) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(
    aes(ymin = mean_cv - sd_cv, ymax = mean_cv + sd_cv),
    width = 0.3, linewidth = 0.6, position = pd
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "AMPAS sampling interval (min)",
    y = "Mean CV (%)",
    colour = "Ventilation rate (ACH)",
    shape  = "Ventilation rate (ACH)"
  ) +
  theme_ampas + scale_ach_colour + scale_ach_shape + ach_legend_tweaks

# =====================================================================
# 30 vs 60 min deposition ratio
# =====================================================================

interval_compare <- ampas_variability %>%
  select(Ventilation_Rate, Location, AMPAS_Interval, mean_conc_norm, sd_conc_norm) %>%
  filter(AMPAS_Interval %in% c(30, 60)) %>%
  mutate(AMPAS_Interval = as.character(AMPAS_Interval)) %>%
  pivot_wider(
    names_from  = AMPAS_Interval,
    values_from = c(mean_conc_norm, sd_conc_norm),
    names_sep   = "_"
  ) %>%
  filter(mean_conc_norm_30 > 0) %>%
  mutate(
    ratio_60_to_30 = mean_conc_norm_60 / mean_conc_norm_30,
    ratio_sd = ratio_with_sd(
      mean_conc_norm_60, sd_conc_norm_60,
      mean_conc_norm_30, sd_conc_norm_30
    ),
    Ventilation_f = factor(Ventilation_Rate, levels = vent_rates)
  )

p_30vs60min_ratio <- ggplot(interval_compare, aes_ach_location()) +
  geom_point(aes(y = ratio_60_to_30), size = 3, position = pd) +
  geom_errorbar(
    aes(
      ymin = pmax(ratio_60_to_30 - ratio_sd, 0),
      ymax = ratio_60_to_30 + ratio_sd
    ),
    width = 0.3, linewidth = 0.6, position = pd
  ) +
  geom_hline(yintercept = 2, linetype = "dashed", colour = "grey20") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = "Plate location",
    y = "Deposition ratio (60 min / 30 min)",
    colour = "Ventilation rate (ACH)",
    shape  = "Ventilation rate (ACH)"
  ) +
  theme_ampas + scale_ach_colour + scale_ach_shape + ach_legend_tweaks

# =====================================================================
# Bland–Altman (60 vs 30 min) using µg cm^-2 h^-1
# =====================================================================

# Plate geometry and extract volume
plate_radius_cm <- 5.5 / 2          # 55 mm plate -> 2.75 cm radius
plate_area_cm2  <- pi * plate_radius_cm^2
extract_vol_L   <- 0.006            # 6 mL extract

# Compute deposition rate per plate in µg cm^-2 h^-1 (raw, not impinger-normalised)
ampas_dep_cm2 <- ampas_55mm %>%
  filter(
    Ventilation_Rate %in% vent_rates,
    AMPAS_Interval %in% c(30, 60),
    is.finite(conc_ppb)
  ) %>%
  mutate(
    # conc_ppb ≈ µg/L
    mass_ug       = conc_ppb * extract_vol_L,          # total µg in extract
    dep_ug_cm2    = mass_ug / plate_area_cm2,          # µg per cm^2
    dep_ug_cm2_h  = dep_ug_cm2 / (AMPAS_Interval / 60),# µg cm^-2 h^-1
    Ventilation_f = factor(Ventilation_Rate, levels = vent_rates),
    Location_f    = factor(Location)
  )

# Mean deposition per location × ACH × interval (µg cm^-2 h^-1)
ba_dep <- ampas_dep_cm2 %>%
  group_by(Ventilation_f, Location_f, AMPAS_Interval) %>%
  summarise(
    dep_mean = mean(dep_ug_cm2_h, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from   = AMPAS_Interval,
    values_from  = dep_mean,
    names_prefix = "min_"
  ) %>%
  rename(dep_30 = min_30, dep_60 = min_60) %>%
  filter(is.finite(dep_30), is.finite(dep_60))

# Bland–Altman data
ba_data <- ba_dep %>%
  mutate(
    mean_ba = (dep_30 + dep_60) / 2,   # mean of 30 & 60
    diff_ba = dep_60 - dep_30          # difference (60 - 30)
  )

# BA stats per ventilation (mean diff + limits of agreement)
ba_stats <- ba_data %>%
  group_by(Ventilation_f) %>%
  summarise(
    mean_diff  = mean(diff_ba, na.rm = TRUE),
    sd_diff    = sd(diff_ba,   na.rm = TRUE),
    loa_lower  = mean_diff - 1.96 * sd_diff,
    loa_upper  = mean_diff + 1.96 * sd_diff,
    .groups    = "drop"
  )

# Long format for horizontal lines (mean & LoA)
ba_lines <- ba_stats %>%
  tidyr::pivot_longer(
    cols      = c(mean_diff, loa_lower, loa_upper),
    names_to  = "type",
    values_to = "y"
  ) %>%
  mutate(
    linetype = dplyr::case_when(
      type == "mean_diff" ~ "mean",
      TRUE                ~ "loa"
    )
  )

# Join stats back for convenience (for y-limits)
ba_plot_data <- ba_data %>%
  left_join(ba_stats, by = "Ventilation_f")

# Common symmetric y-limit across both facets (includes LoA)
y_max <- max(abs(c(
  ba_plot_data$diff_ba,
  ba_stats$loa_lower,
  ba_stats$loa_upper
)), na.rm = TRUE)

# Custom strip labels
vent_labels <- c("1.5" = "1.5 ACH", "10" = "10 ACH")

# Line legend for mean + LoA
line_guide <- scale_linetype_manual(
  name   = NULL,
  values = c("mean" = "dashed", "loa" = "dotted"),
  labels = c(
    "mean" = "Mean difference",
    "loa"  = "95% limits of agreement"
  ),
  guide  = guide_legend(nrow = 2, byrow = TRUE)
)

p_bl_alt_60 <- ggplot() +
  geom_point(
    data = ba_plot_data,
    aes(
      x = mean_ba, y = diff_ba,
      colour = Ventilation_f,
      shape  = Ventilation_f
    ),
    size = 3
  ) +
  geom_hline(
    data = ba_lines,
    aes(yintercept = y, linetype = linetype),
    colour = "black",
    linewidth = 0.7
  ) +
  facet_wrap(
    ~ Ventilation_f,
    labeller = as_labeller(vent_labels),
    scales = "free_x"
  ) +
  scale_y_continuous(limits = c(-y_max, y_max)) +
  scale_ach_colour +
  scale_ach_shape +
  guides(
    colour = "none",
    shape  = "none"
  ) +
  line_guide +
  labs(
    x = expression(
      atop("Mean deposition rate*",
           paste("(", mu, "g cm"^{-2}, " h"^{-1}, ")"))
    ),
    y = expression(
      atop("Difference in deposition rate*",
           paste("(", mu, "g cm"^{-2}, " h"^{-1}, ")"))
    ),
    caption = "* 60 min vs 30 min sample"
  ) +
  theme_ampas +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.spacing.y  = unit(0.1, "cm"),
    legend.spacing.x  = unit(0.1, "cm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.key.height = unit(0.4, "lines"),
    legend.key.width  = unit(0.8, "lines"),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    strip.text = element_text(size = 12),
    plot.caption = element_text(
      hjust = 1,
      vjust = 0,
      size = 10,
      margin = margin(t = -14)
    )
  )

# =====================================================================
# 90 mm vs 55 mm plates
# =====================================================================

pac_subset <- plate_area_comparison %>%
  mutate(
    Type = case_when(
      grepl("55", Type, ignore.case = TRUE) ~ "55mm",
      grepl("90", Type, ignore.case = TRUE) ~ "90mm",
      TRUE ~ Type
    ),
    # Scale the ppb for the 55 mm that was in 6 mL to equivalent in 15 mL
    conc_ppb_scaled = if_else(Type == "55mm", conc_ppb * (6 / 15), conc_ppb)
  ) %>%
  filter(
    Ventilation_Rate %in% vent_rates,
    Type %in% c("55mm", "90mm"),
    is.finite(conc_ppb_scaled)
  )

pac_summary <- pac_subset %>%
  group_by(Ventilation_Rate, AMPAS_Interval, Location, Type) %>%
  summarise(
    n         = n(),
    mean_conc = mean(conc_ppb_scaled),
    sd_conc   = sd(conc_ppb_scaled),
    .groups = "drop"
  )

pac_summary_wide <- pac_summary %>%
  pivot_wider(
    names_from  = Type,
    values_from = c(mean_conc, sd_conc),
    names_sep   = "_"
  ) %>%
  filter(mean_conc_55mm > 0, mean_conc_90mm > 0) %>%
  mutate(
    ratio_90_over_55 = mean_conc_90mm / mean_conc_55mm,
    ratio_sd = ratio_with_sd(
      mean_conc_90mm, sd_conc_90mm,
      mean_conc_55mm, sd_conc_55mm
    ),
    Ventilation_f = factor(Ventilation_Rate, levels = vent_rates)
  )

p_90vs55mm_ratio <- ggplot(pac_summary_wide, aes_ach_location()) +
  geom_point(aes(y = ratio_90_over_55), size = 3, position = pd) +
  geom_errorbar(
    aes(
      ymin = pmax(ratio_90_over_55 - ratio_sd, 0),
      ymax = ratio_90_over_55 + ratio_sd
    ),
    width = 0.3, linewidth = 0.6, position = pd
  ) +
  facet_wrap(
    ~ AMPAS_Interval,
    labeller = as_labeller(c(
      "30" = "Static plate exposure time: 4 h",
      "60" = "Static plate exposure time: 5.5 h"
    ))
  ) +
  geom_hline(yintercept = area_ratio, linetype = "dashed", colour = "grey20") +
  scale_y_continuous(limits = c(0, 6), breaks = 0:6) +
  labs(
    x = "Plate location",
    y = "Deposition ratio (90 mm / 55 mm plates)",
    colour = "Ventilation rate (ACH)",
    shape  = "Ventilation rate (ACH)"
  ) +
  theme_ampas + scale_ach_colour + scale_ach_shape + ach_legend_tweaks

# =====================================================================
# Save images + print
# =====================================================================

save_plot_wide(p_deposition_ampas_scaled, "ampas_deposition_scaled")
save_plot_wide(p_cv,                     "cv_all")
save_plot_narrow(p_cv_interval,          "cv_interval")
save_plot_narrow(p_30vs60min_ratio,      "deposition_ratio_30vs60min")
save_plot_wide(p_bl_alt_60,              "bland_altman_60_vs_30_ug_cm2_h")
save_plot_wide(p_90vs55mm_ratio,         "deposition_ratio_90_vs_55mm")

print(p_cv)
print(p_cv_interval)
print(p_30vs60min_ratio)
print(p_bl_alt_60)
print(p_90vs55mm_ratio)
print(p_deposition_ampas_scaled)
