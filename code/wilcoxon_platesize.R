# =====================================================================
# 10/12/2025, C. A. Higham
# Nonparametric comparison of 90 mm vs 55 mm deposition 
# (Wilcoxon signed-rank test on location-averaged log-ratios)
# ====================================================================

# ----------------------------
# Packages
# ----------------------------
library(readr)
library(dplyr)
library(tidyr)

# ----------------------------
# Paths and constants
# ----------------------------

# Base project directory
base_dir <- "Your_Directory"

# ICP-MS data directory
data_dir <- file.path(base_dir, "results_csv", "deposition_data")

vent_rates <- c(1.5, 10)

# ----------------------------
# Load raw data
# ----------------------------
impinger_conc          <- read_csv(file.path(data_dir, "impinger_conc.csv"),          show_col_types = FALSE)
plate_area_comparison  <- read_csv(file.path(data_dir, "plate_area_comparison.csv"),  show_col_types = FALSE)

# ============================================================
# Impinger normalisation components
# ============================================================

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

impinger_max_by_vr <- impinger_summary_all %>%
  group_by(Ventilation_Rate) %>%
  summarise(
    max_control_ppb_per_h = max(control_ppb_per_h, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# Plate-area dataset (55 mm vs 90 mm)
# - Standardise Type labels (55mm / 90mm)
# - Rescale 55 mm plates to be comparable
# - Compute plate area (m^2)
# ============================================================

pac_subset <- plate_area_comparison %>%
  mutate(
    # Clean up plate type
    Type = case_when(
      grepl("55", Type, ignore.case = TRUE) ~ "55mm",
      grepl("90", Type, ignore.case = TRUE) ~ "90mm",
      TRUE ~ Type
    ),
    
    # Scale 55 mm plate concentrations to be comparable with 90 mm
    conc_ppb_scaled = if_else(Type == "55mm", conc_ppb * (6 / 15), conc_ppb),
    
    # Plate area in m^2 from diameter in mm
    diameter_mm   = as.numeric(gsub("mm", "", Type)),
    diameter_m    = diameter_mm / 1000,
    plate_area_m2 = pi * (diameter_m / 2)^2
  ) %>%
  filter(
    Ventilation_Rate %in% vent_rates,
    Type %in% c("55mm", "90mm"),
    is.finite(conc_ppb_scaled)
  )

# ============================================================
# Join impinger values and normalise
# - Normalise by impinger (ppb/h)
# - Scale by max control per ventilation rate
# - Convert to deposition per m^2 (up to scaling)
# ============================================================

pac_norm <- pac_subset %>%
  left_join(
    impinger_summary_all,
    by = c("Experiment_Index", "Ventilation_Rate")
  ) %>%
  left_join(
    impinger_max_by_vr,
    by = "Ventilation_Rate"
  ) %>%
  filter(
    !is.na(control_ppb_per_h),
    control_ppb_per_h > 0,
    !is.na(max_control_ppb_per_h)
  ) %>%
  mutate(
    conc_norm_raw  = conc_ppb_scaled / control_ppb_per_h,
    conc_norm_imp  = conc_norm_raw * max_control_ppb_per_h,
    conc_norm_area = conc_norm_imp / plate_area_m2   # deposition per m^2 (up to scaling)
  )

# ============================================================
# Build modelling dataset with factors
# ============================================================

mod_dat <- pac_norm %>%
  mutate(
    Type_f        = factor(Type, levels = c("55mm", "90mm")),
    Ventilation_f = factor(Ventilation_Rate, levels = vent_rates),
    Interval_f    = factor(AMPAS_Interval, levels = c(30, 60)),
    Experiment_f  = factor(Experiment_Index),
    Location_f    = factor(Location)
  )

# ============================================================
# Location-level averages and Wilcoxon signed-rank test
# - Average log deposition within Experiment × Location × Type
# - Compute log-ratios (90mm vs 55mm)
# - Test median log-ratio = 0 (i.e. 90mm/55mm ratio = 1)
# ============================================================

# Average log deposition within Experiment × Location × Type
loc_avg <- mod_dat %>%
  group_by(Experiment_f, Location_f, Type_f) %>%
  summarise(
    mean_log_dep = mean(log(conc_norm_area)),
    .groups = "drop"
  )

# Put 55 mm and 90 mm side by side, compute log-ratio and ratio
loc_wide <- loc_avg %>%
  pivot_wider(
    names_from  = Type_f,
    values_from = mean_log_dep
  ) %>%
  filter(!is.na(`55mm`), !is.na(`90mm`)) %>%
  mutate(
    log_ratio = `90mm` - `55mm`,  # log(90) - log(55) = log(90/55)
    ratio     = exp(log_ratio)    # 90mm / 55mm on original scale
  )

# Wilcoxon signed-rank test: is median log_ratio = 0 ?
wilcox_res <- wilcox.test(
  loc_wide$log_ratio,
  mu         = 0,
  conf.int   = TRUE,
  alternative = "two.sided"
)

cat("\n===== Wilcoxon signed-rank test on log-ratios (90mm vs 55mm) =====\n")
print(wilcox_res)

# 7.4 Convert estimate and CI to ratio scale (90mm / 55mm)
ratio_hat <- exp(wilcox_res$estimate)
ratio_CI  <- exp(wilcox_res$conf.int)

cat("\nEstimated median ratio (90mm / 55mm):\n")
print(ratio_hat)

cat("95% CI for ratio (90mm / 55mm):\n")
print(ratio_CI)
