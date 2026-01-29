# =====================================================================
# 10/12/2025, C. A. Higham, OPC analysis script for AMPAS paper
# =====================================================================
# Uses the OPC data taken from the OPC-N3 (Alphasense)
# Reads data at steady state nebuliser conc and background conc
# Calculates conc due to nebuliser, for both 1.5 ACH and 10 ACH
# Plots for both ventilations and a combined plot on log scale for y
# =====================================================================

# =====================================================================
# Packages
# =====================================================================
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)

# =====================================================================
# Plot theme and labels
# =====================================================================
theme_ampas <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

y_lab_conc <- expression("Concentration (m"^{-3}*")")

# =====================================================================
# Paths and output directory
# =====================================================================

# Base project directory
base_dir <- "Your_Directory"

# ICP-MS data directory
data_dir <- file.path(your_dir, "results_csv", "deposition_data")

# Save plot data directory
save_dir <- file.path(your_dir, "output_plots")

# Directory for the csv files for 10 ACH and 1.5 ACH
base_root_10ach  <- file.path(your_dir, "results_csv", "opc_data", "10ACH")
base_root_1_5ach <- file.path(your_dir, "results_csv", "opc_data", "1.5ACH")

exp_ids <- c("exp-1", "exp-2", "exp-3")

# =====================================================================
# Helper functions
# =====================================================================

# Convert counts per second (cps) to concentration per m³
n3_flow_rate <- (0.28 / 1000) / 60   # m³/s

to_per_m3 <- function(cps) {
  cps / n3_flow_rate
}

# Load an OPC CSV, reshape to long format, and compute the conc
load_opc_file <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)
  
  size_cols <- grep("um$", names(raw), value = TRUE)
  
  raw %>%
    pivot_longer(
      cols      = all_of(size_cols),
      names_to  = "Size_um",
      values_to = "cps"
    ) %>%
    mutate(
      Size = as.numeric(sub("um$", "", Size_um)),
      Conc = to_per_m3(cps)
    ) %>%
    select(Size, Conc)
}

# Summarise steady-state (SS) and background (BG) for one experiment
summarise_one_exp <- function(base_root, exp, ach) {
  ss <- load_opc_file(file.path(base_root, exp, "nebuliser-ss.csv")) %>%
    mutate(Type = "SS")
  
  bg <- load_opc_file(file.path(base_root, exp, "background-ss.csv")) %>%
    mutate(Type = "BG")
  
  bind_rows(ss, bg) %>%
    group_by(Size, Type) %>%
    summarise(
      mean_conc = mean(Conc, na.rm = TRUE),
      sd_conc   = sd(Conc, na.rm = TRUE),
      n         = n(),
      .groups   = "drop"
    ) %>%
    pivot_wider(
      names_from  = Type,
      values_from = c(mean_conc, sd_conc, n)
    ) %>%
    mutate(
      Experiment       = exp,
      Ventilation_Rate = ach,
      # SS − BG (truncated at 0)
      Neb_mean_exp     = pmax(mean_conc_SS - mean_conc_BG, 0),
      # Propagated SD for SS − BG
      Neb_sd_exp       = sqrt((sd_conc_SS^2 / n_SS) +
                                (sd_conc_BG^2 / n_BG))
    )
}

# Summarise all three experiments for one ventilation rate
summarise_all_exps <- function(root, ach) {
  map_dfr(exp_ids, ~ summarise_one_exp(root, .x, ach))
}

# =====================================================================
# Build vent_summary_all
# This is the mean and SD (across the 3 experiments) for #/m3 of
# particles at each vent rate and at each particle size
# =====================================================================

exp_10  <- summarise_all_exps(base_root_10ach, 10)
exp_1_5 <- summarise_all_exps(base_root_1_5ach, 1.5)

vent_summary_all <- bind_rows(exp_10, exp_1_5) %>%
  group_by(Ventilation_Rate, Size) %>%
  summarise(
    Neb_mean_vent = mean(Neb_mean_exp),
    Neb_sd_vent   = sqrt(sum(Neb_sd_exp^2)) / n(),  # propagated SD
    .groups       = "drop"
  )

# =====================================================================
# Plot save function
# =====================================================================

save_plot_wide <- function(plot, name) {
  ggsave(
    filename = file.path(save_dir, paste0(name, ".png")),
    plot     = plot,
    width    = 18.5,
    height   = 7,
    units    = "cm",
    dpi      = 400,
    bg       = "transparent"
  )
}

# =====================================================================
# Plotting data frame
# =====================================================================

opc_plot_df <- vent_summary_all %>%
  mutate(
    Size_num      = Size,
    Ventilation_f = factor(Ventilation_Rate, levels = c(1.5, 10))
  )

# The minor breaks for the x scale
minor_breaks_log <- c(seq(0.2, 0.9, 0.1), 2:9, seq(20, 90, 10))

# =====================================================================
# Legend centring helper (keeps your style consistent across plots)
# =====================================================================

theme_legend_centered <- theme(
  legend.position      = c(0.95, 0.95),
  legend.justification = c(1, 1),
  legend.background    = element_rect(fill = alpha("white", 0.7)),
  
  # Center the legend title and the legend entries
  legend.title         = element_text(size = 10, hjust = 0.5),
  legend.text          = element_text(size = 9,  hjust = 0.5),
  
  # Optional: helps the legend contents sit centered as a block
  legend.box.just      = "center"
)

# =====================================================================
# Function to create the individual ventilation rate plots
# =====================================================================

make_plot_separate_ACH <- function(df, ach_value) {
  df_ach <- df %>% filter(Ventilation_Rate == ach_value)
  
  ggplot(
    df_ach,
    aes(
      x      = Size_num,
      y      = Neb_mean_vent,
      colour = Ventilation_f,
      shape  = Ventilation_f
    )
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 3) +
    geom_errorbar( # SD for the plots
      aes(
        ymin = pmax(Neb_mean_vent - Neb_sd_vent, 0),
        ymax = Neb_mean_vent + Neb_sd_vent
      ),
      width     = 0.05,
      linewidth = 0.6
    ) +
    scale_x_log10( # x is on a log scale
      breaks       = c(0.1, 1, 10, 50),
      labels       = c("0.1", "1", "10", "50"),
      minor_breaks = minor_breaks_log
    ) +
    annotation_logticks(sides = "b", linewidth = 0.3) +
    scale_colour_manual(values = c("1.5" = "#009aDE", "10" = "#ff1f5b")) +
    scale_shape_manual(values  = c("1.5" = 16, "10" = 15)) +
    labs(
      colour = "Ventilation rate\n(ACH)",
      shape  = "Ventilation rate\n(ACH)",
      x      = "Particle size (µm)",
      y      = y_lab_conc
    ) +
    theme_ampas +
    theme_legend_centered
}

# =====================================================================
# Create and save individual ventilation plots
# =====================================================================

p_opc_1_5ach <- make_plot_separate_ACH(opc_plot_df, 1.5)
p_opc_10ach  <- make_plot_separate_ACH(opc_plot_df, 10)

print(p_opc_1_5ach)
print(p_opc_10ach)

save_plot_wide(p_opc_1_5ach, "opc_size_spectrum_1_5ach")
save_plot_wide(p_opc_10ach,  "opc_size_spectrum_10ach")

# =====================================================================
# Combined plot: 1.5 ACH and 10 ACH, log–x and log–y
# Create and save the plots
# =====================================================================

# Desired y-axis ticks: 10^3 – 10^7
y_breaks <- 10^(3:7)

p_opc_both <- ggplot(
  opc_plot_df,
  aes(
    x      = Size_num,
    y      = Neb_mean_vent,
    colour = Ventilation_f,
    shape  = Ventilation_f
  )
) +
  geom_line(aes(group = Ventilation_f), linewidth = 0.7) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = pmax(Neb_mean_vent - Neb_sd_vent, 0),
      ymax = Neb_mean_vent + Neb_sd_vent
    ),
    width     = 0.05,
    linewidth = 0.6
  ) +
  scale_x_log10(
    breaks       = c(0.1, 1, 10, 50),
    labels       = c("0.1", "1", "10", "50"),
    minor_breaks = minor_breaks_log
  ) +
  scale_y_log10(
    breaks = y_breaks,
    labels = expression(10^3, 10^4, 10^5, 10^6, 10^7)
  ) +
  annotation_logticks(sides = "b", linewidth = 0.3) +
  scale_colour_manual(values = c("1.5" = "#009aDE", "10" = "#ff1f5b")) +
  scale_shape_manual(values  = c("1.5" = 16, "10" = 15)) +
  labs(
    colour = "Ventilation rate\n(ACH)",
    shape  = "Ventilation rate\n(ACH)",
    x      = "Particle size (µm)",
    y      = y_lab_conc
  ) +
  theme_ampas +
  theme_legend_centered

suppressWarnings(print(p_opc_both))
suppressWarnings(save_plot_wide(p_opc_both, "opc_size_spectrum_both_ach"))
