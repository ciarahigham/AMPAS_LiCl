# **Spatiotemporal Characterisation of Aerosol Deposition under Controlled Ventilation: Performance Evaluation of the AMPAS Device**

Dr Ciara A. Higham<sup>1,*</sup>, Andrew Carnegie<sup>1</sup>, Dr Waseem Hiwar<sup>1</sup>, Dr Louise Fletcher<sup>1</sup>, Prof. Catherine J. Noakes<sup>1</sup>, Dr Marco-Felipe King<sup>1</sup>

<sup>1</sup> School of Civil Engineering, University of Leeds, Woodhouse Lane, Leeds, LS2 9JT, United Kingdom. 

<sup>2</sup> AirSentry Limited, Hermes House, Fire Fly Avenue, Swindon, SN2 2GA, United Kingdom.


<sup>*</sup> Corresponding email: ciara.higham@hotmail.co.uk

Supporting code and data for "Spatiotemporal Characterisation of Aerosol Deposition under Controlled Ventilation: Performance Evaluation of the AMPAS Device".

doi:

Analysis was performed in **R 4.5.1** using **RStudio 2025.05.1**.

---

## Repository Contents

## Repository Contents

### `experiment_details/`
- `.xlsx` file containing metadata for each experiment, including temperature and relative humidity.

### `results_csv/deposition_data/`
Raw ICP-MS lithium concentration datasets:
- `ampas_55mm.csv` – Li concentration extracted from 55 mm AMPAS plates.  
- `impinger_conc.csv` – Impinger Li concentrations and final volumes.  
- `initial_li_conc.csv` – Initial Li nebuliser solution concentrations.  
- `plate_area_comparison.csv` – Li concentrations used for the 90 mm vs 55 mm plate comparison.

### `results_csv/opc_data/`
OPC-N3 particle count data for each ventilation rate and experiment. Each file contains particle counts (each row count per second) across OPC size bins (columns).
- `background-ss.csv` – Background steady-state particle counts.  
- `nebuliser-ss.csv` – Nebuliser steady-state particle counts.

---

## Code Descriptions

### `code/opc_concentrations.R`
- Converts OPC-N3 particle counts to particle concentrations (m⁻³) for background and nebuliser steady-state conditions.
- Computes net nebuliser particle concentrations by subtracting background levels.
- Aggregates results across triplicate experiments for each ventilation rate.
- Produces size-resolved particle concentration spectra for **1.5 ACH** and **10 ACH**, including a combined comparison plot.
- Saves all figures as `.png` files in the `output_plots/` directory.

---

### `code/experimental_output_code.R`
Implements the AMPAS deposition analysis described in the Methods:

- Reads AMPAS plate (`ampas_55mm.csv`), impinger (`impinger_conc.csv`), and plate comparison (`plate_area_comparison.csv`) ICP-MS data.
- Normalises AMPAS concentrations using impinger concentrations to account for variability in airborne concentration between experiments.
- Aggregates triplicate samples at each grid location to obtain mean, standard deviation, and coefficient of variation (CV) of the normalised deposition for each ventilation rate and sampling interval.
- Converts normalised deposition to a deposition rate by dividing by the AMPAS sampling interval and scales by the global maximum of (mean + SD) to generate spatial deposition plots.
- Computes and plots:
  - CV at each location and mean ± SD CV across locations for each sampling interval and ventilation rate.
  - Deposition ratios between 60 min and 30 min AMPAS samples (60/30), with propagated uncertainty.
  - A Bland–Altman analysis comparing 60 min vs 30 min deposition rates, summarised per location and ventilation rate (mean difference and 95% limits of agreement shown as horizontal lines).
  - Deposition ratios between 90 mm and 55 mm plates, correcting for extraction volume and including the theoretical plate area ratio as a reference line.
- Saves all figures as `.png` files in the `output_plots/` directory.

---

### `code/lmer_model.R`
- Reconstructs the normalised deposition rate Ḋ (Ddot) from AMPAS plate and impinger data.
- Applies a natural log-transform to Ḋ and assesses normality and homoscedasticity using Shapiro–Wilk tests and Q–Q plots.
- Fits linear mixed-effects models with ventilation rate, sampling interval, and their interaction as fixed effects, and grid location and experiment index as random intercepts.
- Uses maximum-likelihood likelihood ratio tests to simplify the model, then refits the final reduced model with REML.
- Performs Type III Wald tests for fixed effects and generates diagnostic plots for residuals and random effects (saved to `output_plots/`).

---

### `code/wilcoxon_platesize.R`
- Normalises AMPAS deposition data using impinger concentrations and plate area.
- Computes location-averaged deposition for 90 mm and 55 mm AMPAS plates.
- Compares plate sizes using a Wilcoxon signed-rank test on log-transformed deposition ratios.
- Reports the estimated median deposition ratio (90 mm / 55 mm) with 95% confidence intervals.

---
