# -----------------------------------------------------------------------------
# 03_create_growth_stage_rankings.R 
# Creates final biomarker rankings and INDIVIDUAL plots from Python results.
# -----------------------------------------------------------------------------
# Author:      Thomas Martini Jørgensen <tmjq@dtu.dk>
#
# Copyright (c) 2026, Thomas Martini Jørgensen / DTU Compute, Technical University of Denmark
#
# This script is part of the code accompanying the paper:
# "Identification of robust bacterial biomarkers in the wheat root microbiome for predicting soil legacy effects and crop development "
# by Frederik Bak, Thomas Martini Jørgensen, Inês Nunes, Veronika Hansen, and  Mette Haubjerg Nicolaisen
#
# GitHub Repository: https://github.com/BakDK/Soil_legacy_2026 
#


# This code is licensed under the MIT License. See the LICENSE file in the
# root directory of this repository for the full license text.
#
# -----------------------------------------------------------------------------


# ===================================================================
# --- 1. LOAD REQUIRED LIBRARIES ---
# ===================================================================
# Set a non-interactive backend for plotting when run from a script
options(bitmapType = 'cairo') 
library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(viridis)

# ===================================================================
# --- SCRIPT EXECUTION ORCHESTRATOR ---
# ===================================================================

# --- 1. Argument Parsing ---
# This script is designed to be called with two directory arguments,
# as the original data is no longer needed for this reduced analysis.
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("FATAL: This script requires exactly 2 arguments: <python_results_dir> <final_output_dir>", call. = FALSE)
}

python_input_dir <- args[1]
output_dir <- args[2]

message("--- R Script: Creating Final Rankings & INDIVIDUAL Plots ---")
message("Reading Python results from: ", python_input_dir)
message("Saving final plots and CSVs to: ", output_dir)

# Ensure the final output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- 2. Load and Prepare Core Data ---
plot_title_prefix <- "RandomForest"
targeted_csv_path <- file.path(python_input_dir, "rf_targeted_shap_summary.csv")

if (!file.exists(targeted_csv_path)) {
  stop("FATAL: Targeted SHAP summary file from Python not found in '", python_input_dir, "'")
}

targeted_summary <- read_csv(targeted_csv_path, show_col_types = FALSE)

class_order <- c("early", "middle", "late")
targeted_summary <- targeted_summary %>% 
  mutate(Class = factor(Class, levels = class_order))

# --- 3. Perform Definitive "Filter-Then-Rank" Analysis ---
message("\n--- Performing definitive biomarker ranking ---")

abundance_percentile_threshold <- 0.0
fold_change_threshold <- 1.5
n_runs_total <- 20

qualified_biomarkers <- targeted_summary %>%
  filter(Avg_Fold_Change_TP_vs_nonTP >= fold_change_threshold) %>%
  group_by(Class) %>%
  filter(Avg_Mean_Abund_in_High_Abund_TPs >= quantile(Avg_Mean_Abund_in_High_Abund_TPs, probs = abundance_percentile_threshold, na.rm = TRUE)) %>%
  ungroup()

final_ranking_summary <- qualified_biomarkers %>%
  mutate(
    inter_run_robustness = Num_Runs_Found / n_runs_total,
    intra_run_robustness = Avg_Proportion_of_TPs_with_Feature,
    total_robustness_factor = inter_run_robustness * intra_run_robustness,
    robust_potency = Avg_Potency_SHAP * total_robustness_factor,
    robust_avg_fold_change = Avg_Fold_Change_TP_vs_nonTP * (Num_Runs_Found / n_runs_total)
  ) %>%
  arrange(Class, desc(robust_potency))

# --- 4. Save the Required Ranking CSV Files ---
write_csv(final_ranking_summary, file.path(output_dir, "rf_final_biomarker_ranking_FULL.csv"))

top_15_final_ranking <- final_ranking_summary %>% 
  group_by(Class) %>% 
  slice_max(order_by = robust_potency, n = 15) %>% 
  ungroup()
write_csv(top_15_final_ranking, file.path(output_dir, "rf_top_15_final_biomarkers.csv"))

message("Saved final biomarker ranking summaries (FULL and TOP 15) to CSV.")

# --- 5. Create and Save Individual Plots and R Objects ---
message("\n--- Generating and saving individual plots for each class ---")

# Initialize an empty list to store the ggplot objects
individual_plots_list <- list()

for (current_class in class_order) {
  
  message(paste("  - Generating plot for '", current_class, "' class", sep=""))
  
  # Filter the top 15 data for the current class
  class_data <- top_15_final_ranking %>% filter(Class == current_class)
  
  if (nrow(class_data) == 0) {
    message(paste("    - No data found for class '", current_class, "'. Skipping plot.", sep=""))
    next
  }
  
  # Create the individual plot (logic from your post-processing script)
  individual_plot <- ggplot(class_data, 
                           aes(x = reorder(Feature, robust_potency), 
                               y = robust_potency,
                               fill = robust_avg_fold_change)) +
    geom_col() +
    coord_flip() +
    scale_fill_viridis_c(
        name = "Robust Fold Change", 
        trans = "log10",
        limits = c(NA, 1000), 
        oob = scales::squish, 
        breaks = c(1.5, 10, 100, 1000),
        labels = c("1.5x", "10x", "100x", "≥1000x")
    ) +
    labs(
      title = paste("Top", current_class, "Biomarkers Ranked by Model Importance"),
      subtitle = paste("Data from", plot_title_prefix, "analysis"),
      x = NULL, 
      y = "Robust Potency (True Average SHAP Strength)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")

  # Save the individual PNG file
  output_filename_png <- paste0("rf_final_biomarker_ranking_plot_", current_class, ".png")
  output_filepath_png <- file.path(output_dir, output_filename_png)
  ggsave(output_filepath_png, individual_plot, width = 8, height = 6, dpi = 300, bg = "white")
  message(paste("    - Saved PNG plot to:", output_filepath_png))
  
  # Add the plot object to our list, naming it by its class
  individual_plots_list[[current_class]] <- individual_plot
}

# --- 6. Save the List of Plot Objects to a Single RDS File ---
rds_filepath <- file.path(output_dir, "rf_individual_ranking_plots.rds")
saveRDS(individual_plots_list, file = rds_filepath)
message(paste("\nSaved list of individual plot objects to:", rds_filepath))

message("\n\nAll R analyses for this fold are complete!")