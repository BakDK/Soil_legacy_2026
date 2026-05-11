# -----------------------------------------------------------------------------
# 03_Legacy_Biomarker_Rankings.R - Analyze SHAP results from Random Forest models
# trained on legacy fertilizer data to identify and rank potential biomarkers.
# -----------------------------------------------------------------------------
#
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




# ===================================================================
# --- 1. LOAD REQUIRED LIBRARIES ---
# ===================================================================
library(ggplot2)
library(dplyr)
library(readr)
library(viridis)
library(scales) 

# ===================================================================
# --- 2. CORE ANALYSIS FUNCTION ---
# ===================================================================
#' Analyzes SHAP results to find and plot top-ranked positive biomarkers.
#'
analyze_biomarkers <- function(analysis_prefix,
                               base_dir,
                               output_dir,
                               fold_change_threshold = 1.5,
                               n_top_features = 15,
                               n_runs_total = 20) {

  message(paste("\n--- Analyzing:", analysis_prefix, "---"))

  
  targeted_csv_path <- file.path(base_dir, paste0("rf_targeted_shap_summary_", analysis_prefix, ".csv"))
  if (!file.exists(targeted_csv_path)) {
    warning(paste("Targeted summary file not found, skipping:", basename(targeted_csv_path)))
    return(NULL)
  }
  targeted_summary <- readr::read_csv(targeted_csv_path, show_col_types = FALSE)
  if (nrow(targeted_summary) == 0) {
    message(paste("Targeted summary for", analysis_prefix, "is empty."))
    return(NULL)
  }
  class_order <- sort(unique(targeted_summary$Class))
  targeted_summary <- targeted_summary %>%
    mutate(Class = factor(Class, levels = class_order))
  
  qualified_biomarkers <- targeted_summary %>%
    filter(Avg_Fold_Change_TP_vs_nonTP >= fold_change_threshold)
  if (nrow(qualified_biomarkers) == 0) {
    message("No features met the fold change threshold.")
    return(NULL)
  }
  final_ranking_summary <- qualified_biomarkers %>%
    mutate(
      inter_run_robustness = Num_Runs_Found / n_runs_total,
      intra_run_robustness = Avg_Proportion_of_TPs_with_Feature,
      total_robustness_factor = inter_run_robustness * intra_run_robustness,
      robust_potency = Avg_Potency_SHAP * total_robustness_factor,
      robust_avg_fold_change = Avg_Fold_Change_TP_vs_nonTP * (Num_Runs_Found / n_runs_total)
    ) %>%
    arrange(Class, desc(robust_potency))

  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  ranked_csv_path <- file.path(output_dir, paste0("final_biomarker_ranking_", analysis_prefix, ".csv"))
  write_csv(final_ranking_summary, ranked_csv_path)
  message(paste("Saved final biomarker ranking to:", basename(ranked_csv_path)))

  # --- Step 5: Loop through each class to create individual plots ---
  for (current_class in class_order) { 
      class_data <- final_ranking_summary %>%
      filter(Class == current_class) %>%
      slice_max(order_by = robust_potency, n = n_top_features)

    if (nrow(class_data) == 0) {
      message(paste("  - No top biomarkers to plot for class:", current_class))
      next
    }
    
   
    plot_max_fc <- max(class_data$robust_avg_fold_change, na.rm = TRUE)
    if (plot_max_fc > 1000) {
      scale_limits <- c(NA, 1000)
      scale_breaks <- c(1.5, 10, 100, 1000)
      scale_labels <- c("1.5x", "10x", "100x", "≥1000x")
    } else {
      scale_limits <- c(NA, NA)
      scale_breaks <- pretty_breaks(n = 4)
      scale_labels <- waiver()
    }


    individual_plot <- ggplot(class_data,
                           aes(x = reorder(Feature, robust_potency),
                               y = robust_potency,
                               fill = robust_avg_fold_change)) +
      geom_col() +
      coord_flip() +
      scale_fill_viridis_c(
        name = "Robust Fold Change",
        trans = "log10",
        limits = scale_limits,
        oob = scales::squish,
        breaks = scale_breaks,
        labels = scale_labels
      ) +
      labs(
        title = paste("Top Positive Biomarkers for:", analysis_prefix),
        subtitle = paste("Class:", current_class, "| Pre-filtered by Fold Change (>=", fold_change_threshold, ")"),
        x = NULL,
        y = "Robust Potency (True Average SHAP Strength)"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right")
    
    # --- SAVING LOGIC: Save both PNG and RDS files ---

    # 1. Save the visual PNG file
    png_path <- file.path(output_dir, paste0("final_biomarker_plot_", analysis_prefix, "_", current_class, ".png"))
    ggsave(png_path, individual_plot, width = 8, height = 7, dpi = 300, bg = "white")
    message(paste("  - Saved PNG plot for class", current_class, "to:", basename(png_path)))

    # 2. Save the R plot object as an RDS file
    rds_path <- file.path(output_dir, paste0("final_biomarker_plot_", analysis_prefix, "_", current_class, ".rds"))
    saveRDS(individual_plot, file = rds_path)
    message(paste("  - Saved RDS plot object for class", current_class, "to:", basename(rds_path)))
    # ===========================
  }
}


# ===================================================================
# --- SCRIPT EXECUTION ORCHESTRATOR 
# ===================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript 03_create_rankings.R <input_dir> <output_dir>", call. = FALSE)
}

input_dir <- args[1]
output_dir <- args[2]

message("--- R Script: Creating Final Rankings ---")
message("Reading Python results from: ", input_dir)
message("Saving final plots and CSVs to: ", output_dir)

targeted_files <- list.files(
  path = input_dir,
  pattern = "^rf_targeted_shap_summary_.*\\.csv$",
  full.names = FALSE
)

if (length(targeted_files) == 0) {
  stop("No 'rf_targeted_shap_summary' files found in the input directory: ", input_dir)
}

prefixes_to_analyze <- sub(
  pattern = "\\.csv$", 
  replacement = "", 
  sub(pattern = "^rf_targeted_shap_summary_", replacement = "", targeted_files)
)

message("\nFound the following prefixes to analyze: ", paste(prefixes_to_analyze, collapse=", "))

for (prefix in prefixes_to_analyze) {
  tryCatch({
    analyze_biomarkers(
      analysis_prefix = prefix,
      base_dir = input_dir,
      output_dir = output_dir
    )
  }, error = function(e) {
    message(paste("An error occurred during analysis for", prefix, ":", e$message))
  })
}

message("\n\nAll R analyses complete!")