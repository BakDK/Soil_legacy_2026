# -----------------------------------------------------------------------------
# 01GeneratingDataSetsForLegacyAnalysis.R
# -----------------------------------------------------------------------------
#
# Author:      Thomas Martini Jørgensen <tmjq@dtu.dk>
#
# Copyright (c) 2025, Thomas Martini Jørgensen / DTU Compute, Technical University of Denmark
#
# Copyright (c) 2026, Thomas Martini Jørgensen / DTU Compute, Technical University of Denmark
#
# This script is part of the code accompanying the paper:
# "Identification of robust bacterial biomarkers in the wheat root microbiome for predicting soil legacy effects and crop development "
# by Frederik Bak, Thomas Martini Jørgensen, Inês Nunes, Veronika Hansen, and  Mette Haubjerg Nicolaisen
#
# GitHub Repository: https://github.com/BakDK/Soil_legacy_2026 
#
# This script generates the data sets from a phyloseq object that will be used  to train random forest models performs to predict fertilizer treatment. 
# It creates a base training dataset from early time points ("T0", "T1", "3DBS"), splits it into training and testing sets. The split into training and test can be done using ranodom split or based on blocks.
# In addition dividual test sets are created for later time points to evaluate if the # fertilizer "legacy" is still predictable as the plant grows.
# All data is processed # at the genus level.
#
# This code is licensed under the MIT License. See the LICENSE file in the
# root directory of this repository for the full license text.
#




library(tidyverse)
library(phyloseq)
library(stringr)

#' Create and Split Datasets to Predict Fertilizer Legacy Over Time
#'
#' This function prepares datasets from a phyloseq object to predict fertilizer treatment
#' based on microbial composition. It creates a base training dataset from early
#' time points ("T0", "T1", "3DBS"), splits it into training and testing sets,
#' and then creates individual test sets for later time points to evaluate if the
#' fertilizer "legacy" is still predictable as the plant grows. All data is processed
#' at the genus level.
#'
#' @param rds_filepath Path to the phyloseq object .rds file.
#' @param split_method Method for splitting the base data. Either "random" for an 80/20 split
#'   or "block" for splitting by specific block numbers.
#' @param train_ratio The proportion of data for training when split_method is "random".
#' @param train_blocks A numeric vector of block numbers for the training set.
#' @param test_blocks A numeric vector of block numbers for the test set.
#' @param output_dir The directory where the resulting CSV files will be saved.
#'
#' @return A list containing all generated data frames: train_data, main_test_data,
#'   and a list of individual_test_sets.
#'
create_fertilizer_legacy_datasets <- function(
    rds_filepath, 
    split_method = "random", # "random" or "block"
    train_ratio = 0.8,
    train_blocks = NULL, 
    test_blocks = NULL,
    output_dir = "data/training_and_test_data/"
) {
    
    # --- 0. Input Validation and Setup ---
    message("--- Starting Data Preparation for Fertilizer Legacy Prediction ---")
    if (!split_method %in% c("random", "block")) {
        stop("Error: split_method must be either 'random' or 'block'.")
    }
    if (split_method == "block" && (is.null(train_blocks) || is.null(test_blocks))) {
        stop("Error: For 'block' split method, you must provide train_blocks and test_blocks.")
    }
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # --- 1. Import Phyloseq Object and Pre-process Metadata ---
    message(paste("Reading phyloseq object from:", rds_filepath))
    if (!file.exists(rds_filepath)) stop(paste("Error: RDS file not found at", rds_filepath))
    physeq_obj <- readRDS(rds_filepath)
    if (!inherits(physeq_obj, "phyloseq")) stop("Object is not a phyloseq object.")

    meta_df <- data.frame(sample_data(physeq_obj))
    sample_id_source <- if ("Sample" %in% colnames(meta_df)) meta_df$Sample else rownames(meta_df)
    matches <- str_match(sample_id_source, "^[^.]+\\.[A-Za-z](\\d+)")
    meta_df$Block <- as.numeric(matches[, 2])
    if(all(is.na(meta_df$Block))) {
        warning("Block extraction failed. Please double-check the regex and your sample names.")
    }
    sample_data(physeq_obj) <- sample_data(meta_df)

    # --- 2. Filter for Base Training/Testing Data from Early Time Points ---
    base_times <- c("T0", "T1", "3DBS")
    message(paste("Filtering samples to create base dataset from Time points:", paste(base_times, collapse = ", ")))
    
    # Get the names of the samples that match our criteria first
    samples_to_keep <- sample_names(physeq_obj)[get_variable(physeq_obj, "Time") %in% base_times]
    # Then, create the new phyloseq object by keeping only those samples
    physeq_base <- prune_samples(samples_to_keep, physeq_obj)
    
    if (nsamples(physeq_base) == 0) {
        stop("No samples found for the specified base time points. Halting execution.")
    }
    message(paste("Base dataset contains", nsamples(physeq_base), "samples."))

    # --- 3. Process to Genus Level with NA Handling ---
    tax_df_base <- as.data.frame(tax_table(physeq_base))
    na_genus_indices <- which(is.na(tax_df_base$Genus))
    
    if (length(na_genus_indices) > 0) {
        message(paste("Found", length(na_genus_indices), "ASVs with NA Genus. Creating unique ASV-level labels..."))
        original_asv_names <- rownames(tax_df_base)[na_genus_indices]
        new_na_names <- paste0("Genus_NA_Rep_", original_asv_names)
        tax_df_base$Genus[na_genus_indices] <- new_na_names
        tax_table(physeq_base) <- tax_table(as.matrix(tax_df_base))
    }
    
    physeq_glommed <- tax_glom(physeq_base, taxrank = "Genus", NArm = FALSE)
    final_taxa_names <- make.unique(as.character(data.frame(tax_table(physeq_glommed))$Genus))
    taxa_names(physeq_glommed) <- final_taxa_names

    # --- 4. Create Full Machine Learning Data Frame ---
    otu_tab <- otu_table(physeq_glommed)
    features <- if (taxa_are_rows(otu_tab)) data.frame(t(otu_tab)) else data.frame(otu_tab)
    
    meta_filtered <- data.frame(sample_data(physeq_glommed))
    response_var <- factor(meta_filtered$Fertilizer)
    block_number <- as.factor(meta_filtered$Block)
    ml_data_full <- data.frame(features, response_var, block_number)
    
    # --- 5. Split Base Data into Training and Main Test Set ---
    message(paste("Splitting base data using the '", split_method, "' method.", sep=""))
    
    if (split_method == "random") {
        set.seed(42) 
        train_indices <- sample(1:nrow(ml_data_full), size = floor(train_ratio * nrow(ml_data_full)))
        test_indices <- setdiff(1:nrow(ml_data_full), train_indices)
        
        
    } else if (split_method == "block") {
        train_indices <- which(ml_data_full$block_number %in% train_blocks)
        test_indices <- which(ml_data_full$block_number %in% test_blocks)
        if (length(intersect(train_blocks, test_blocks)) > 0) {
            warning("Training and testing blocks overlap. This may lead to data leakage.")
        }
       
    }
    
    if (length(train_indices) == 0 || length(test_indices) == 0) {
        stop("Data splitting resulted in an empty train or test set. Check your parameters.")
    }
    
    train_data <- ml_data_full[train_indices, ]
    main_test_data <- ml_data_full[test_indices, ]
    train_data$block_number <- NULL
    main_test_data$block_number <- NULL
    
    # --- 6. Create Individual Test Sets for Future Time Points ---
    future_times <- c("11DAS", "14DAS", "17DAS", "31DAS", "158DAS", "187DAS", "215DAS", "227DAS", "256DAS", "318DAS")
    message("\n--- Creating individual test sets for future time points ---")
    
    individual_test_sets <- list()
    train_feature_names <- colnames(train_data)[colnames(train_data) != "response_var"]
    
    for (time_point in future_times) {
        message(paste("Processing:", time_point))
        
        samples_to_keep_future <- sample_names(physeq_obj)[get_variable(physeq_obj, "Time") == time_point]
        physeq_time <- prune_samples(samples_to_keep_future, physeq_obj)
    

        if (nsamples(physeq_time) == 0) {
            warning(paste("No samples found for", time_point, ". Skipping."))
            next
        }
        
        meta_time <- data.frame(sample_data(physeq_time))

        tax_df_time <- as.data.frame(tax_table(physeq_time))
        na_genus_indices_time <- which(is.na(tax_df_time$Genus))
        if (length(na_genus_indices_time) > 0) {
            original_asv_names_time <- rownames(tax_df_time)[na_genus_indices_time]
            new_na_names_time <- paste0("Genus_NA_Rep_", original_asv_names_time)
            tax_df_time$Genus[na_genus_indices_time] <- new_na_names_time
            tax_table(physeq_time) <- tax_table(as.matrix(tax_df_time))
        }
        physeq_glommed_time <- tax_glom(physeq_time, taxrank = "Genus", NArm = FALSE)
        final_taxa_names_time <- make.unique(as.character(data.frame(tax_table(physeq_glommed_time))$Genus))
        taxa_names(physeq_glommed_time) <- final_taxa_names_time
        
        otu_time <- otu_table(physeq_glommed_time)
        features_time <- if (taxa_are_rows(otu_time)) data.frame(t(otu_time)) else data.frame(otu_time)
        
        aligned_df <- as.data.frame(matrix(0, nrow = nrow(features_time), ncol = length(train_feature_names)))
        colnames(aligned_df) <- train_feature_names
        rownames(aligned_df) <- rownames(features_time)
        
        common_cols <- intersect(train_feature_names, colnames(features_time))
        aligned_df[, common_cols] <- features_time[, common_cols]
        
        aligned_df$response_var <- factor(meta_time$Fertilizer)
        
        individual_test_sets[[time_point]] <- aligned_df
    }

    # --- 7. Save All Datasets to Files ---
    message("\n--- Saving all datasets to CSV files ---")
    
    write.csv(train_data, file.path(output_dir, "train_data_fertilizer.csv"), row.names = TRUE)
    message(paste("Saved train_data_fertilizer.csv with", nrow(train_data), "samples."))

    main_test_filename <- "test_data_main_fertilizer.csv"
    write.csv(main_test_data, file.path(output_dir, main_test_filename), row.names = TRUE)
    message(paste("Saved", main_test_filename, "with", nrow(main_test_data), "samples."))
    
    for (time_point in names(individual_test_sets)) {
        df_to_save <- individual_test_sets[[time_point]]
        filename <- file.path(output_dir, paste0("test_set_fertilizer_", time_point, ".csv"))
        write.csv(df_to_save, filename, row.names = TRUE)
        message(paste("Saved", basename(filename), "with", nrow(df_to_save), "samples."))
    }
    
    message(paste("\nData preparation complete. All files saved in '", output_dir, "' directory.", sep=""))
    
    return(list(
        train_data = train_data, 
        main_test_data = main_test_data, 
        individual_test_sets = individual_test_sets
    ))
}


# ===================================================================
# --- SCRIPT EXECUTION: Define parameters and call the function ---
# ===================================================================


# ===================================================================
# --- SCRIPT: 01_GeneratingDataSetsForLegacyAnalysis.R ---
#
# Description: Reads a phyloseq object and generates training/test
#              datasets. It creates a base dataset from early time 
#              points and individual test sets for future time points.
#              Splitting can be done randomly or by specific blocks.
#
# Usage:
#   Rscript src/r/01_GeneratingDataSetsForLegacyAnalysis.R <path_to_rds> <split_method> <output_dir> <train_blocks> <test_blocks>
#
# Arguments:
#   path_to_rds:  The full path to the input .rds file (phyloseq object).
#   split_method: The method for splitting. Must be 'random' or 'block'.
#   output_dir:   The directory where output CSV files will be saved.
#   train_blocks: Comma-separated string of block IDs for training (e.g., "1,3"). 
#                 Use "none" if split_method is 'random'.
#   test_blocks:  Comma-separated string of block IDs for testing (e.g., "2,4"). 
#                 Use "none" if split_method is 'random'.
#
# Examples:
#   # For Block-based split:
#   Rscript src/r/01_GeneratingDataSetsForLegacyAnalysis.R \
#     data/input_files/phyloseq_genus.rds block outputs/legacy/block_split/ "1,3" "2,4"
#
#   # For Random split:
#   Rscript src/r/01_GeneratingDataSetsForLegacyAnalysis.R \
#     data/input_files/phyloseq_genus.rds random outputs/legacy/random_split/ "none" "none"
# ===================================================================

# --- 1. Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)

# Validate the number of arguments (Now expecting 5)
if (length(args) != 5) {
  stop("Usage: Rscript GeneratingDataSetsForLegacyAnalysis.r <path_to_rds> <split_method> <output_dir> <train_blocks> <test_blocks>", call. = FALSE)
}

# Assign arguments
phyloseq_file_path <- args[1]
split_method       <- args[2]
output_dir         <- args[3]
train_blocks_raw   <- args[4] # e.g., "1,3"
test_blocks_raw    <- args[5] # e.g., "2,4"

# --- 2. Input Validation ---
if (!file.exists(phyloseq_file_path)) {
  stop("Input file not found: ", phyloseq_file_path, call. = FALSE)
}

if (!(split_method %in% c("random", "block"))) {
  stop("Invalid split_method. Must be 'random' or 'block'.", call. = FALSE)
}

# Helper to convert "1,3" into c(1, 3)
parse_blocks <- function(block_string) {
  if (block_string == "none" || block_string == "") return(NULL)
  return(as.numeric(unlist(strsplit(block_string, ","))))
}

blocks_for_training <- parse_blocks(train_blocks_raw)
blocks_for_testing  <- parse_blocks(test_blocks_raw)

# --- 3. Prepare Output Directory ---
if (!dir.exists(output_dir)) {
  message("Output directory not found. Creating it now: ", output_dir)
  dir.create(output_dir, recursive = TRUE)
}

# --- 4. Conditional Script Execution ---
message(">>> Reading input file: ", phyloseq_file_path)
message(">>> Running with split method: ", toupper(split_method))

if (split_method == "random") {
  output_list <- create_fertilizer_legacy_datasets(
    rds_filepath = phyloseq_file_path,
    split_method = "random",
    output_dir = output_dir
  )
} else if (split_method == "block") {
  message(">>> Training Blocks: ", paste(blocks_for_training, collapse=", "))
  message(">>> Testing Blocks:  ", paste(blocks_for_testing, collapse=", "))
  
  output_list <- create_fertilizer_legacy_datasets(
    rds_filepath = phyloseq_file_path,
    split_method = "block",
    train_blocks = blocks_for_training,
    test_blocks  = blocks_for_testing,
    output_dir   = output_dir
  )
}

message("\n--- Dataset generation complete. ---")
#message("Output files are saved in: ", output_dir)