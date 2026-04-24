
# -----------------------------------------------------------------------------
# 01_generate_growth_stage_datasets.R - Prepares training and test sets for
# growth stage prediction based on block cross-validation.
# -----------------------------------------------------------------------------
#
# Author:      Thomas Martini Jørgensen <tmjq@dtu.dk>
#
# Copyright (c) 2026, Thomas Martini Jørgensen / DTU Compute, Technical University of Denmark
#
# This script is part of the code accompanying the paper:
# "Identification of robust bacterial biomarkers in the wheat root microbiome for predicting soil legacy and crop development "
# by Frederik Bak, Thomas Martini Jørgensen, Inês Nunes, Veronika Hansen, and  Mette Haubjerg Nicolaisen
#
# GitHub Repository: https://https://github.com/BakDK/Future_cropping
#
# This code is licensed under the MIT License. See the LICENSE file in the
# root directory of this repository for the full license text.
#
# -----------------------------------------------------------------------------

# --- >>> LOAD REQUIRED LIBRARIES <<< ---
library(tidyverse)
library(dplyr)
library(phyloseq)
library(stringr)

# ===================================================================
# --- CORE ANALYSIS FUNCTION ---
## ===================================================================
prepare_dataset_phyloseq_from_file <- function(
    rds_filepath, 
    fertilizer_types, 
    train_blocks, 
    test_blocks,
    output_dir  
  ) {
    
    # --- 0. Import the phyloseq object ---
    message(paste("Reading phyloseq object from:", rds_filepath))
    if (!file.exists(rds_filepath)) stop(paste("Error: RDS file not found at", rds_filepath))
    physeq_obj <- readRDS(rds_filepath)
    if (!inherits(physeq_obj, "phyloseq")) stop("Object is not a phyloseq object.")

    # --- 1. Add/Modify Metadata ---
    meta_df <- data.frame(sample_data(physeq_obj))
    meta_df <- meta_df %>%
      mutate(TimePeriod = case_when(
          Time %in% c("11DAS", "14DAS", "17DAS", "31DAS") ~ "early",
          Time %in% c("158DAS", "187DAS", "215DAS", "227DAS") ~ "middle",
          Time %in% c("256DAS", "318DAS") ~ "late",
          TRUE ~ NA_character_
        ))

    sample_id_source <- if ("Sample" %in% colnames(meta_df)) meta_df$Sample else rownames(meta_df)
    matches <- str_match(sample_id_source, "^[^.]+\\.[A-Za-z](\\d+)")
    meta_df$Block <- as.numeric(matches[, 2])

    if(all(is.na(meta_df$Block))) {
        warning("Block extraction failed. Please double-check the regex and your sample names.")
    }
    sample_data(physeq_obj) <- sample_data(meta_df)

    # --- 2. Filter Samples (The Robust dplyr Method) ---
    message(paste("Filtering samples for fertilizer type(s):", paste(fertilizer_types, collapse = ", ")))
    
    sample_ids_to_keep <- data.frame(sample_data(physeq_obj)) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      filter(
        !is.na(.data[["TimePeriod"]]),
        .data[["TimePeriod"]] %in% c("early", "middle", "late"),
        .data[["Fertilizer"]] %in% fertilizer_types # Use the argument, not a global variable
      ) %>%
      pull(Sample_ID)
      
    if (length(sample_ids_to_keep) == 0) {
      stop(paste("No samples remained after filtering for Fertilizers  =  ", paste(fertilizer_types, collapse = ", ")))
    }
    physeq_filtered <- prune_samples(sample_ids_to_keep, physeq_obj)
    message(paste("Samples remaining after filtering:", nsamples(physeq_filtered)))

    # Steps 3-5:  NA Genus handling logic and taxonomic agglomeration
    tax_df <- as.data.frame(tax_table(physeq_filtered))
    na_genus_indices <- which(is.na(tax_df$Genus))
    if (length(na_genus_indices) > 0) {
      message(paste("Found", length(na_genus_indices), "ASVs with NA Genus. Creating unique ASV-level labels..."))
      original_asv_names <- rownames(tax_df)[na_genus_indices]
      new_na_names <- paste0("Genus_NA_Rep_", original_asv_names)
      tax_df$Genus[na_genus_indices] <- new_na_names
      tax_table(physeq_filtered) <- tax_table(as.matrix(tax_df))
    } else {
      message("No ASVs with NA Genus found.")
    }
    physeq_glommed <- tax_glom(physeq_filtered, taxrank = "Genus", NArm = FALSE)
    final_tax_table <- data.frame(tax_table(physeq_glommed))
    final_taxa_names <- make.unique(as.character(final_tax_table$Genus))
    taxa_names(physeq_glommed) <- final_taxa_names

    # --- 6. Prepare Data for Machine Learning  ---
    otu_tab_obj <- otu_table(physeq_glommed)
    Features_filtered <- if (otu_tab_obj@taxa_are_rows) data.frame(t(otu_tab_obj)) else data.frame(otu_tab_obj)
    
    total_abundance <- colSums(Features_filtered)
    n_taxa_to_keep <- min(8000, ncol(Features_filtered))
    top_taxa_names <- names(sort(total_abundance, decreasing = TRUE))[1:n_taxa_to_keep]
    Features_filtered_red <- Features_filtered[, top_taxa_names, drop = FALSE]
    
    meta_filtered <- data.frame(sample_data(physeq_filtered))
    response_var <- factor(meta_filtered$TimePeriod, levels = c("early", "middle", "late"))
    block_number <- as.factor(meta_filtered$Block)
    
    ml_data_full <- data.frame(Features_filtered_red, response_var, block_number)
    
    train_indices <- which(as.character(ml_data_full$block_number) %in% train_blocks)
    test_indices <- which(as.character(ml_data_full$block_number) %in% test_blocks)
    
    if (length(train_indices) == 0 || length(test_indices) == 0) {
      stop("Data splitting resulted in an empty train or test set. Check block numbers.")
    }
    
    train_data <- ml_data_full[train_indices, ]
    test_data <- ml_data_full[test_indices, ]
    train_data$block_number <- NULL
    test_data$block_number <- NULL

   
    # Create the output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Define descriptive output filenames
    train_path <- file.path(output_dir, "growth_stage_train_data.csv")
    test_path <- file.path(output_dir, "growth_stage_test_data.csv")
    
    message("Saving data for Python...")
    write.csv(train_data, train_path, row.names = FALSE)
    write.csv(test_data, test_path, row.names = FALSE)
    message(paste("  - Saved training data to:", train_path))
    message(paste("  - Saved test data to:", test_path))
    # =============================

    return(list(train_data = train_data, test_data = test_data, full_data = ml_data_full))
}


# ===================================================================
# --- SCRIPT EXECUTION ORCHESTRATOR ---
# 

# --- 1. Argument Parsing ---
# This script is designed to be called with five command-line arguments.
# Lists of items (blocks, fertilizers) should be passed as a single,
# comma-separated string, e.g., "1,3".
#
# USAGE:
# Rscript <script_name>.R <rds_path> <output_dir> <train_blocks> <test_blocks> <fertilizers>
#
# EXAMPLE:
# Rscript src/r/01...R data/Phyl_rare_genus.rds outputs/growth_stage_cv13-24 "1,3" "2,4" "M1P1,N1K1,N1P2K2"
#
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("FATAL: This script requires exactly 5 arguments: <rds_path> <output_dir> <train_blocks> <test_blocks> <fertilizers>", call. = FALSE)
}

# Assign arguments to meaningful variables
input_rds_path <- args[1]
output_dir_path <- args[2]
train_blocks_str <- args[3]
test_blocks_str <- args[4]
fertilizers_str <- args[5]

message("--- Pipeline settings received ---")
message("Input RDS file: ", input_rds_path)
message("Output directory: ", output_dir_path)
message("Training blocks (raw): ", train_blocks_str)
message("Testing blocks (raw): ", test_blocks_str)
message("Fertilizer types (raw): ", fertilizers_str)

# --- 2. Process and Validate Arguments ---
if (!file.exists(input_rds_path)) {
  stop("FATAL: Input RDS file not found at '", input_rds_path, "'")
}

# Parse the comma-separated strings into vectors
train_blocks_vec <- strsplit(train_blocks_str, ",")[[1]]
test_blocks_vec <- strsplit(test_blocks_str, ",")[[1]]
fertilizers_vec <- strsplit(fertilizers_str, ",")[[1]]

# --- 3. Call the Main Function ---
# Use a tryCatch block for robust error handling
tryCatch({
  result_list <- prepare_dataset_phyloseq_from_file(
    rds_filepath = input_rds_path,
    fertilizer_types = fertilizers_vec,
    train_blocks = train_blocks_vec,
    test_blocks = test_blocks_vec,
    output_dir = output_dir_path
  )
}, error = function(e) {
  stop("An error occurred during data preparation: ", e$message)
})

message("\nData preparation complete.")