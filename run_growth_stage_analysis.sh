#!/bin/bash

# =============================================================================
# MASTER SCRIPT FOR A SINGLE GROWTH STAGE ANALYSIS FOLD
# =============================================================================
# This script runs the full pipeline for a specific cross-validation fold:
#   1. Generates training/test datasets using R.
#   2. Trains models and runs SHAP analysis using Python.
#   3. Creates final rankings and plots using R.
#
# USAGE:
#   bash run_growth_stage_analysis.sh <input_file> <train_blocks> <test_blocks> <fertilizers>
#
# EXAMPLE for the 1,3 vs 2,4 fold:
#   bash run_growth_stage_analysis.sh Input_files/Phyloseq_genus.rds "1,3" "2,4" "M1P1,N1K1,N1P2K2"
#
# =============================================================================

set -e

# --- 1. Argument Validation ---
if [ "$#" -ne 4 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: bash run_growth_stage_analysis.sh <input_file> <train_blocks> <test_blocks> <fertilizers>"
    exit 1
fi

INPUT_FILE=$1
TRAIN_BLOCKS=$2
TEST_BLOCKS=$3
FERTILIZERS=$4

echo "--- Starting Full Pipeline for Growth Stage Analysis Fold ---"
echo "Training on blocks: $TRAIN_BLOCKS"
echo "Testing on blocks:  $TEST_BLOCKS"

# --- 2. Construct Dynamic Output Directory Names ---
# We strip commas to create a clean suffix for the folder name
TRAIN_CLEAN=$(echo $TRAIN_BLOCKS | tr -d ',')
TEST_CLEAN=$(echo $TEST_BLOCKS | tr -d ',')

# Consistent naming convention: Outputs/growth_stage/trainXX_testXX
OUTPUT_ROOT="Outputs/growth_stage/train${TRAIN_CLEAN}_test${TEST_CLEAN}"

# Define the subdirectories for each step of the pipeline
DATASET_DIR="$OUTPUT_ROOT/01_datasets"
MODEL_DIR="$OUTPUT_ROOT/02_models_and_shap"
RANKING_DIR="$OUTPUT_ROOT/03_final_rankings"

echo "[INFO] All outputs for this fold will be saved in: $OUTPUT_ROOT"

# Create all necessary directories
mkdir -p "$DATASET_DIR"
mkdir -p "$MODEL_DIR"
mkdir -p "$RANKING_DIR"


# =============================================================================
# --- 3. RUN PIPELINE STEPS ---
# =============================================================================

# --- Step 1: Generate Datasets using R ---
echo ""
echo "[INFO] Step 1/3: Generating training and test datasets..."
Rscript src/r/01_generate_growth_stage_datasets.R \
    "$INPUT_FILE" \
    "$DATASET_DIR" \
    "$TRAIN_BLOCKS" \
    "$TEST_BLOCKS" \
    "$FERTILIZERS"

# --- Step 2: Train Models and Analyze using Python ---
echo ""
echo "[INFO] Step 2/3: Training models and running SHAP analysis..."
python src/python/02_train_growth_stage_models.py \
    "$DATASET_DIR" \
    "$MODEL_DIR"

# --- Step 3: Create Final Rankings and Plots using R ---
echo ""
echo "[INFO] Step 3/3: Creating final rankings and plots..."
Rscript src/r/03_create_growth_stage_rankings.R \
    "$MODEL_DIR" \
    "$RANKING_DIR"

echo ""
echo "--- Full pipeline for this fold has completed successfully! ---"
echo "Results are located in: $OUTPUT_ROOT"