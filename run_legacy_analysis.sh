#!/bin/bash

# =============================================================================
# MASTER SCRIPT FOR THE LEGACY BIOMARKER ANALYSIS
# =============================================================================
# USAGE:
#   bash run_legacy_analysis.sh <input_file> <split_method> <train_blocks> <test_blocks>
#
# EXAMPLES:
#   bash run_legacy_analysis.sh Input_files/phyloseq_genus.rds random "none" "none"
#   bash run_legacy_analysis.sh Input_files/phyloseq_genus.rds block "1,3" "2,4"
# =============================================================================

set -e

# --- 1. Argument Validation ---
if [ "$#" -ne 4 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: bash run_legacy_analysis.sh <input_file> <split_method> <train_blocks> <test_blocks>"
    exit 1
fi

INPUT_FILE=$1
SPLIT_METHOD=$2
TRAIN_BLOCKS=$3
TEST_BLOCKS=$4

# Validate the split method argument
if [[ "$SPLIT_METHOD" != "random" && "$SPLIT_METHOD" != "block" ]]; then
    echo "Error: Invalid split_method. Must be 'random' or 'block'."
    exit 1
fi

# Validate that the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found at '$INPUT_FILE'"
    exit 1
fi

echo "--- Starting legacy analysis pipeline ---"
echo "Input file:   $INPUT_FILE"
echo "Split method: $SPLIT_METHOD"

# --- 2. Define Output Directories DYNAMICALLY ---
if [ "$SPLIT_METHOD" == "block" ]; then
    # Clean commas for folder naming (e.g., "1,3" becomes "13")
    TRAIN_CLEAN=$(echo $TRAIN_BLOCKS | tr -d ',')
    TEST_CLEAN=$(echo $TEST_BLOCKS | tr -d ',')
    OUTPUT_ROOT="Outputs/legacy/block_train${TRAIN_CLEAN}_test${TEST_CLEAN}"
else
    OUTPUT_ROOT="Outputs/legacy/random_split"
fi

DATASET_DIR="$OUTPUT_ROOT/01_datasets"
MODEL_DIR="$OUTPUT_ROOT/02_models_and_shap"
RANKING_DIR="$OUTPUT_ROOT/03_biomarker_rankings"

echo "[INFO] Results will be saved in: $OUTPUT_ROOT"

mkdir -p "$DATASET_DIR"
mkdir -p "$MODEL_DIR"
mkdir -p "$RANKING_DIR"

# --- 3. Run Pipeline Steps ---

echo ""
echo "[INFO] Step 1/3: Generating datasets..."
# We pass all 5 required arguments to the R script here
Rscript src/r/01_GeneratingDataSetsForLegacyAnalysis.R \
    "$INPUT_FILE" \
    "$SPLIT_METHOD" \
    "$DATASET_DIR" \
    "$TRAIN_BLOCKS" \
    "$TEST_BLOCKS"

echo ""
echo "[INFO] Step 2/3: Training models and running SHAP analysis..."
python src/python/02_RandomForestTrainingAndShapleyAnalysisForLegacyData.py \
    "$DATASET_DIR" \
    "$MODEL_DIR"

echo ""
echo "[INFO] Step 3/3: Generating final rankings and plots..."
Rscript src/r/03_Legacy_Biomarker_Ranking.R \
    "$MODEL_DIR" \
    "$RANKING_DIR"

# --- 4. Completion Message ---
echo ""
echo "---"
echo "Pipeline finished successfully!"
echo "All outputs are located in: $OUTPUT_ROOT"