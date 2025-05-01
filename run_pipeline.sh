#!/usr/bin/env bash
# -----------------------------------------------
# Runs the pipeline.py script with given arguments
#
# Usage:
#   ./run_pipeline.sh <input_data_dir> <output_dir>
#
# If arguments are not provided, defaults will be used:
#   input_data_dir = example_data/
#   output_dir = results/
# -----------------------------------------------

set -e  # Exit immediately if any command fails

DATA_DIR=${1:-example_data}
OUT_DIR=${2:-results}

# Optional: uncomment the following if using a conda environment
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate nbase-seic

python pipeline.py --data_dir "$DATA_DIR" --out_dir "$OUT_DIR"