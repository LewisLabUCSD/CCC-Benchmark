#!/usr/bin/env bash

# Root Directory
ROOT=$1

# Output Directory
OUTPUT=$ROOT$"/outputs/"

# ARGS
NSAMPLES=$2
CONTEXT=$3

# Environment
source ~/.bashrc
conda activate cellchat

# Tensor-cell2cell
python $ROOT$"/pbmc_samples/timing_src/time_tensor_cell2cell.py" \
 $ROOT \
 $OUTPUT \
 $NSAMPLES \
 $CONTEXT

# &>> $OUTPUT$"/tensor_cell2cell.out"

# Stop environment
conda deactivate