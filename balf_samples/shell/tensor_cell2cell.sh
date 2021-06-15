#!/usr/bin/env bash

# Root Directory
ROOT="/home/earmingo/CCC-Benchmark/balf_samples/"

# Output Directory
OUTPUT=$ROOT$"/outputs/"

# Environment
source ~/.bashrc
conda activate cellchat

# Tensor-cell2cell
python $ROOT$"/timing_src/time_tensor_cell2cell.py" \
 $ROOT \
 $OUTPUT &>> $OUTPUT$"/tensor_cell2cell.out"

# Stop environment
conda deactivate