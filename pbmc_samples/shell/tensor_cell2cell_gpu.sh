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
conda activate cell2cell_gpu

# Tensor-cell2cell
python $ROOT$"/pbmc_samples/timing_src/time_tensor_cell2cell_gpu.py" \
 $ROOT \
 $OUTPUT \
 $NSAMPLES \
 $CONTEXT

##  /home/earmingo/.conda/envs/cell2cell_gpu/bin/python
# &>> $OUTPUT$"/tensor_cell2cell.out"

# Stop environment
conda deactivate