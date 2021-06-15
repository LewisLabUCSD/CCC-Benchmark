#!/usr/bin/env bash

# Root Directory
ROOT="/home/earmingo/CCC-Benchmark/balf_samples/"

# Output Directory
OUTPUT=$ROOT$"/outputs/"

# Environment
source ~/.bashrc
conda activate cellchat

# CellChat
Rscript $ROOT$"/timing_src/time_cellchat_bycontext.r" \
 $ROOT \
 $OUTPUT &>> $OUTPUT$"/cellchat_context.out"

# Stop environment
conda deactivate