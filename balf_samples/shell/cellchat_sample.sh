#!/usr/bin/env bash

# Root Directory
ROOT="/home/earmingo/CCC-Benchmark/"

# Output Directory
OUTPUT=$ROOT$"/outputs/"

# Environment
source ~/.bashrc
conda activate cellchat

# CellChat
Rscript $ROOT$"/balf_samples/timing_src/time_cellchat_bysample.r" \
 $ROOT \
 $OUTPUT &>> $OUTPUT$"/cellchat_sample.out"

# Stop environment
conda deactivate