#!/usr/bin/env bash

# Root Directory
ROOT="/home/earmingo/CCC-Benchmark/"

# Output Directory
OUTPUT=$ROOT$"/outputs/"

# Environment
source ~/.bashrc
conda activate cellchat

# CellChat
Rscript $ROOT$"/balf_samples/test/CellChat-Tutorial.r" \
 $ROOT \
 $OUTPUT &>> $OUTPUT$"/test_cellchat.out"

# Stop environment
conda deactivate