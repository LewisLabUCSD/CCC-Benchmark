#!/usr/bin/env bash

# Root Directory
ROOT=$1

# Output Directory
OUTPUT=$ROOT$"/outputs/"

# ARGS
NSAMPLES=$2
CONTEXT=$3
SEED="8"

# Environment
source ~/.bashrc
conda activate cellchat

# CellChat
Rscript $ROOT$"/pbmc_samples/timing_src/time_cellchat.r" \
 $"--number="$NSAMPLES \
 $"--group="$CONTEXT \
 $"--root="$ROOT \
 $"--outputs="$OUTPUT \
 $"--seed="$SEED

# &>> $OUTPUT$"/cellchat.out"

# Stop environment
conda deactivate