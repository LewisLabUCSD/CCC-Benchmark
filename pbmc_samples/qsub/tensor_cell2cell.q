#!/bin/bash
#PBS -m abe
#PBS -M earmingol@eng.ucsd.edu
#PBS -l pmem=128gb
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe

# Inputs
ROOT=${folder} #"/home/earmingo/CCC-Benchmark/"
NSAMPLES=${samples} #"12"
CONTEXT=${context} #"TRUE"

# Script
START_TIME=$(date +%s)

echo "Starting the analysis"
$ROOT$"/pbmc_samples/shell/tensor_cell2cell.sh" $ROOT $NSAMPLES $CONTEXT

END_TIME=$(date +%s)
RUN_TIME=$(($END_TIME-$START_TIME))
echo $"Analysis completed in $RUN_TIME sec"