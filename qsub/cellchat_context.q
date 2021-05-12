#!/bin/bash
#PBS -N CellChat-Contexts
#PBS -m abe
#PBS -M earmingol@eng.ucsd.edu
#PBS -l pmem=32gb
#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o cellchat_context.out

START_TIME=$(date +%s)

echo "Starting the analysis"
/home/earmingo/CCC-Benchmark/shell/cellchat_context.sh

END_TIME=$(date +%s)
RUN_TIME=$(($END_TIME-$START_TIME))
echo $"Analysis completed in $RUN_TIME sec"