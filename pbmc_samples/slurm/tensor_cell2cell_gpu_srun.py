#!/usr/bin/env python
# coding: utf-8
import subprocess

## INPUTS
# Sample cases : [3, 6, 12, 24, 36, 48, 60]
#contexts = [False, True]
#samples = [3, 6, 12, 24, 36, 48, 60]
contexts = [False]
samples = [3]
benchmark_folder = '/home/earmingo/CCC-Benchmark/'

count = 1
for context in contexts:
    if context:
        suffix = 'agg'
        context_input = 'TRUE'
    else:
        suffix = 'ind'
        context_input = 'FALSE'
    for sample in samples:
        srun_command = 'sbatch --job-name=Tensor_cell2cell-{} '.format(count) \
                       + '--output={}outputs/gpu_tensor_cell2cell-{}-{}.out '.format(benchmark_folder, sample, suffix) \
                       + '--error={}outputs/gpu_tensor_cell2cell-{}-{}.err '.format(benchmark_folder, sample, suffix) \
                       + '-p gpu --gres=gpu:1 --exclusive ' \
                       + "--export=folder='{}',samples='{}',context='{}' ".format(benchmark_folder, sample, context_input) \
                       + benchmark_folder + '/pbmc_samples/slurm/tensor_cell2cell.sub'

        exit_status = subprocess.call(srun_command, shell=True)
        if exit_status == 1:
            print("Job {0} failed to submit".format(srun_command))
        count +=1