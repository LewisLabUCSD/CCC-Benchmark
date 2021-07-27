#!/usr/bin/env python
# coding: utf-8
import subprocess

## INPUTS
# Sample cases : [3, 6, 12, 24, 36, 48, 60]
contexts = [False, True]
samples = [3, 6, 12, 24, 36, 48, 60]
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
        qsub_command = 'qsub -N Tensor_cell2cell-{} '.format(count) \
                       + '-o {}outputs/tensor_cell2cell-{}-{}.out '.format(benchmark_folder, sample, suffix) \
                       + '-e {}outputs/tensor_cell2cell-{}-{}.err '.format(benchmark_folder, sample, suffix) \
                       + "-v folder='{}',samples='{}',context='{}' ".format(benchmark_folder, sample, context_input) \
                       + benchmark_folder + '/pbmc_samples/qsub/tensor_cell2cell.q'

        exit_status = subprocess.call(qsub_command, shell=True)
        if exit_status is 1:
            print("Job {0} failed to submit".format(qsub_command))
        count +=1