# Timing analysis for Tensor-cell2cell

## Environment setup

[See instructions for creating the conda environment](./env_setup/README.md)

## Testing CellChat

Add the test_cellchat.q to the queue:

- ```qsub qsub/test_cellchat.q```

## Timing the tools

Run independently:

- ```qsub qsub/tensor_cell2cell.q```
- ```qsub qsub/cellchat_context.q```
- ```qsub qsub/cellchat_sample.q```
