Running the bash script will set up a conda environment with all necessary packages to run the CellChat package (only tested on Linux)

Use the following command:

```
$ bash -i setup_cellchat_env.sh -n <env_name> # if no name is provided, environment name will default to "cellchat"
$ conda activate <env_name>
```

Once in the activated environment, timing scripts should be able to run