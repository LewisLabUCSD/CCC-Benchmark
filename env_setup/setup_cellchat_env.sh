#!/bin/bash
# take in environment name: taken from https://www.shellscript.sh/tips/getopt/
NAME="cellchat"
usage()
{
  echo "Usage: bash -i setup_cellchat_env.sh [ -n | --name ENV_NAME ]"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -o n: --long name: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -n | --name) NAME="$2" ; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

# initialize and check conda environment
# echo "NAME   : $NAME"
conda create -y -n "$NAME" python=3.8.8
conda activate "$NAME"
conda info|egrep "conda version|active env"

ACT_ENV="$(conda info|egrep "active environment")"
ACT_ENV=(${ACT_ENV// : / })
ACT_ENV=${ACT_ENV[2]}  
if [[ "$ACT_ENV" != "$NAME" ]]; then
  echo "The environment $NAME has not been activated"
  usage
fi

# begin package installs
conda install -y -c conda-forge r-base r-devtools r-cairo r-systemfonts r-irkernel r-igraph umap-learn
Rscript setup_cellchat_env.r
conda install -y -c conda-forge r-rhpcblasctl=0.20_137
pip install 'cell2cell==0.4.4'
conda deactivate
echo "Complete, activate environment using: conda activate $NAME"
