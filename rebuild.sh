#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
printf "Current directory: %s\n" $SCRIPT_DIR
DNALOADER_DIR="$SCRIPT_DIR/dnaloader"
if [[ ":$PYTHONPATH:" != *":$DNALOADER_DIR:"* ]]; then
    export PYTHONPATH=$DNALOADER_DIR:$PYTHONPATH
fi
echo $PYTHONPATH
printf "Creating enviroment\n"
conda env create -f "environment.yml"
conda activate bwenv
python "$SCRIPT_DIR/setup.py" build_ext --inplace