#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
printf "%s\n" $SCRIPT_DIR
DNALOADER_DIR="$SCRIPT_DIR/dnaloader"
printf "%s\n" $DNALOADER_DIR
export PYTHONPATH="$DNALOADER_DIR:$PYTHONPATH"
python "$SCRIPT_DIR/setup.py" build_ext --inplace