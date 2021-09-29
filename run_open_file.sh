#!/bin/bash/

set -e
set -x
PWD=$(pwd)
export PATH="/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/bin:"$PATH
export PYTHONDIR="/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/bin/"
python3 main.py
