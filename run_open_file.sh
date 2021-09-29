#!/bin/bash/

set -e
set -x
PWD=$(pwd)
export PATH="ems/elsc-labs/segev-i/moria.fridman/anaconda3/bin:"$PATH
export PYTHONDIR="ems/elsc-labs/segev-i/moria.fridman/anaconda3/bin/"
python3 open_file.py
