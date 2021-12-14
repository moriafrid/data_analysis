#!/bin/bash


#pass 1 argument = offspringsize
#pass 2 argument = number of genarations
#pass 3 argument = seed number
#pass 4 argument = cell to optimize
#pass 5 argument = size of ipcluster
#pass 6 argument = Ra_100 || Ra_250
#pass 7 argument = spine type human_spine|| mouse_spine || shaft_spine
#pass 8 argument = Rneck normal_neck || float of Ra for neck

set -e
set -x

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

OFFSPRING_SIZE=$1
MAX_NGEN=$2


export PATH="/ems/elsc-labs/segev-i/yoni.leibner/anaconda2/bin:$PATH"  ## <<-- Change this to the path that you use
export PATH="/ems/elsc-labs/segev-i/yoni.leibner/neuron/nrn/x86_64/bin:$PATH"
export PATH="/ems/elsc-labs/segev-i/yoni.leibner/neuron/nrn/lib/python:$PATH"
export PYTHONDIR="/ems/elsc-labs/segev-i/yoni.leibner/anaconda2/bin/"

export profile=yoni_active_models$3$4  #Note the profile name, you will need to use it in the python script
ipython profile create --parallel --profile=${profile}

# ipcluster stop --profile=${profile} &
sleep 50
ipcluster start --profile=${profile} --n=$5 &
sleep 120

echo "Launching job"
python MOO_fit_syn_response_same_weights_first_cm_g_pas.py $3 ${OFFSPRING_SIZE} ${MAX_NGEN} $4 ${profile} $6 $7 $8
