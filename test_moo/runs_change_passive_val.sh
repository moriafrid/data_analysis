#!/bin/bash

#pass 1 argument = RM
#pass 2 argument = RA
#pass 3 argument = CM
#pass 4 argument = shrinkage_by
#pass 5 argumant = resize_dend_by
#pass 6 argument = passive_vel_name


# Write output as following (%j is JOB_ID)
#SBATCH -o outputs/output-%j.out
#SBATCH -e errors/error-%j.err
#SBATCH --mem 60000
#SBATCH -t 2-0
#SBATCH -c 1


set -x

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS
#
#export PATH="/ems/elsc-labs/segev-i/moria.fridman/anaconda3/bin:$PATH"  ## <<-- Change this to the path that you use
#export PYTHONDIR="/ems/elsc-labs/segev-i/moria.fridman/anaconda3/envs/project/bin/python"

export profile=''  #Note the profile name, you will need to use it in the python script
ipython profile create --parallel --profile=${profile}

##pass ipcluster stop --profile=${profile} &
#pass sleep 50
#pass ipcluster start --profile=${profile} --n=$1 &
#pass sleep 60

echo "Launching job"
conda init
conda activate project
python MOO_get_parameters.py $1 $2 $3 $4 $5 $6 ${profile}
#conda deactivate
