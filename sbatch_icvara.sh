#!/bin/bash

#SBATCH --time 00-00:15:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --array 1-8


env | grep SLURM

module load gcc python

DIR=/users/ibarbier
source $DIR/myfirstvenv/bin/activate


IN=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/input.txt)
echo IN = $IN


python3 abc_smc.py $IN $SLURM_CPUS_PER_TASK

deactivate
