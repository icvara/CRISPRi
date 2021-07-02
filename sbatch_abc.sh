#!/bin/bash

#SBATCH --time 00-00:15:00
#SBATCH --node 1
#SBATCH --ntask 1
#SBATCH --cpus-per-task 12
#SBATCH --array 1-8


env | grep SLURM

module load gcc python

DIR=/work/CTR/CI/DCSR/rfabbret/default/eorliac/Support/Icvara/

VENV=$DIR/ABC
source $VENV/bin/activate


IN=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/input.txt)
echo IN = $IN

cd $DIR/CRISPRi
python abc_smc.py $IN $SLURM_CPUS_PER_TASK

deactivate
