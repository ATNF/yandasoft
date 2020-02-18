#!/bin/bash -l
#SBATCH --partition=workq
#SBATCH --clusters=galaxy
# Using the default account
# No reservation requested
#SBATCH --time=00:10:00
#SBATCH --ntasks=7
#SBATCH --ntasks-per-node=8
#SBATCH --job-name=test00
#SBATCH --mail-user=stephen.ord@csiro.au
#SBATCH --mail-type=ALL
#SBATCH --export=ALL
#SBATCH --output=/group/askap/sord/yandasoft/tests/data/simulation/synthregression/galaxy.out
#SBATCH --chdir=/group/askap/sord/yandasoft/tests/data/simulation/synthregression

export PATH=${PATH}:/group/askap/sord/yanda/bin
source /group/askap/sord/yandasoft/deploy/pawsey/ModuleFile
module unload askapdata
module load askapdata
export MEASURESDATA=${AIPSPATH}/data
srun -n 1 /group/askap/sord/yanda/bin/csimulator -c wtermtest_template.in
srun -n 2 /group/askap/sord/yanda/bin/imager -c sbatch_test.in

