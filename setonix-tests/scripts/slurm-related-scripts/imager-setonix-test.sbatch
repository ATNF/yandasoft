#!/bin/bash -l
#SBATCH --partition=workq
#SBATCH --clusters=galaxy
#SBATCH --account=askap
# No further constraints applied
# No reservation requested
# No email notifications sent
#SBATCH --time=24:00:00
#SBATCH --nodes=10
#SBATCH --job-name=contIMG
#SBATCH --export=NONE
#SBATCH --output="imaging-test-setonix-%j.out"

# Need to load the slurm module directly
module load slurm
# Using user-defined askapsoft module
module use /group/askap/modulefiles
module load askapdata
module unload askapsoft
module load askapsoft

export OMP_NUM_THREADS=20
export MPICH_RANK_REORDER_METHOD=0

# Run the imager, calibrating if not the first time.
log="imager-setonix-test_${SLURM_JOB_ID}.log"
parset="imager-setonix-test.in"
echo "--- Imaging with imager ---" > "$log"
#export OMP_NUM_THREADS=20
NCORES=55
srun --export=ALL --bcast=/tmp/imager --cpus-per-task=20 --ntasks=1 --ntasks-per-node=1 imager  -c "$parset" : --cpus-per-task=1 --ntasks=54 --ntasks-per-node=6 imager  -c "$parset" >> "$log"

