#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --partition=gpu_cuda
#SBATCH --qos=gpu
#SBATCH --gres=gpu:nvidia_a100_80gb_pcie_1g.10gb:2
#SBATCH --account=a_omara

# Requesting the use of 2 nodes on a GPU node with --gres=gpu:[type]:[number] increases performance, but increases wait time in queue. Users choice.

# Load libraries and the GROMACS module
# Kudos to Angus
module load gcc/11.3.0
module load cmake/3.24.3-gcccore-11.3.0
module load python/3.11.3-gcccore-12.3.0

source /scratch/project/omara_scratch/gromacs/gromacs_2024.3_a100/bin/GMXRC


GMX='srun gmx'
GMXMDRUN=$GMX' mdrun -ntomp 2 -ntmpi 32 -nb gpu -bonded gpu -pin on -maxh 9.95'


# === Error Function ===
# Kudos to Ada

errexit () {
    errstat=$?
    if [ $errstat != 0 ]; then
        # A brief nap so slurm kills us in normal termination
        sleep 5
        echo "Job returned error status $errstat - stopping job sequence $SLURM_JOB_NAME at job $SLURM_JOB_ID"
        exit $errstat
    fi
}
echo "Error function setup working"


# === GROMPP and MDRUN ===
# Add your GROMPP and MDRUN commands here!
