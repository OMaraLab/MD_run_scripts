#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --account=pawsey0420

PROD_TPR="/path/to/system.tpr"

module load gromacs/2021.4 # note that setonix gromacs was double precision as of benchmarking 2022
export OMP_NUM_THREADS=1
srun -N 1 -n 128 gmx_mpi_d mdrun -s ${PROD_TPR}