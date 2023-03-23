#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --account=pawsey0420


module load gromacs/2021.4
export OMP_NUM_THREADS=1
srun -N 1 -n 128 gmx_mpi_d mdrun -s TKTK.tpr