#!/bin/bash
#SBATCH --ntasks-per-node=96
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --account=a_omara

# Load the GROMACS module
module load gcc/11.3.0
source /scratch/project/omara_scratch/gromacs/gromacs_2024.3_CPU_MPI/bin/GMXRC

GMX='gmx_mpi'
GMXMDRUN='srun gmx_mpi mdrun -pin on -dlb yes -maxh 9.95'
# NOTE: -ntmpi flag is NOT supported by my compilation, but -ntomp flag is


# Add production run GROMPP and MDRUN commands for your system here!
