#!/bin/bash
#SBATCH --ntasks-per-node=75
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=10:00:00
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --account=a_omara

export OMP_NUM_THREADS=2

# Load the GROMACS module
module load gcc/11.3.0
source /scratch/project/omara_scratch/gromacs/gromacs_2024.3_CPU_MPI/bin/GMXRC

GMX='gmx_mpi'
GMXMDRUN='srun gmx_mpi mdrun -ntomp 2 -pin on -dlb yes -maxh 0.08 -dd 5 5 3 -dds 0.95'
# NOTE: -ntmpi flag is NOT supported by my compilation


# Add production run GROMPP and MDRUN commands for your system here!
