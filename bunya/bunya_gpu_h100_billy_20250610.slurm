#!/bin/bash -l

#SBATCH --job-name=meta-test                # Job name
#SBATCH --qos=sxm                           # QoS name
#SBATCH --partition=gpu_sxm                 # Partition name
#SBATCH --ntasks-per-node=2               # Tasks per node
#SBATCH --nodes=1
#SBATCH --mem=250G                          # Memory per node
#SBATCH --gres=gpu:h100:1                 # Request 1 H100 GPU
#SBATCH --time=00:06:00                     # Walltime limit
#SBATCH --account=a_omara                   # Account name
#SBATCH --output=NavEquil.out           # Save standard output
#SBATCH --error=NavEquil.err            # Save error messages

# This is an example script by Dr Billy Williams-Noonan.

# Load libraries used during compilation

module load cmake/3.26.3-gcccore-12.3.0
module load foss/2023a
module load cuda/12.2.0
module load python/3.11.3-gcccore-12.3.0

export MPIEXEC=srun
 
# Load Plumed Module

source /sw/local/rocky8/xeonsp4_h100/rcc/software/Plumed/plumed-2.9.2/sourceme.sh
 
# Load the GROMACS module

source /sw/local/rocky8/xeonsp4_h100/rcc/software/Gromacs/gromacs-2023.5/bin/GMXRC
 
gmx_mpi grompp -f minim.mdp -o em1.tpr -c ions.gro -p Nav_1.8.top -maxwarn 3 > grompp_em1.out
srun -n 1 gmx_mpi mdrun -v -deffnm em1 -ntomp 1 -maxh 1.95 -nb gpu
