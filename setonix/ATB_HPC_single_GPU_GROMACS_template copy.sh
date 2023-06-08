#!/bin/bash --login
#SBATCH --job-name=<JOB_NAME>
#SBATCH --partition=gpu
#SBATCH --nodes=1              #1 nodes in this example
#SBATCH --ntasks-per-node=1    #1 tasks for the 1 GPUs in this job
#SBATCH --gpus-per-node=1      #1 GPUs in this job
#SBATCH --sockets-per-node=1   #Use the 1 slurm-sockets in this job
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --account=<ACCOUNT>-gpu

module load PrgEnv-cray
module load rocm craype-accel-amd-gfx90a

module load gcc/12.1.0
module load gromacs-amd-gfx90a/2022.3.amd1_174
export GMX_MAXBACKUP=-1

unset OMP_NUM_THREADS

srun -l -u -c 8 gmx mdrun <ARGS> -nb gpu -bonded gpu -pin on -update gpu -ntomp 8 -ntmpi 1
