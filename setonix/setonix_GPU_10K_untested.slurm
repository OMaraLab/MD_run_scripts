"GPU01_TEST10K" 23L, 661B                                                                                            3,33          All
#!/bin/bash --login

#SBATCH --account=pawsey-0420-gpu
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --time=1:00:00

# test script for running on gpus on setonix
# currently fails with error message:
#
#   sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified
#
# so I guess we need to wait for a gpu account string to be switched on


gromacs-amd-gfx90a/2022.3.amd1_174

export OMP_NUM_THREADS=1

GMX="srun -c 8 gmx"
GMXMDRUN="srun -c 8 gmx mdrun -nsteps 10000"

if [ ! -f ${SLURM_JOB_NAME}.tpr ]; then
    $GMX grompp -f ${SLURM_JOB_NAME}.mdp -c ${SLURM_JOB_NAME}_start.gro -o ${SLURM_JOB_NAME}.tpr -p ${SLURM_JOB_NAME}.top  -n ${SLURM_JOB_NAME}.ndx -maxwarn 1 &> ${SLURM_JOB_NAME}_grompp_${SLURM_JOB_ID}.txt
fi

$GMXMDRUN -v -deffnm ${SLURM_JOB_NAME} -cpi ${SLURM_JOB_NAME}.cpt || errexit
