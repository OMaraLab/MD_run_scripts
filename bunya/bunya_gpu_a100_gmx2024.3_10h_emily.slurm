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
#
## GROMPP if there's no TPR file (eg, this is the first submission)
if [ ! -f ${SLURM_JOB_NAME}.tpr ]; then
    $GMX grompp -f ${SLURM_JOB_NAME}.mdp -c ${SLURM_JOB_NAME}_start.gro -o ${SLURM_JOB_NAME}.tpr -p ${SLURM_JOB_NAME}.top  -n ${SLURM_JOB_NAME}.ndx -maxwarn 2 &> ${SLURM_JOB_NAME}_grompp_${SLURM_JOB_ID}.txt
fi
#####  RUN MD FOR 9.95 HOURS  #####

$GMXMDRUN -v -deffnm ${SLURM_JOB_NAME} -cpi ${SLURM_JOB_NAME}.cpt || errexit

#####  CHECK IF JOB IS DONE; IF NOT DONE RESUBMIT THIS SCRIPT  #####

# Check the log file for the number of steps completed
steps_done=$(awk '/Statistics over/ {print $3}' ${PBS_JOBNAME}.log | tail -1)
# Check the mdp file for the number of steps we want
steps_wanted=$(awk '/nsteps *=/ {print $3}' ${PBS_JOBNAME}.mdp | tail -1)
# Resubmit if we need to
if (( steps_done < steps_wanted )); then
	echo "Job ${PBS_JOBID}: ${steps_done}/${steps_wanted} steps completed."
	echo "Resubmitting job sequence..."
	qsub $PBS_JOBNAME
else
	echo "Simulation Complete"
fi

#####  END  #####
