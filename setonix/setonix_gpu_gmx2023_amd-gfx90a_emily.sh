#!/bin/bash --login

#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --sockets-per-node=1
#SBATCH --time=10:00:00
#SBATCH --account=pawsey0420-gpu

# === INFO ===
## Built for production runs using Setonix
## Runs jobs in 10 hour blocks with checkpointing and resubmission
## All initial files should have same name, e.g.: .top, .mdp, _start.gro

echo $SLURM_JOB_NODELIST

# Load libraries and the GROMACS module
module load gromacs-amd-gfx90a/2023
module list

GMX='srun -l -u -c 8 gmx_mpi'
GMXMDRUN=$GMX' mdrun -nb gpu -bonded gpu -pin on -update gpu -ntomp 8 -maxh 9.95'


# === ERROR FUNCTION ===
# Kudos to Ada

errexit ()
{
    errstat=$?
    if [ $errstat != 0 ]; then
        # A brief nap so slurm kills us in normal termination
        sleep 5
        echo "Job returned error status $errstat - stopping job sequence $SLURM_JOB_NAME at job $SLURM_JOB_ID"
        exit $errstat
    fi
}


# === GROMPP AND MDRUN ===

# GROMPP if there's no TPR file (e.g., this is the first submission)

if [ ! -f ${SLURM_JOB_NAME}.tpr ]; then
    $GMX grompp -f ${SLURM_JOB_NAME}.mdp -c ${SLURM_JOB_NAME}_start.gro -o ${SLURM_JOB_NAME}.tpr -p ${SLURM_JOB_NAME}.top -maxwarn 2 &> ${SLURM_JOB_NAME}_grompp_${SLURM_JOB_ID}.txt
fi

# Run MD for 9.95 hours

$GMXMDRUN -v -deffnm ${SLURM_JOB_NAME} -cpi ${SLURM_JOB_NAME}.cpt || errexit


# === CHECKPOINT + RESUBMISSION ===
# Check the log file for the number of steps completed
steps_done=`perl -n -e'/Statistics over (\d+) steps using (\d+) frames/ && print $1' ${SLURM_JOB_NAME}.log`

# Check the mdp file for the number of steps we want
steps_wanted=`perl -n -e'/nsteps\s*=\s*(\d+)/ && print $1' ${SLURM_JOB_NAME}.mdp`

# Resubmit if we need to
if (( steps_done < steps_wanted )); then
    echo "Job ${SLURM_JOB_NAME} terminated with ${SLURM_JOB_NAME}/${SLURM_JOB_NAME} steps finished."
    echo "Submitting next job in sequence ${SLURM_JOB_NAME}."
    sbatch ${SLURM_JOB_NAME}
fi
