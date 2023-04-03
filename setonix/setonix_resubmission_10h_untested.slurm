#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --exclusive
#SBATCH --time=10:00:00
#SBATCH --account=pawsey0420




## General-purpose resubmit script for GROMACS jobs on Setonix

## this runs jobs in ten hour blocks with checkpinting and resubmission, 
## edit $GMXMDRUN to vary the runtime
## nsteps is set in the mdp
## Starting structure, .top, .ndx and .mdp should have the same name as the
## script, and all be in the same folder. Output will also have 
## the same name.
## Eg:  if script name is GlyT2_POPC_CHOL_r1, mdp is GlyT2_POPC_CHOL_r1.mdp
##
## IMPORTANT: SCRIPT NAME SHOULD NOT END WITH .sh


# Define error function so we can see the error code given when something
# important crashes
errexit ()
{
    errstat=$?
    if [ $errstat != 0 ]; then
        # A brief nap so slurm kills us in normal termination
        # Prefer to be killed by slurm if slurm detected some resource excess
        sleep 5
        echo "Job returned error status $errstat - stopping job sequence $SLURM_JOB_NAME at job $SLURM_JOB_ID"
        exit $errstat
    fi
}


module load gromacs/2021.4

GMX='srun -N 1 -n 128 gmx_mpi_d'
GMXMDRUN='srun -N 1 -n 128 gmx_mpi_d mdrun -maxh 9.95'

## GROMPP if there's no TPR file (eg, this is the first submission)

if [ ! -f ${SLURM_JOB_NAME}.tpr ]; then
    $GMX grompp -f ${SLURM_JOB_NAME}.mdp -c ${SLURM_JOB_NAME}_start.gro -o ${SLURM_JOB_NAME}.tpr -p ${SLURM_JOB_NAME}.top  -n ${SLURM_JOB_NAME}.ndx -maxwarn 1 &> ${SLURM_JOB_NAME}_grompp_${SLURM_JOB_ID}.txt
fi

export OMP_NUM_THREADS=1

# run MD

$GMXMDRUN -v -deffnm ${SLURM_JOB_NAME} -cpi ${SLURM_JOB_NAME}.cpt || errexit


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


## clean up linebreaks in output file

sed -i 's|\r|\n|g' slurm-${SLURM_JOB_ID}.out

## backup logfile for debugging

cp ${SURM_JOB_NAME}.log lastlog.txt

#####  END  #####
