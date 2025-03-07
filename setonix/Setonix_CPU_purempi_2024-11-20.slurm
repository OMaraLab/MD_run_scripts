#!/bin/bash
#SBATCH --account=pawsey0420
#SBATCH --ntasks=512 # Varies based on system size, 512 for ~500,000 particle system.
#SBATCH --ntasks-per-node=64 # 128, 32, 16 etc. 64 ideal here - doubles the memory per core (recommended to leave --exclusive off to allow the rest of the node to be utilised with jobs - the memory per core will decrease as the node fills up reducing overall performance.
##SBATCH --exclusive ## Only use if you are not utilising a full node (i.e 64 ntasks per node) and NEED the speed.
#SBATCH --time=10:00:00

echo $SLURM_JOB_NODELIST  # print node id to the end of output file for benchmarking / identifying troublesome nodes
echo "CLUSTERID=SETONIX-CPU" # print which cluster was used, so if a system is run on multiple systems I can separate the log files for benchmarking

# Load modules
module load gromacs/2023.3-mx33m5h
module list

export OMP_NUM_THREADS=1 # explicitly specify for GROMACS

GMX="gmx_mpi"
GMXMDRUN="srun $GMX mdrun -maxh 9.95"

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


#####  START MD WITH A GROMPP  #####

## GROMPP if there's no TPR file (eg, this is the first submission)
if [ ! -f ${SLURM_JOB_NAME}.tpr ]; then
    $GMX grompp -f ${SLURM_JOB_NAME}.mdp -c ${SLURM_JOB_NAME}_start.gro -o ${SLURM_JOB_NAME}.tpr -p ${SLURM_JOB_NAME}.top  -n ${SLURM_JOB_NAME}.ndx -maxwarn 2 &> ${SLURM_JOB_NAME}_grompp_${SLURM_JOB_ID}.txt
fi

#####  RUN MD FOR 9.95 HOURS  #####

$GMXMDRUN -v -deffnm ${SLURM_JOB_NAME} -cpi ${SLURM_JOB_NAME}.cpt || errexit

#####  CHECK IF JOB IS DONE; IF NOT DONE RESUBMIT THIS SCRIPT  #####


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

# ## backup logfile for live debugging

# cp ${SLURM_JOB_NAME}.log lastlog.txt

### append with some info onthe run for debugging and benchmarking


#####  END  #####
