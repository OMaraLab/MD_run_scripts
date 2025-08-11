#!/bin/bash
#SBATCH --account=pawsey0420
#SBATCH --ntasks=512 # This will vary based on your system size. 512 is based on a ~500,000 particle system
#SBATCH --ntasks-per-node=64 # ntasks-per-node x cpus-per-task = 128
#SBATCH --cpus-per-task=2 # must match OMP_NUM_THREADS - 2 OMP threads ideal for UA/AA systems based on 2024 benchmarking.
#SBATCH --exclusive # Only necessary if you a
#SBATCH --time=10:00:00
#SBATCH -m block:block:block # At resource level - place 64 tasks and respective threads next to each other

# Submission script based on the Pawsey Uptake Project conducted in 2024 by the O'Mara Group in conjunction with Dr Deva Deeptimahanti and Dr Emily Kahl (Pawsey Specialists)

echo $SLURM_JOB_NODELIST  # print node id to the end of output file for benchmarking / identifying troublesome nodes
echo "CLUSTERID=SETONIX-CPU" # print which cluster was used, so if a system is run on multiple systems I can separate the log files for benchmarking

# Load modules and switch to 2024 Pawsey environment to access custom GROMACS compile
module unload gcc-native/14.2
module switch pawseyenv/2025.08 pawseyenv/2024.05
module load gcc/12.2.0
module load libfabric/1.15.2.0
module load gromacs/2023.3-mx33m5h
module list

export OMP_NUM_THREADS=2 # explicitly specify for GROMACS

GMX="gmx_mpi"
GMXMDRUN="srun -c $OMP_NUM_THREADS -m block:block:block GMX mdrun -maxh 9.95" # for hybrid -c is recommended for job placement

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
